"""Regenerate nested_summary.csv from the per-target checkpoints.

The summary is derived, not maintained by hand. It was maintained by hand once, and
drifted: ALK, EGFR and KIT were still recorded as "flat only" after their nested runs
had completed, so the file disagreed with nested_<target>.json and with Table S8 of the
manuscript. verify_new.py now re-derives the table and compares, so that cannot recur.

Run from this directory: python3 make_nested_summary.py
"""
import json, glob, os, csv
import numpy as np
import pandas as pd

ROOT = os.environ.get('SAVDR_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ORDER = ['ABL1_ENST00000318560', 'ALK', 'KIT', 'SMO', 'NT5C2', 'EGFR',
         'MET', 'ESR1_ENST00000206249', 'IDH1_ENST00000415913',
         'FLT3', 'BTK_ENST00000308731']

MEASURED = 'measured'
INNER_FAIL = 'all 6 configurations fail on the inner folds'
NO_TRAIN = 'no configuration trains, even on the flat grid'

COLS = ['target', 'n', 'n_pos', 'n_configs_trainable', 'n_configs_untrainable',
        'grid_min', 'grid_median', 'grid_max', 'spread', 'max_minus_median',
        'best_config', 'nested_auc', 'optimism', 'stage']


def load(t):
    """The full checkpoint if the nested run produced one, else the flat-only record.

    Where both exist the full one wins, except for flat_best_config: the earlier flat-only
    run recorded it and the full checkpoint does not always carry it forward.
    """
    full, flat = f'nested_{t}.json', f'nested_{t}_flat.json'
    if os.path.exists(full):
        j, src = json.load(open(full)), full
        if not j.get('flat_best_config') and os.path.exists(flat):
            j = dict(j, flat_best_config=json.load(open(flat)).get('flat_best_config'))
        return j, src
    if os.path.exists(flat):
        return json.load(open(flat)), flat
    raise SystemExit(f'no checkpoint for {t}')


def row(t):
    j, src = load(t)
    fa = {k: v for k, v in (j.get('flat_auc') or {}).items() if v is not None}
    nt = j.get('flat_not_trainable') or []
    # n and n_pos come from the curated variant table, which is what every other
    # verifier counts. The checkpoint's own n_pos + n_neg is smaller: it counts only the
    # variants scored at the intersection of the six grid configurations.
    v = pd.read_csv(f'{ROOT}/inputs/variants/{t}.csv')
    if int(v.y.sum()) != j['n_pos']:
        raise SystemExit(f'{t}: variant table has {int(v.y.sum())} positives, '
                         f'checkpoint says {j["n_pos"]}')
    r = {'target': t, 'n': len(v), 'n_pos': int(v.y.sum()),
         'n_configs_trainable': len(fa), 'n_configs_untrainable': len(nt),
         'best_config': j.get('flat_best_config') or '',
         'nested_auc': '', 'optimism': ''}
    if fa:
        g = list(fa.values())
        mn, md, mx = min(g), float(np.median(g)), max(g)
        r.update(grid_min=round(mn, 4), grid_median=round(md, 4), grid_max=round(mx, 4),
                 spread=round(mx - mn, 4), max_minus_median=round(mx - md, 4))
    else:
        r.update(grid_min='', grid_median='', grid_max='', spread='', max_minus_median='')

    if j.get('nested_auc') is not None:
        r['nested_auc'] = round(j['nested_auc'], 4)
        # The stored optimism is the difference of the unrounded AUCs, which is not always
        # the difference of the two rounded columns: ALK's 0.0010 is 0.9287500 - 0.9277462
        # = 0.0010038, while the rounded columns give 0.0011. Keep the stored value.
        r['optimism'] = round(j['optimism'], 4)
        r['stage'] = MEASURED
    elif j.get('nested_computable') is False or 'declines to train on the inner' in (j.get('note') or ''):
        r['stage'] = INNER_FAIL
    elif not fa:
        r['stage'] = NO_TRAIN
    else:
        raise SystemExit(f'{t}: {src} has neither a nested result nor a stated reason')
    return r


if __name__ == '__main__':
    rows = [row(t) for t in ORDER]
    with open('nested_summary.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, COLS)
        w.writeheader()
        w.writerows(rows)
    m = [r for r in rows if r['stage'] == MEASURED]
    print(f'{len(rows)} targets: {len(m)} measured, '
          f'{sum(r["stage"] == INNER_FAIL for r in rows)} not computable on the inner folds, '
          f'{sum(r["stage"] == NO_TRAIN for r in rows)} with no trainable configuration')
    print('optimism range: '
          f'{min(r["optimism"] for r in m):.4f} to {max(r["optimism"] for r in m):.4f}')
