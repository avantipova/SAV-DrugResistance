"""The one test that could yield a POSITIVE claim: discrimination WITHIN a residue position.

Everything measured so far strips the between-position shortcut (region margins, LOTO) and finds
nothing left. But those tests still pool variants across positions, so a model could in
principle hold real substitution-level signal and have it drowned out.

The strict question: at a position that is known to harbour resistance, does the model rank the
resistant substitution above the non-resistant substitutions AT THAT SAME POSITION? That is the
"which substitution", not "where", question -- exactly what a clinician looking at a novel
variant at a known hotspot needs.

Metric: within-position concordance. For every position carrying at least one positive and at
least one negative, form all (positive, negative) pairs at that position and score
(concordant + 0.5 * tied) / pairs. Pooling over positions gives a position-stratified AUC.

Built-in sanity check: `pos_only` and `nearest_pos` assign the SAME score to every variant at a
given position, so their within-position concordance must come out at exactly 0.500. If it does
not, the metric is implemented wrong. Nothing else in this analysis is trustworthy unless those
two land on 0.500 exactly.
"""
import os as _os

def _savdr_root():
    """Locate the repository root by walking up to the .savdr-root marker."""
    d = _os.path.dirname(_os.path.abspath(__file__))
    while True:
        if _os.path.exists(_os.path.join(d, '.savdr-root')):
            return d
        parent = _os.path.dirname(d)
        if parent == d:
            raise SystemExit('cannot locate repository root: .savdr-root marker not found')
        d = parent

_SAVDR_ROOT = _os.environ.get('SAVDR_ROOT') or _savdr_root()

import glob, re, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # baselines, first_record, ci, delong, AAS/KD/VOL/CHG, U


def nearest_pos_scores(d, seed):
    y = d.y.values; s = np.zeros(len(y))
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(y.reshape(-1, 1), y, d.pos.values):
        hot = np.unique(d.pos.values[tr][y[tr] == 1])
        s[te] = 0.0 if len(hot) == 0 else -np.abs(d.pos.values[te][:, None] - hot[None, :]).min(1)
    return s


def within_position(y, pos, s):
    """Returns (concordance, n_pairs, n_positions_used)."""
    conc = tot = nposit = 0.0
    for p in np.unique(pos):
        m = pos == p
        sp, sn = s[m & (y == 1)], s[m & (y == 0)]
        if len(sp) == 0 or len(sn) == 0: continue
        nposit += 1
        d = sp[:, None] - sn[None, :]
        conc += (d > 0).sum() + 0.5 * (d == 0).sum()
        tot += d.size
    return (conc / tot if tot else np.nan), int(tot), int(nposit)


def boot_within(y, pos, s, B=600, seed=0):
    """Bootstrap over POSITIONS (the unit of analysis here), not over variants."""
    rng = np.random.default_rng(seed)
    usable = [p for p in np.unique(pos)
              if ((pos == p) & (y == 1)).any() and ((pos == p) & (y == 0)).any()]
    if len(usable) < 2: return (np.nan, np.nan)
    out = []
    for _ in range(B):
        pick = rng.choice(usable, len(usable))
        conc = tot = 0.0
        for p in pick:
            m = pos == p
            d = s[m & (y == 1)][:, None] - s[m & (y == 0)][None, :]
            conc += (d > 0).sum() + 0.5 * (d == 0).sum()
            tot += d.size
        if tot: out.append(conc / tot)
    return tuple(np.percentile(out, [2.5, 97.5])) if out else (np.nan, np.nan)


groups = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m: continue
    key = (m.group(1), int(m.group(2)), int(m.group(3)), bool(m.group(4)))
    groups.setdefault(key, []).append((int(m.group(5)), f))

rows = []
for (name, L, MNA, region), runs in sorted(groups.items()):
    if region: continue                   # the full table already contains every variant
    acc = {}
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        sc = baselines(d, name, seed)
        sc['MultiPASS'] = np.nan_to_num(d.multipass_pa.values, nan=0.0)
        sc['nearest_pos'] = nearest_pos_scores(d, seed)
        for k, v in sc.items():
            acc.setdefault(k, []).append(stats.rankdata(v) / len(v))
    y, pos = d.y.values, d.pos.values
    avg = {k: np.mean(v, 0) for k, v in acc.items()}
    for k, v in avg.items():
        c, npairs, nposit = within_position(y, pos, v)
        lo, hi = boot_within(y, pos, v)
        rows.append(dict(target=name, model=k, n_positions=nposit, n_pairs=npairs,
                         within_pos_auc=round(c, 4) if c == c else np.nan,
                         ci_lo=round(lo, 3) if lo == lo else np.nan,
                         ci_hi=round(hi, 3) if hi == hi else np.nan))

df = pd.DataFrame(rows)
df.to_csv(f'{_SAVDR_ROOT}/multipass/within_position.csv', index=False)

print('=' * 96)
print('WITHIN-POSITION DISCRIMINATION (full setting, seed-averaged out-of-fold scores)')
print('=' * 96)
w = df.pivot_table(index='target', columns='model', values='within_pos_auc')
cols = [c for c in ['MultiPASS', 'full_gbm', 'subst_lr', 'pos_only', 'nearest_pos'] if c in w.columns]
npos = df[df.model == 'MultiPASS'].set_index('target')[['n_positions', 'n_pairs']]
print(pd.concat([npos, w[cols]], axis=1).to_string())

print('\nSANITY CHECK -- these two MUST be exactly 0.500 (one score per position):')
for k in ('pos_only', 'nearest_pos'):
    if k not in w.columns: continue
    v = w[k].dropna()
    print(f'  {k:12s} values {sorted(set(v.round(6)))}  '
          + ('OK' if np.allclose(v, 0.5) else '*** METRIC IS WRONG ***'))

print('\nMultiPASS, with a 95% CI bootstrapped over positions:')
mp = df[df.model == 'MultiPASS'].set_index('target')
for t in mp.index:
    r = mp.loc[t]
    if r.within_pos_auc != r.within_pos_auc:
        print(f'  {t:22s} no position carries both a positive and a negative'); continue
    verdict = 'excludes 0.5' if (r.ci_lo > 0.5 or r.ci_hi < 0.5) else 'includes 0.5'
    print(f'  {t:22s} {r.within_pos_auc:.3f}  CI [{r.ci_lo:.3f}, {r.ci_hi:.3f}]  {verdict}'
          f'   ({r.n_positions} positions, {r.n_pairs} pairs)')

tot = df[df.model == 'MultiPASS'].n_pairs.sum()
print(f'\ntotal usable positions across all six targets: '
      f'{int(df[df.model == "MultiPASS"].n_positions.sum())}, {int(tot)} pairs')
print('mean within-position AUC:')
print(w[cols].mean().round(3).to_string())
