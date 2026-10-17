"""Three additions that need no new MultiPASS compute, on the existing fold score matrices.

1. `nearest_pos` -- the harshest possible position baseline. Not a fitted model at all: for a
   held-out variant, score = -(distance in residues to the nearest RESISTANT position in that
   fold's training set). One integer per variant, no wild-type residue, no substituted residue,
   no sequence, no training. If this matches MultiPASS, "the model has learned where resistance
   mutations sit" stops being an inference about pos_only's GBM and becomes a direct reading.

2. Operating points -- precision@top-10/20/50, and sensitivity at 95% and 99% specificity.
   AUC and AUPRC are threshold-free; a tool with a web front end is used at a threshold.

3. Calibration of MultiPASS's Pa -- Pa is presented as a probability of activity, so test it as
   one: Brier score, and a 10-bin reliability curve, on the raw Pa averaged over seeds (NOT the
   rank-averaged score, which is deliberately not a probability).
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
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # helpers: midrank, delong, baselines, ci, first_record


def nearest_pos_scores(d, seed):
    """Out-of-fold score = -(residue distance to nearest training-fold positive position)."""
    y, g = d.y.values, d.pos.values
    s = np.zeros(len(y))
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(y.reshape(-1, 1), y, g):
        hot = np.unique(d.pos.values[tr][y[tr] == 1])
        if len(hot) == 0:
            s[te] = 0.0
            continue
        s[te] = -np.abs(d.pos.values[te][:, None] - hot[None, :]).min(1)
    return s


def op_points(y, s):
    o = np.argsort(-s)
    ys = y[o]
    out = {}
    for k in (10, 20, 50):
        out[f'prec_top{k}'] = round(float(ys[:k].mean()), 3) if len(ys) >= k else np.nan
    neg = s[y == 0]
    for spec in (0.95, 0.99):
        thr = np.quantile(neg, spec)
        out[f'sens_at_spec{int(spec*100)}'] = round(float((s[y == 1] > thr).mean()), 3)
    return out


groups = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m: continue
    key = (m.group(1), int(m.group(2)), int(m.group(3)), bool(m.group(4)))
    groups.setdefault(key, []).append((int(m.group(5)), f))

rows, oprows, calrows, relrows = [], [], [], []
for (name, L, MNA, region), runs in sorted(groups.items()):
    acc, per_seed, pa_raw = {}, {}, []
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        mp = np.nan_to_num(d.multipass_pa.values, nan=0.0)
        sc = baselines(d, name, seed)
        sc['MultiPASS'] = mp
        sc['nearest_pos'] = nearest_pos_scores(d, seed)
        y = d.y.values
        pa_raw.append(mp)
        for k, v in sc.items():
            acc.setdefault(k, []).append(stats.rankdata(v) / len(v))
            per_seed.setdefault(k, []).append(roc_auc_score(y, v))
    avg = {k: np.mean(v, 0) for k, v in acc.items()}
    setting = 'region' if region else 'full'
    mp_avg = avg['MultiPASS']

    for k, v in avg.items():
        lo, hi = ci(y, v)
        line = dict(target=name, setting=setting, model=k,
                    n_pos=int(y.sum()), n_neg=int((1 - y).sum()), seeds=len(runs),
                    auc_mean_over_seeds=round(float(np.mean(per_seed[k])), 3),
                    auc_sd_over_seeds=round(float(np.std(per_seed[k])), 3),
                    auc_pooled=round(roc_auc_score(y, v), 3),
                    ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                    auprc=round(average_precision_score(y, v), 3))
        if k != 'MultiPASS':
            _, p = delong(y, mp_avg, v)
            line['delong_p_vs_multipass'] = round(p, 4)
        rows.append(line)
        oprows.append(dict(target=name, setting=setting, model=k,
                           n_pos=int(y.sum()), n_neg=int((1 - y).sum()), **op_points(y, v)))

    # calibration of the raw Pa, averaged over seeds
    pa = np.mean(pa_raw, 0)
    calrows.append(dict(target=name, setting=setting, n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                        prevalence=round(float(y.mean()), 4),
                        pa_mean=round(float(pa.mean()), 4),
                        pa_mean_on_pos=round(float(pa[y == 1].mean()), 4),
                        pa_mean_on_neg=round(float(pa[y == 0].mean()), 4),
                        brier=round(float(brier_score_loss(y, np.clip(pa, 0, 1))), 4),
                        brier_baserate=round(float(brier_score_loss(y, np.full(len(y), y.mean()))), 4)))
    edges = np.linspace(0, 1, 11)
    bins = np.clip(np.digitize(pa, edges) - 1, 0, 9)
    for b in range(10):
        msk = bins == b
        if not msk.any(): continue
        relrows.append(dict(target=name, setting=setting, bin=f'{edges[b]:.1f}-{edges[b+1]:.1f}',
                            n=int(msk.sum()), mean_pa=round(float(pa[msk].mean()), 4),
                            observed_rate=round(float(y[msk].mean()), 4)))

df = pd.DataFrame(rows)
mask = df.delong_p_vs_multipass.notna()
p = df.loc[mask, 'delong_p_vs_multipass'].values
order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)

O = f'{_SAVDR_ROOT}/multipass'
df.to_csv(f'{O}/summary_with_nearestpos.csv', index=False)
pd.DataFrame(oprows).to_csv(f'{O}/operating_points.csv', index=False)
pd.DataFrame(calrows).to_csv(f'{O}/calibration.csv', index=False)
pd.DataFrame(relrows).to_csv(f'{O}/calibration_reliability.csv', index=False)

print('=' * 104)
print('1. nearest_pos -- untrained, one feature: -(distance to nearest training-fold resistant position)')
print('=' * 104)
w = df.pivot_table(index=['target', 'setting'], columns='model', values='auc_pooled')
cols = [c for c in ['MultiPASS', 'nearest_pos', 'pos_only', 'full_gbm', 'subst_lr'] if c in w.columns]
print(w[cols].round(3).to_string())
npd = df[df.model == 'nearest_pos']
print(f'\nnearest_pos vs MultiPASS: mean delta = {(w["MultiPASS"] - w["nearest_pos"]).mean():+.3f} AUC, '
      f'MultiPASS ahead in {int((w["MultiPASS"] > w["nearest_pos"]).sum())} of {len(w)} cells')
print(f'DeLong vs MultiPASS, Holm-corrected across {m} tests:')
print(npd[['target', 'setting', 'auc_pooled', 'delong_p_vs_multipass', 'delong_p_holm']].to_string(index=False))
print(f'significant after Holm: {int((npd.delong_p_holm < 0.05).sum())} of {len(npd)}')

print('\n' + '=' * 104)
print('2. OPERATING POINTS')
print('=' * 104)
op = pd.DataFrame(oprows)
print(op[op.model.isin(['MultiPASS', 'nearest_pos', 'pos_only'])].to_string(index=False))

print('\n' + '=' * 104)
print('3. CALIBRATION of raw Pa (averaged over seeds)')
print('=' * 104)
print(pd.DataFrame(calrows).to_string(index=False))
print('\nreliability curve:')
print(pd.DataFrame(relrows).to_string(index=False))
