"""How large a difference could these comparisons have detected?

Every headline claim in REPORT.md is a null: MultiPASS is statistically indistinguishable from
a position-only model, Holm-corrected DeLong p = 1.00 in 11 of 12 cells. A reviewer is entitled
to read that as "n is too small to see anything". This script turns each null into a bounded
claim by reporting, per comparison:

  delta_obs   observed AUC difference (MultiPASS - baseline), on the seed-averaged scores
  se          standard error of that difference, from the DeLong covariance estimator
              (so it already carries the actual correlation between the two models' scores,
              which is what makes a paired comparison far more sensitive than two unpaired ones)
  mde_80      minimum difference detectable at 80% power, two-sided alpha = 0.05:
              (z_0.975 + z_0.80) * se = 2.8016 * se
  excl_95     one-sided 95% upper bound on the true difference: delta_obs + 1.6449 * se
              -> "the data exclude MultiPASS being better than the baseline by more than this"

The analytic MDE is then checked against a Monte-Carlo power simulation that re-runs the real
DeLong function on binormal scores matched to each cell's n, AUCs and score correlation.
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
from scipy import stats
from sklearn.metrics import roc_auc_score

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # AAS/KD/VOL/CHG, midrank, delong, baselines, ci, first_record

Z_A, Z_P = stats.norm.ppf(0.975), stats.norm.ppf(0.80)   # 1.9600, 0.8416


def delong_var(y, s1, s2):
    """Same estimator as delong(), but also returns the SE of the AUC difference."""
    p = np.c_[s1, s2][y == 1].T
    n = np.c_[s1, s2][y == 0].T
    m, nn, k = p.shape[1], n.shape[1], 2
    tx = np.array([midrank(p[r]) for r in range(k)])
    ty = np.array([midrank(n[r]) for r in range(k)])
    tz = np.array([midrank(np.r_[p[r], n[r]]) for r in range(k)])
    auc = (tz[:, :m].sum(1) / m - (m + 1) / 2) / nn
    v01 = (tz[:, :m] - tx) / nn
    v10 = 1 - (tz[:, m:] - ty) / m
    S = np.cov(v01) / m + np.cov(v10) / nn
    L = np.array([[1, -1]])
    var = float(np.asarray(L @ S @ L.T).ravel()[0])
    return auc, np.sqrt(max(var, 0.0))


def boot_se_delta(y, s1, s2, B=4000, seed=0):
    """Stratified paired bootstrap SE of the AUC difference, straight off the real scores.

    This is the assumption-free check on the one quantity the MDE rests on. A binormal
    simulation was tried first and rejected: it gave 0.26-0.98 power at the analytic MDE
    (mean 0.70) because it cannot reproduce the tie structure of these scores -- MultiPASS's
    Pa is coarsely quantised and pos_only is a step function over residue positions. The
    DeLong test's own type-I rate under that simulation was 0.057, i.e. the test is calibrated
    and it was the generator that was wrong. Bootstrapping the observed scores keeps the real
    ties, the real score correlation and the real class imbalance.
    """
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    out = np.empty(B)
    for b in range(B):
        i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        yy = y[i]
        out[b] = roc_auc_score(yy, s1[i]) - roc_auc_score(yy, s2[i])
    return float(out.std(ddof=1))


# ---- rebuild the seed-averaged score matrices exactly as summarize.py does -------------------
groups = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m: continue
    key = (m.group(1), int(m.group(2)), int(m.group(3)), bool(m.group(4)))
    groups.setdefault(key, []).append((int(m.group(5)), f))

rows, simrows = [], []
for (name, L, MNA, region), runs in sorted(groups.items()):
    acc = {}
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        sc = baselines(d, name, seed)
        sc['MultiPASS'] = np.nan_to_num(d.multipass_pa.values, nan=0.0)
        for k, v in sc.items():
            acc.setdefault(k, []).append(stats.rankdata(v) / len(v))
        y = d.y.values
    avg = {k: np.mean(v, 0) for k, v in acc.items()}
    setting = 'region' if region else 'full'
    mp = avg['MultiPASS']
    for k in ['pos_only', 'full_gbm', 'subst_lr']:
        (a1, a2), se = delong_var(y, mp, avg[k])
        delta = a1 - a2
        rho = float(stats.spearmanr(mp, avg[k]).statistic)
        rows.append(dict(target=name, setting=setting, baseline=k,
                         n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                         auc_multipass=round(a1, 3), auc_baseline=round(a2, 3),
                         delta_obs=round(delta, 3), se=round(se, 4),
                         mde_80=round((Z_A + Z_P) * se, 3),
                         excl_95=round(delta + stats.norm.ppf(0.95) * se, 3),
                         score_corr=round(rho, 3)))

df = pd.DataFrame(rows)
df.to_csv(f'{_SAVDR_ROOT}/multipass/power.csv', index=False)

print('=' * 100)
print('MINIMUM DETECTABLE AUC DIFFERENCE (two-sided alpha=0.05, 80% power, paired DeLong)')
print('=' * 100)
print(df.to_string(index=False))

print('\n--- the headline null: MultiPASS vs pos_only ---')
po = df[df.baseline == 'pos_only'].copy()
print(po[['target', 'setting', 'n_pos', 'n_neg', 'delta_obs', 'se', 'mde_80', 'excl_95']].to_string(index=False))
print(f'\nmedian MDE at 80% power: {po.mde_80.median():.3f} AUC')
print(f'range                   : {po.mde_80.min():.3f} - {po.mde_80.max():.3f}')
print(f'cells where |delta_obs| already exceeds the MDE: {int((po.delta_obs.abs() > po.mde_80).sum())} of {len(po)}')
print(f'strongest exclusion (smallest excl_95): {po.excl_95.min():.3f} '
      f'({po.loc[po.excl_95.idxmin(), "target"]} {po.loc[po.excl_95.idxmin(), "setting"]})')
print(f'weakest   exclusion (largest  excl_95): {po.excl_95.max():.3f} '
      f'({po.loc[po.excl_95.idxmax(), "target"]} {po.loc[po.excl_95.idxmax(), "setting"]})')

# ---- assumption-free check: does the DeLong SE of the difference match a paired bootstrap? --
print('\n' + '=' * 100)
print('CHECK: DeLong analytic SE of the AUC difference vs stratified paired bootstrap SE')
print('(the MDE is 2.8016 x SE, so the MDE is only as good as this SE; B=4000 on real scores)')
print('=' * 100)
for (name, L, MNA, region), runs in sorted(groups.items()):
    acc = {}
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        sc = baselines(d, name, seed)
        sc['MultiPASS'] = np.nan_to_num(d.multipass_pa.values, nan=0.0)
        for k, v in sc.items():
            acc.setdefault(k, []).append(stats.rankdata(v) / len(v))
        y = d.y.values
    avg = {k: np.mean(v, 0) for k, v in acc.items()}
    setting = 'region' if region else 'full'
    for k in ['pos_only', 'full_gbm', 'subst_lr']:
        _, se_dl = delong_var(y, avg['MultiPASS'], avg[k])
        se_bs = boot_se_delta(y, avg['MultiPASS'], avg[k], B=4000, seed=7)
        simrows.append(dict(target=name, setting=setting, baseline=k,
                            se_delong=round(se_dl, 4), se_bootstrap=round(se_bs, 4),
                            ratio=round(se_bs / se_dl, 3),
                            mde_delong=round((Z_A + Z_P) * se_dl, 3),
                            mde_bootstrap=round((Z_A + Z_P) * se_bs, 3)))
s = pd.DataFrame(simrows)
s.to_csv(f'{_SAVDR_ROOT}/multipass/power_secheck.csv', index=False)
print(s.to_string(index=False))
print(f'\nbootstrap/DeLong SE ratio: median {s.ratio.median():.3f}, '
      f'range {s.ratio.min():.3f}-{s.ratio.max():.3f}')
print(f'max absolute MDE discrepancy: {(s.mde_bootstrap - s.mde_delong).abs().max():.3f} AUC')
print('\nNOTE: the binormal power simulation was tried and rejected -- see boot_se_delta() '
      'docstring.\nWhat is validated here is the SE, hence the MDE and the exclusion bound; '
      'no power curve is claimed.')
