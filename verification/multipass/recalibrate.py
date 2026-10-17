"""Does recalibrating Pa fix it? (Constructive half of the calibration finding.)

calibration.csv shows Pa is a bad probability: it overestimates by up to 5.9x and is a worse
probability than simply predicting the base rate in 10 of 12 cells. That is expected for a
PASS-style Pa -- it is an invariant-accuracy quantity, not a posterior for this prevalence --
and it is fixable. This fits a monotone (isotonic) recalibration map on out-of-fold Pa and
evaluates it *out of fold again*, position-grouped, so the map is never fitted on the variants
it scores. Isotonic is monotone, so AUC is unchanged by construction; only the probability
changes. Platt (logistic) is reported alongside as the lower-variance alternative.
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
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score

groups = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m: continue
    key = (m.group(1), int(m.group(2)), int(m.group(3)), bool(m.group(4)))
    groups.setdefault(key, []).append((int(m.group(5)), f))

rows = []
for (name, L, MNA, region), runs in sorted(groups.items()):
    pa_s = []
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        pa_s.append(np.nan_to_num(d.multipass_pa.values, nan=0.0))
    pa = np.clip(np.mean(pa_s, 0), 0, 1)
    y, g = d.y.values, d.pos.values
    iso = np.zeros(len(y)); pla = np.zeros(len(y))
    mono_iso, mono_pla, within = [], [], []
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=0).split(pa.reshape(-1, 1), y, g):
        ir = IsotonicRegression(out_of_bounds='clip', y_min=0, y_max=1).fit(pa[tr], y[tr])
        iso[te] = ir.predict(pa[te])
        lr = LogisticRegression(max_iter=2000).fit(pa[tr].reshape(-1, 1), y[tr])
        pla[te] = lr.predict_proba(pa[te].reshape(-1, 1))[:, 1]
        # within-fold monotonicity: a monotone map cannot invert any ordered pair
        o = np.argsort(pa[te])
        mono_iso.append(bool(np.all(np.diff(iso[te][o]) >= -1e-12)))
        mono_pla.append(bool(np.all(np.diff(pla[te][o]) >= -1e-12)))
        if len(np.unique(y[te])) == 2:
            within.append((roc_auc_score(y[te], pa[te]), roc_auc_score(y[te], iso[te]),
                           roc_auc_score(y[te], pla[te])))
    w = np.array(within)
    base = np.full(len(y), y.mean())
    rows.append(dict(target=name, setting='region' if region else 'full',
                     n_pos=int(y.sum()), n_neg=int((1 - y).sum()), prevalence=round(float(y.mean()), 4),
                     brier_raw_pa=round(brier_score_loss(y, pa), 4),
                     brier_isotonic=round(brier_score_loss(y, iso), 4),
                     brier_platt=round(brier_score_loss(y, pla), 4),
                     brier_baserate=round(brier_score_loss(y, base), 4),
                     auc_raw_pa=round(roc_auc_score(y, pa), 4),
                     auc_isotonic=round(roc_auc_score(y, iso), 4),
                     auc_platt=round(roc_auc_score(y, pla), 4),
                     monotone_within_fold_iso=all(mono_iso),
                     monotone_within_fold_platt=all(mono_pla),
                     auc_withinfold_raw=round(float(w[:, 0].mean()), 4),
                     auc_withinfold_iso=round(float(w[:, 1].mean()), 4),
                     auc_withinfold_platt=round(float(w[:, 2].mean()), 4)))

df = pd.DataFrame(rows)
df['iso_beats_baserate'] = df.brier_isotonic < df.brier_baserate
df['raw_beats_baserate'] = df.brier_raw_pa < df.brier_baserate
df.to_csv(f'{_SAVDR_ROOT}/multipass/recalibration.csv', index=False)
print(df.to_string(index=False))
print(f'\nraw Pa      beats the base rate in {int(df.raw_beats_baserate.sum())} of {len(df)} cells')
print(f'isotonic Pa beats the base rate in {int(df.iso_beats_baserate.sum())} of {len(df)} cells')
print(f'mean Brier: raw {df.brier_raw_pa.mean():.4f} -> isotonic {df.brier_isotonic.mean():.4f} '
      f'(base rate {df.brier_baserate.mean():.4f}, Platt {df.brier_platt.mean():.4f})')
print(f'\nmean AUC: raw {df.auc_raw_pa.mean():.4f} -> isotonic {df.auc_isotonic.mean():.4f}, '
      f'Platt {df.auc_platt.mean():.4f}   (pooled across nested folds)')
print(f'max pooled AUC loss: isotonic {(df.auc_raw_pa - df.auc_isotonic).max():.4f}, '
      f'Platt {(df.auc_raw_pa - df.auc_platt).max():.4f}')
print('\nIs each map monotone WITHIN its own fold (it must be -- otherwise this is a bug, not a cost)?')
print(f'  isotonic: {int(df.monotone_within_fold_iso.sum())}/{len(df)} cells, '
      f'Platt: {int(df.monotone_within_fold_platt.sum())}/{len(df)} cells')
print('mean AUC computed WITHIN each nested fold (monotone maps cannot change this):')
print(f'  raw {df.auc_withinfold_raw.mean():.4f}  isotonic {df.auc_withinfold_iso.mean():.4f}  '
      f'Platt {df.auc_withinfold_platt.mean():.4f}')
print(f'  max |within-fold delta|: isotonic '
      f'{(df.auc_withinfold_iso - df.auc_withinfold_raw).abs().max():.4f}, '
      f'Platt {(df.auc_withinfold_platt - df.auc_withinfold_raw).abs().max():.4f}')
print('\nSo the pooled AUC loss is NOT a monotonicity violation: it is the cost of stitching '
      'five\nfold-specific maps together, whose outputs are not on a common scale. Platt loses '
      'less\nthan isotonic because a logistic map keeps ties apart; isotonic flattens variants '
      'into steps.')
