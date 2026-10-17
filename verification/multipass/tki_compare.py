"""MultiPASS against the TKI benchmark of SPLDExtraTrees (Zhu et al., Brief Bioinform 2022).

SPLDExtraTrees predicts the change in inhibitor binding affinity caused by a kinase mutation,
from 146 structure-derived features (Rosetta REF15, FoldX, PLIP, AutoDock Vina). It is NOT in the
manuscript's Table S3 and has no published code, but its test set -- the "TKI dataset", 144
measured ddG values over ABL kinase point mutations and eight tyrosine-kinase inhibitors -- is
redistributed inside MdrDB, so our scores can be put beside its published correlations.

Published performance on this TKI dataset (their Table, mean over 20 repetitions):
    scenario 1  trained on Platinum only, no kinase data   RMSE 0.73  Pearson 0.54  AUPRC 0.50
    scenario 2  Platinum + 2 TKI samples per inhibitor     RMSE 0.70  Pearson 0.60  AUPRC 0.55
    scenario 3  trained and tested on the TKI dataset      RMSE 0.58  Pearson 0.74  AUPRC 0.56
    Rosetta REF15, for reference                           RMSE 0.72  Pearson 0.67  AUPRC 0.53

Scenario 1 is the fair comparison for us: neither model has seen kinase ddG data. Scenario 3 is
their best case and is cross-validated on the very set it is scored on.

**A data-quality finding about MdrDB that this analysis had to fix.** MdrDB spreads the TKI set
over four UniProt accessions depending on the template structure used: the axitinib, bosutinib
and dasatinib rows carry P00519 (human ABL1), imatinib and ponatinib carry **P00520, which is
mouse Abl1**, nilotinib carries **Q835X8**, and the erlotinib and gefitinib rows carry
**UNDEFINED**. All 144 rows are in fact human ABL1 substitutions: every one of them was checked
against the human sequence here and all 144 wild-type residues match. Filtering MdrDB by human
accession -- the conservative thing to do, and what `mdrdb_compare.py` does -- therefore keeps
only 68 of the 144 TKI rows and 3 of the 8 inhibitors. This script uses all 144, justified by the
wild-type check rather than by the accession label.
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

import glob, sys, numpy as np, pandas as pd
from scipy import stats
exec(open('summarize.py').read().split('rows = []')[0])

T = pd.read_csv(f'{_SAVDR_ROOT}/external/mdrdb/tki_verified.csv')
print(f'TKI rows, wild-type-verified against human ABL1: {len(T)}; '
      f'{T["mut"].nunique()} distinct mutations; drugs {sorted(T.drug.unique())}')

# our per-drug Pa for ABL1, rank-averaged over the five fold-assignment seeds, nothing refitted
per = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/drug_ABL1_ENST00000318560_*_s*/scores_drug.csv')):
    s = pd.read_csv(f)
    idx = s.aa.str.replace(r'^p\.', '', regex=True)
    for c in s.columns:
        if not c.lower().startswith('pa::'): continue
        name = c.split('::', 1)[1].strip().lower()
        v = pd.Series(stats.rankdata(np.nan_to_num(s[c].values, nan=0.0)) / len(s), index=idx)
        per.setdefault(name, []).append(v)
PD = {k: pd.concat(v, axis=1).mean(axis=1) for k, v in per.items()}
print(f'per-drug ABL1 scores available for: {sorted(PD)} ({len(per["imatinib"])} seeds each)')

# and the binary resistance Pa from the main full-setting run
bin_per = []
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_ABL1_ENST00000318560_17_9_s*/scores.csv')):
    s = pd.read_csv(f)
    bin_per.append(pd.Series(stats.rankdata(np.nan_to_num(s.multipass_pa.values, nan=0.0)) / len(s),
                             index=s.aa.str.replace(r'^p\.', '', regex=True)))
BIN = pd.concat(bin_per, axis=1).mean(axis=1)

T['pa_drug'] = [PD[r.drug].get(r['mut'], np.nan) if r.drug in PD else np.nan
                for _, r in T.iterrows()]
T['pa_binary'] = [BIN.get(m, np.nan) for m in T['mut']]
T.to_csv(f'{_SAVDR_ROOT}/multipass/tki_overlap.csv', index=False)


def corr(x, y, label, B=4000, seed=0):
    m = ~(pd.isna(x) | pd.isna(y))
    x, y = np.asarray(x, float)[m], np.asarray(y, float)[m]
    if len(x) < 6: return None
    sp, pe = stats.spearmanr(x, y), stats.pearsonr(x, y)
    rng = np.random.default_rng(seed)
    bs = [stats.spearmanr(x[i], y[i]).statistic
          for i in (rng.integers(0, len(x), len(x)) for _ in range(B))]
    bs = [v for v in bs if not np.isnan(v)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(panel=label, n=len(x), spearman=round(float(sp.statistic), 3),
                spearman_p=float(f'{sp.pvalue:.3g}'), spearman_ci=f'{lo:.3f}..{hi:.3f}',
                pearson=round(float(pe.statistic), 3), pearson_p=float(f'{pe.pvalue:.3g}'))


rows = []
md = T[T.pa_drug.notna()]
rows.append(corr(md.pa_drug, md.ddg, 'drug-matched Pa, 4 shared inhibitors'))
for dr, g in md.groupby('drug'):
    rows.append(corr(g.pa_drug, g.ddg, f'  {dr} (drug-matched Pa)'))
rows.append(corr(T.pa_binary, T.ddg, 'binary resistance Pa, all 8 inhibitors'))
sub = T[T.drug.isin(PD)]
rows.append(corr(sub.pa_binary, sub.ddg, 'binary resistance Pa, the same 4 inhibitors'))
oth = T[~T.drug.isin(PD)]
rows.append(corr(oth.pa_binary, oth.ddg, 'binary Pa, the 4 inhibitors we do not label'))
R = pd.DataFrame([r for r in rows if r])
pd.set_option('display.width', 240)
print('\n' + '=' * 110)
print("Our scores against the measured ddG of the TKI dataset")
print('=' * 110)
print(R.to_string(index=False))
R.to_csv(f'{_SAVDR_ROOT}/multipass/tki_compare.csv', index=False)

print('\nSPLDExtraTrees on the same 144 ddG values (published, not recomputed):')
print('  scenario 1 (no kinase data in training)  Pearson 0.54   <- the fair comparison')
print('  scenario 2 (+2 TKI samples per drug)     Pearson 0.60')
print('  scenario 3 (trained on the TKI set)      Pearson 0.74   <- their best case')
print('  Rosetta REF15                            Pearson 0.67')
print('\nNote: they predict ddG by regression and report Pearson over all 144 rows; our score is')
print('a resistance probability, so Spearman is the honest statistic for us and Pearson is given')
print('only for comparability. Our per-drug scores exist for 4 of their 8 inhibitors.')
