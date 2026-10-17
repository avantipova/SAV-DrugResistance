"""Independent verification of the MdrDB comparison.

Each check recomputes a reported number by a different route, or tests an invariant that would
break if `mdrdb_compare.py` were wrong. The numbering is checked hardest, because MdrDB turned
out to be internally inconsistent for ABL1.
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

import hashlib, sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from scipy import stats

ok = fail = 0
def chk(name, got, want):
    global ok, fail
    good = got == want
    ok, fail = ok + good, fail + (not good)
    print(f'  [{"PASS" if good else "FAIL"}] {name:64s} got={got!r:>24} want={want!r}')

TSV = f'{_SAVDR_ROOT}/external/mdrdb/core.tsv'
M = pd.read_csv(f'{_SAVDR_ROOT}/multipass/mdrdb_overlap.csv')
d = pd.read_csv(TSV, sep='\t', low_memory=False)
THR = 1.36
src = open('summarize.py').read()
exec(src.split('rows = []')[0])

print('\n1. The downloaded file')
raw = open(TSV, 'rb').read()
chk('byte size matches what the server sent', len(raw), 5032052)
print(f'         sha256: {hashlib.sha256(raw).hexdigest()}')
chk('row count', len(d), 4292)
chk('column count', len(d.columns), 167)
chk('proteins / drugs / mutations in the core set',
    (d.UNIPROT_ID.nunique(), d.DRUG.nunique(), d.MUTATION.nunique()), (240, 410, 2503))
chk('the paper\'s headline 240 proteins and 2503 mutations reproduce',
    (d.UNIPROT_ID.nunique(), d.MUTATION.nunique()), (240, 2503))
chk('every row carries a ddG', int(d['DDG.EXP'].isna().sum()), 0)
chk('sources present', sorted(d.DATASET.unique()),
    ['AIMMS', 'Depmap', 'GDSC', 'KinaseMD', 'Platinum', 'RET', 'TKI'])

print('\n2. The slice we used')
OURS = {'P00519': 'ABL1_ENST00000318560', 'Q9UM73': 'ALK', 'P00533': 'EGFR',
        'P10721': 'KIT', 'Q06187': 'BTK_ENST00000308731', 'P36888': 'FLT3',
        'P08581': 'MET', 'P03372': 'ESR1_ENST00000206249', 'Q99835': 'SMO',
        'P49902': 'NT5C2', 'O75874': 'IDH1_ENST00000415913'}
s = d[(d.UNIPROT_ID.isin(OURS)) & (d.TYPE == 'Single Substitution')]
chk('single-substitution rows on our targets', len(s), 647)
chk('they span 9 of our 11 proteins (no NT5C2, no IDH1)',
    sorted({OURS[a].split('_')[0] for a in s.UNIPROT_ID.unique()}),
    ['ABL1', 'ALK', 'BTK', 'EGFR', 'ESR1', 'FLT3', 'KIT', 'MET', 'SMO'])
chk('matched rows', len(M), 411)
chk('rows above MdrDB\'s threshold, before matching', int((s['DDG.EXP'] > THR).sum()), 39)

print('\n3. Numbering: every matched row\'s wild-type residue must match the reference sequence')
bad = 0
for tgt, sub in M.groupby('target'):
    seq = first_record(f'{_SAVDR_ROOT}/inputs/fasta/{tgt}.faa')
    for _, r in sub.iterrows():
        p = int(r['mut'][1:-1])
        if not (1 <= p <= len(seq)) or seq[p - 1] != r['mut'][0]:
            bad += 1
chk('wild-type mismatches among matched rows', bad, 0)
chk('the offset applied is always 0 or +19', sorted(M.offset_used.unique()), [0, 19])
chk('the +19 rows are ABL1 only', sorted(M[M.offset_used == 19].gene.unique()), ['ABL1'])
n19 = M[M.offset_used == 19].mdrdb_mut.nunique()
# nine distinct ABL1 substitutions resolve only on isoform b. (An earlier version of this check
# expected five: that count came from a wild-type map built out of our variant table, which
# covers only positions carrying a COSMIC missense variant. Against the full reference sequence
# the correct number is nine.)
chk('ABL1 rows numbered on isoform b (+19): distinct mutations', n19, 9)
print(f'         those mutations: {sorted(M[M.offset_used == 19].mdrdb_mut.unique())}')
chk('MdrDB is therefore internally inconsistent for ABL1',
    bool(0 in set(M[M.gene == "ABL1"].offset_used) and 19 in set(M[M.gene == "ABL1"].offset_used)),
    True)
# the shift must be exactly the length difference between the two ABL1 isoforms
a_can = first_record(f'{_SAVDR_ROOT}/inputs/fasta/ABL1_ENST00000318560.faa')
a_iso = first_record(f'{_SAVDR_ROOT}/inputs/fasta/ABL1.faa')
chk('ABL1 isoform lengths are 1130 and 1149', (len(a_can), len(a_iso)), (1130, 1149))
chk('their difference is the +19 offset', len(a_iso) - len(a_can), 19)
# and a rescued mutation must sit on the isoform-b sequence at its own number
for v in sorted(M[M.offset_used == 19].mdrdb_mut.unique()):
    p = int(v[1:-1])
    chk(f'  {v} matches isoform b at position {p}', a_iso[p - 1], v[0])

print('\n4. The collapse to one row per mutation')
agg = (M.groupby(['gene', 'target', 'mut'])
        .agg(ddg_max=('ddg', 'max'), y_mdrdb=('y_mdrdb', 'max'), y_ours=('y_ours', 'first'),
             loto_pa=('loto_pa', 'first'), cv_rank=('cv_rank', 'first')).reset_index())
chk('distinct (gene, mutation) pairs', len(agg), 168)
chk('resistant by MdrDB\'s threshold after collapsing', int(agg.y_mdrdb.sum()), 18)
chk('y_mdrdb is exactly ddg_max > 1.36',
    int((agg.y_mdrdb != (agg.ddg_max > THR).astype(int)).sum()), 0)
chk('ddg_max really is the maximum over that mutation\'s rows',
    bool(all(agg.set_index(['gene', 'mut']).ddg_max.loc[k] == g.ddg.max()
             for k, g in M.groupby(['gene', 'mut']))), True)
chk('one score per mutation: no mutation carries two different LOTO values',
    int(M.groupby(['gene', 'mut']).loto_pa.nunique().max()), 1)

print('\n5. Correlations and AUCs, recomputed independently')
for col, n, want_rho, want_auc in [('cv_rank', 120, 0.317, 0.710), ('loto_pa', 168, 0.220, 0.716)]:
    a = agg[agg[col].notna()]
    x = np.nan_to_num(a[col].values, nan=0.0)
    y = a.ddg_max.values
    chk(f'{col}: n', len(a), n)
    # Spearman as Pearson on ranks, not via scipy.spearmanr
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    rho_manual = float(np.corrcoef(rx, ry)[0, 1])
    chk(f'{col}: Spearman from ranks equals scipy',
        round(rho_manual, 3), round(float(stats.spearmanr(x, y).statistic), 3))
    chk(f'{col}: Spearman matches the report', round(rho_manual, 3), want_rho)
    yb = (y > THR).astype(int)
    u = stats.mannwhitneyu(x[yb == 1], x[yb == 0], alternative='greater').statistic
    chk(f'{col}: AUC from Mann-Whitney U equals sklearn',
        round(u / (yb.sum() * (1 - yb).sum()), 3), round(roc_auc_score(yb, x), 3))
    chk(f'{col}: AUC matches the report', round(roc_auc_score(yb, x), 3), want_auc)

print('\n6. Label agreement with ours, and what the score adds over the label')
ct = pd.crosstab(agg.y_mdrdb, agg.y_ours).reindex(index=[0, 1], columns=[0, 1]).fillna(0).values
chk('cross-tabulation', [int(v) for v in ct.ravel()], [122, 28, 7, 11])
n = ct.sum(); po = (ct[0, 0] + ct[1, 1]) / n
pe = ((ct[0].sum() * ct[:, 0].sum()) + (ct[1].sum() * ct[:, 1].sum())) / n ** 2
chk('raw agreement', round(po, 3), 0.792)
chk('Cohen kappa', round((po - pe) / (1 - pe), 3), 0.280)
chk('kappa is far below DRSP\'s 0.995, i.e. this source is independent',
    (po - pe) / (1 - pe) < 0.5, True)
yb = (agg.ddg_max > THR).astype(int).values
chk('our binary label alone, as a predictor of MdrDB resistance',
    round(roc_auc_score(yb, agg.y_ours.values.astype(float)), 3), 0.712)
chk('MultiPASS LOTO, same task',
    round(roc_auc_score(yb, np.nan_to_num(agg.loto_pa.values, nan=0.0)), 3), 0.716)
chk('the score does not beat the label by more than 0.01',
    round(roc_auc_score(yb, np.nan_to_num(agg.loto_pa.values, nan=0.0))
          - roc_auc_score(yb, agg.y_ours.values.astype(float)), 3) <= 0.01, True)

print('\n7. Our labels against the measured ddG')
p1 = agg[agg.y_ours == 1].ddg_max.values
p0 = agg[agg.y_ours == 0].ddg_max.values
chk('group sizes', (len(p1), len(p0)), (39, 129))
chk('means', (round(float(p1.mean()), 3), round(float(p0.mean()), 3)), (0.838, 0.169))
chk('medians', (round(float(np.median(p1)), 3), round(float(np.median(p0)), 3)), (0.680, 0.147))
u = stats.mannwhitneyu(p1, p0, alternative='greater')
chk('Mann-Whitney U', int(u.statistic), 3498)
chk('p-value', float(f'{u.pvalue:.3g}'), 0.000112)
# the same claim by a label-permutation test, which assumes nothing about the distribution
rng = np.random.default_rng(0)
allv = np.concatenate([p1, p0]); k = len(p1)
obs = p1.mean() - p0.mean()
cnt = sum((lambda z: z[:k].mean() - z[k:].mean())(rng.permutation(allv)) >= obs
          for _ in range(20000))
chk('permutation test agrees (p < 0.001 over 20000 shuffles)', cnt / 20000 < 0.001, True)
print(f'         permutation p = {cnt / 20000:.5f}, observed difference in means {obs:.3f}')

print('\n8. Source-level sanity: measurement-based subsets should be less noisy than GDSC')
for ds in ['TKI', 'AIMMS', 'GDSC']:
    keep = M[M.dataset == ds][['gene', 'target', 'mut']].drop_duplicates()
    a = agg.merge(keep, on=['gene', 'target', 'mut'])
    a = a[a.loto_pa.notna()]
    rho = float(stats.spearmanr(np.nan_to_num(a.loto_pa.values, nan=0.0), a.ddg_max.values).statistic)
    print(f'         {ds:6s} n={len(a):4d}  rho={rho:+.3f}')
chk('GDSC is the largest source in the slice',
    M.dataset.value_counts().idxmax(), 'GDSC')
chk('GDSC ddG is centred near zero (IC50-derived, noisy)',
    abs(float(M[M.dataset == "GDSC"].ddg.median())) < 0.2, True)
chk('measurement-based subsets correlate better than GDSC',
    all(float(stats.spearmanr(
            np.nan_to_num(agg.merge(M[M.dataset == ds][['gene','target','mut']].drop_duplicates(),
                                    on=['gene','target','mut']).loto_pa.values, nan=0.0),
            agg.merge(M[M.dataset == ds][['gene','target','mut']].drop_duplicates(),
                      on=['gene','target','mut']).ddg_max.values).statistic)
        > float(stats.spearmanr(
            np.nan_to_num(agg.merge(M[M.dataset == 'GDSC'][['gene','target','mut']].drop_duplicates(),
                                    on=['gene','target','mut']).loto_pa.values, nan=0.0),
            agg.merge(M[M.dataset == 'GDSC'][['gene','target','mut']].drop_duplicates(),
                      on=['gene','target','mut']).ddg_max.values).statistic)
        for ds in ['TKI', 'AIMMS']), True)

print(f'\n{"=" * 112}\n{ok} passed, {fail} failed')
sys.exit(1 if fail else 0)
