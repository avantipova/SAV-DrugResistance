"""Independent verification of the TKI-benchmark comparison."""
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

ok = fail = 0
def chk(name, got, want):
    global ok, fail
    good = got == want
    ok, fail = ok + good, fail + (not good)
    print(f'  [{"PASS" if good else "FAIL"}] {name:62s} got={got!r:>22} want={want!r}')

T = pd.read_csv(f'{_SAVDR_ROOT}/multipass/tki_overlap.csv')
d = pd.read_csv(f'{_SAVDR_ROOT}/external/mdrdb/core.tsv', sep='\t', low_memory=False)
src = d[(d.DATASET == 'TKI') & (d.TYPE == 'Single Substitution')]

print('\n1. The TKI slice of MdrDB')
chk('rows in MdrDB', len(src), 144)
chk('all of them carried over', len(T), 144)
chk('distinct substitutions', T['mut'].nunique(), 31)
chk('eight inhibitors', sorted(T.drug.unique()),
    ['axitinib', 'bosutinib', 'dasatinib', 'erlotinib', 'gefitinib', 'imatinib', 'nilotinib',
     'ponatinib'])
chk('MdrDB spreads them over four accessions', sorted(T.acc.unique()),
    ['P00519', 'P00520', 'Q835X8', 'UNDEFINED'])
chk('imatinib and ponatinib rows carry the MOUSE accession',
    sorted(T[T.acc == 'P00520'].drug.unique()), ['imatinib', 'ponatinib'])
chk('nilotinib rows carry Q835X8', sorted(T[T.acc == 'Q835X8'].drug.unique()), ['nilotinib'])
chk('erlotinib and gefitinib rows carry no accession',
    sorted(T[T.acc == 'UNDEFINED'].drug.unique()), ['erlotinib', 'gefitinib'])
chk('only 68 rows would survive a human-accession filter', int((T.acc == 'P00519').sum()), 68)

print('\n2. Drug identity, resolved by PubChem CID rather than by ligand code')
CID = {'5291': 'imatinib', '644241': 'nilotinib', '3062316': 'dasatinib', '6450551': 'axitinib',
       '24826799': 'ponatinib', '176870': 'erlotinib', '123631': 'gefitinib',
       '5328940': 'bosutinib'}
s = src.copy()
s['cid'] = s.CID.astype(str).str.replace(r'\.0$', '', regex=True)
chk('every TKI row has a CID in the mapping', int((~s.cid.isin(CID)).sum()), 0)
lig = s.assign(name=s.cid.map(CID))[['DRUG', 'name']].drop_duplicates()
chk('ligand code to drug name', dict(zip(lig.DRUG, lig.name)),
    {'AXI': 'axitinib', '1N1': 'dasatinib', 'DB8': 'bosutinib', '0LI': 'ponatinib',
     'STI': 'imatinib', 'NIL': 'nilotinib', 'AQ4': 'erlotinib', 'IRE': 'gefitinib'})
# imatinib, PubChem CID 5291, canonical SMILES begins with the 4-methylphenyl-benzamide head
chk('CID 5291 carries the imatinib SMILES',
    s[s.cid == '5291'].SMILES.iloc[0],
    'CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5')
chk('CID 5291 appears only as the ligand code STI',
    sorted(s[s.cid == '5291'].DRUG.unique()), ['STI'])

print('\n3. Numbering: every row checked against the HUMAN ABL1 sequence')
can = first_record(f'{_SAVDR_ROOT}/inputs/fasta/ABL1_ENST00000318560.faa')
iso = first_record(f'{_SAVDR_ROOT}/inputs/fasta/ABL1.faa')
chk('sequence lengths', (len(can), len(iso)), (1130, 1149))
bad = sum(not (1 <= int(r['mut'][1:-1]) <= len(can)
               and can[int(r['mut'][1:-1]) - 1] == r['mut'][0]) for _, r in T.iterrows())
chk('wild-type mismatches against human ABL1 after the shift', bad, 0)
chk('offsets used', sorted(T.offset.unique()), [0])
chk('so the TKI set is numbered on the canonical 1130-residue sequence',
    bool((T.offset == 0).all()), True)
tab = pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/ABL1_ENST00000318560.csv')
have = set(tab.aa.str.replace(r'^p\.', '', regex=True))
chk('rows present in our own variant table', int(T['mut'].isin(have).sum()), 143)
chk('the single absent substitution', sorted(set(T.loc[~T['mut'].isin(have), 'mut'])), ['E292L'])
chk('all 31 mutations lie in the kinase domain (242-493)',
    (int(T['mut'].str[1:-1].astype(int).min()), int(T['mut'].str[1:-1].astype(int).max())),
    (244, 486))

print('\n4. Our scores')
chk('per-drug scores present for 84 rows', int(T.pa_drug.notna().sum()), 84)
chk('those are the four inhibitors we label',
    sorted(T[T.pa_drug.notna()].drug.unique()),
    ['bosutinib', 'dasatinib', 'imatinib', 'nilotinib'])
chk('21 rows per labelled inhibitor',
    sorted(T[T.pa_drug.notna()].groupby('drug').size().unique()), [21])
chk('binary score present for 143 rows', int(T.pa_binary.notna().sum()), 143)
# per-drug Pa must equal the five-seed rank average the main pipeline produces
per = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/drug_ABL1_ENST00000318560_*_s*/scores_drug.csv')):
    sc = pd.read_csv(f)
    idx = sc.aa.str.replace(r'^p\.', '', regex=True)
    for c in sc.columns:
        if c.lower().startswith('pa::'):
            per.setdefault(c.split('::', 1)[1].strip().lower(), []).append(
                pd.Series(stats.rankdata(np.nan_to_num(sc[c].values, nan=0.0)) / len(sc), index=idx))
chk('five seeds per drug', sorted({len(v) for v in per.values()}), [5])
PD = {k: pd.concat(v, axis=1).mean(axis=1) for k, v in per.items()}
mism = 0
for _, r in T[T.pa_drug.notna()].iterrows():
    if abs(PD[r.drug].loc[r['mut']] - r.pa_drug) > 1e-12: mism += 1
chk('every per-drug score equals the five-seed rank average', mism, 0)

print('\n5. Correlations, recomputed independently')
R = pd.read_csv(f'{_SAVDR_ROOT}/multipass/tki_compare.csv')
md = T[T.pa_drug.notna()]
rx, ry = stats.rankdata(md.pa_drug.values), stats.rankdata(md.ddg.values)
chk('pooled drug-matched Spearman from ranks equals scipy',
    round(float(np.corrcoef(rx, ry)[0, 1]), 3),
    round(float(stats.spearmanr(md.pa_drug, md.ddg).statistic), 3))
chk('pooled drug-matched Spearman matches the report',
    round(float(stats.spearmanr(md.pa_drug, md.ddg).statistic), 3),
    float(R[R.panel.str.startswith('drug-matched')].spearman.iloc[0]))
for dr, want in [('dasatinib', 0.507), ('nilotinib', -0.406), ('imatinib', 0.171),
                 ('bosutinib', -0.038)]:
    g = md[md.drug == dr]
    chk(f'  {dr}', round(float(stats.spearmanr(g.pa_drug, g.ddg).statistic), 3), want)
chk('dasatinib is the only one nominally positive at p<0.05',
    sum(1 for dr in ['dasatinib', 'nilotinib', 'imatinib', 'bosutinib']
        if (lambda g: stats.spearmanr(g.pa_drug, g.ddg).pvalue < 0.05
            and stats.spearmanr(g.pa_drug, g.ddg).statistic > 0)(md[md.drug == dr])), 1)
chk('the pooled CI includes zero',
    bool(float(R[R.panel.str.startswith('drug-matched')].spearman_ci.iloc[0].split('..')[0]) < 0
         < float(R[R.panel.str.startswith('drug-matched')].spearman_ci.iloc[0].split('..')[1])),
    True)
chk('binary Pa over all 143 rows is slightly negative',
    round(float(stats.spearmanr(T.pa_binary, T.ddg, nan_policy='omit').statistic), 3), -0.129)
chk('our best pooled correlation is far below their worst published Pearson (0.54)',
    bool(max(abs(float(x)) for x in R.pearson) < 0.54), True)

print(f'\n{"=" * 108}\n{ok} passed, {fail} failed')
sys.exit(1 if fail else 0)
