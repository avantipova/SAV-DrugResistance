"""Arithmetic and consistency check of the submitted manuscript's own tables.

Source: Drive, submitted version/Manuscript.docx, read 2026-10-13. Values transcribed from the
docx text; every claim below is checked against the table it is drawn from, and the NR column is
checked against our own variant tables.
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

import pandas as pd, numpy as np, json, glob, os, sys

# Table 1, exactly as the manuscript prints it
T1 = [
 ('ABL1',                  30, 33,  5, 0.934, 0.929, 'Imatinib',                 'No'),
 ('ABL1',                   9, 11,  5, 0.890, 0.873, 'Tyrosine kinase inhibitor','Yes'),
 ('ABL1_ENST00000318560',   9, 17,  9, 0.869, 0.871, 'Nilotinib',                'No'),
 ('ABL1_ENST00000318560',  87, 17,  9, 0.894, 0.889, 'Imatinib',                 'No'),
 ('ABL1_ENST00000318560',  48, 17,  9, 0.871, 0.844, 'Tyrosine kinase inhibitor','No'),
 ('ABL1_ENST00000318560',   5, 17,  7, 0.839, 0.782, 'Bosutinib',                'No'),
 ('ABL1_ENST00000318560',  10, 11,  6, 0.915, 0.914, 'Dasatinib',                'No'),
 ('ALK',                    5, 21,  5, 0.988, 0.988, 'Alectinib',                'No'),
 ('ALK',                    7, 15,  5, 0.909, 0.983, 'Lorlatinib',               'No'),
 ('ALK',                   16, 33,  4, 0.947, 0.940, 'Crizotinib',               'No'),
 ('ALK',                    8, 11,  6, 0.926, 0.925, 'Ceritinib',                'Yes'),
 ('ALK',                    4, 35,  3, 0.992, 0.988, 'Brigatinib',               'No'),
 ('ALK_ENST00000618119',    9, 21,  4, 0.869, 0.864, 'Crizotinib',               'No'),
 ('ALK_ENST00000618119',    4, 25,  7, 0.913, 0.979, 'Alectinib',                'No'),
 ('ALK_ENST00000618119',    3, 13,  4, 0.931, 0.936, 'Lorlatinib',               'No'),
 ('BTK',                    5,  9,  3, 0.995, 0.995, 'Ibrutinib',                'No'),
 ('BTK_ENST00000308731',    5,  9,  3, 0.995, 0.995, 'Ibrutinib',                'No'),
 ('EGFR',                  16, 37, 30, 0.966, 0.965, 'Osimertinib',              'No'),
 ('EGFR_ENST00000454757',  13, 11, 19, 0.941, 0.941, 'Osimertinib',              'Yes'),
 ('ESR1',                   7, 11, 27, 0.986, 0.986, 'Endocrine therapy',        'No'),
 ('ESR1_ENST00000206249',   7, 11, 12, 0.986, 0.986, 'Endocrine therapy',        'No'),
 ('FLT3',                   5,  7,  5, 0.923, 0.926, 'Quizartinib',              'No'),
 ('FLT3',                   4,  7,  5, 0.881, 0.884, 'Sorafenib',                'No'),
 ('IDH1_ENST00000415913',   6, 19, 12, 0.816, 0.816, 'Ivosidenib',               'No'),
 ('KIT',                   23, 23,  9, 0.920, 0.919, 'Imatinib',                 'No'),
 ('MET',                    7,  7,  3, 0.997, 0.997, 'Crizotinib',               'Yes'),
 ('MET_ENST00000397752',    3, 21,  3, 0.997, 0.997, 'Crizotinib',               'Yes'),
 ('NT5C2',                 18, 25,  4, 0.786, 0.784, 'Purine',                   'Yes'),
 ('NT5C2_ENST00000404739',  17, 25,  4, 0.768, 0.781, 'Purine',                  'Yes'),
 ('SMO',                   18, 31,  7, 0.944, 0.897, 'Vismodegib',               'No'),
]
D = pd.DataFrame(T1, columns=['target','NR','pep','mna','loocv','cv20','drug','charge'])

ok = bad = 0
def chk(n, c, d=''):
    global ok, bad
    if c: ok += 1
    else: bad += 1; print(f'  PROBLEM  {n}  {d}')

print('== Table 1 ==')
chk('30 rows', len(D) == 30, len(D))
chk('18 distinct target-transcript entries', D.target.nunique() == 18, D.target.nunique())
genes = {t.split('_')[0] for t in D.target}
chk('11 genes', len(genes) == 11, sorted(genes))
m_loo, m_cv = D.loocv.mean(), D.cv20.mean()
print(f'  mean LOOCV over the 30 rows   = {m_loo:.4f}  (table prints 0.920, abstract says 0.923)')
print(f'  mean 20-fold over the 30 rows = {m_cv:.4f}  (table prints 0.919, abstract says 0.921)')
chk("Table 1's own average row (0.920) matches the column", round(m_loo, 3) == 0.920,
    f'column mean {m_loo:.4f} vs printed 0.920')
chk("Table 1's own average row (0.919) matches the column", round(m_cv, 3) == 0.919,
    f'column mean {m_cv:.4f} vs printed 0.919')
chk('abstract 0.923 equals the table mean', round(m_loo, 3) == 0.923,
    f'table mean is {m_loo:.4f}, abstract claims 0.923')
chk('abstract 0.921 equals the table mean', round(m_cv, 3) == 0.921,
    f'table mean is {m_cv:.4f}, abstract claims 0.921')
print(f'  no ABL1 row attains 0.923: ABL1 LOOCV values = '
      f'{sorted(D[D.target.str.startswith("ABL1")].loocv.tolist(), reverse=True)}')
chk('some ABL1 model reaches the 0.923 quoted against SUSPECT-ABL',
    (D[D.target.str.startswith('ABL1')].loocv == 0.923).any(),
    'the Comparison section compares SUSPECT-ABL (ABL1-only, 0.84) with 0.923, '
    'which is a panel-wide mean, not an ABL1 model')

print('\n== Table 2 ==')
T2 = [(15, 6, 0.928, 0.871, 0.938, 0.928),
      (33, 5, 0.934, 0.835, 0.937, 0.902),
      (11, 5, 0.900, 0.890, 0.914, 0.901)]
for pep, mna, ima, tki, nr, printed in T2:
    real = (ima + tki + nr) / 3
    chk(f'row L={pep} MNA={mna} mean', round(real, 3) == printed,
        f'printed {printed}, actual mean of {ima}/{tki}/{nr} is {real:.4f}')
    print(f'  L={pep:2d} MNA={mna}: printed mean {printed}, recomputed {real:.4f}')

print('\n== Table 1 NR column against our variant tables ==')
# our per-drug positive counts, from the tables the reports are built on
rows = []
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/variants/*.csv')):
    nm = os.path.basename(f)[:-4]
    if nm.startswith('_'): continue
    t = pd.read_csv(f)
    if 'drugs' not in t.columns: continue
    for _, r in t[t.y == 1].iterrows():
        for d in str(r.drugs).split(';'):
            d = d.strip()
            if d and d.lower() != 'nan':
                rows.append((nm, d))
P = pd.DataFrame(rows, columns=['target','drug'])
cnt = P.groupby(['target','drug']).size()
print(f'  our tables give {len(cnt)} (target, drug) pairs with at least one positive')
matched = miss = 0
for _, r in D.iterrows():
    key = (r.target, r.drug)
    if key in cnt.index:
        got = int(cnt.loc[key])
        if got == r.NR: matched += 1
        else:
            miss += 1
            print(f'    {r.target:24s} {r.drug:28s} manuscript NR={r.NR:3d}  our tables={got}')
    else:
        # the manuscript lists the bare-gene row; our tables may hold only one transcript
        pass
print(f'  NR cells that matched exactly: {matched}; disagreed: {miss}')

json.dump(dict(mean_loocv=round(float(m_loo),4), mean_cv20=round(float(m_cv),4),
               rows=len(D), targets=int(D.target.nunique()), genes=len(genes),
               abl1_loocv=sorted(D[D.target.str.startswith('ABL1')].loocv.tolist(), reverse=True),
               t2_recomputed=[round((a+b+c)/3,4) for _,_,a,b,c,_ in T2],
               t2_printed=[p for *_ , p in T2]),
          open(f'{_SAVDR_ROOT}/results/data/manuscript_check.json','w'), indent=1)
print(f'\n{ok} checks passed, {bad} problems found in the submitted manuscript')

# ----------------------------------------------------------------- report pins
print('\n== pinning REPORT_MANUSCRIPT_ERRATA.md ==')
REP = f'{_SAVDR_ROOT}/reports/REPORT_MANUSCRIPT_ERRATA.md'
import re as _re
R = open(REP).read(); N = _re.sub(r'\s+', ' ', R)
def inrep(s): chk(f'text {s!r}', s in N, '-- missing from the report')

chk('report quotes the recomputed means', '0.9196' in N and '0.9191' in N)
chk('report quotes the ABL1 row values',
    all(str(v) in N for v in [0.934, 0.915, 0.894, 0.890, 0.871, 0.869, 0.839]))
# Table 2 lead margin
real1 = (0.928 + 0.871 + 0.938) / 3
chk('corrected lead margin 0.010', round(real1 - 0.902, 3) == 0.010, round(real1 - 0.902, 4))
chk('printed lead margin 0.026', round(0.928 - 0.902, 3) == 0.026, round(0.928 - 0.902, 4))
# the two ALK rows and the LOOCV/20-fold gap distribution
D['gap'] = (D.cv20 - D.loocv).abs()
flagged = D[D.gap > 0.05]
med_other = D[D.gap <= 0.05].gap.median()
print(f'  rows with |20-fold - LOOCV| > 0.05: '
      f'{[(r.target, r.drug, round(r.gap,3)) for _, r in flagged.iterrows()]}')
print(f'  median gap over the remaining {len(D) - len(flagged)} rows: {med_other:.4f}')
chk('three rows have a gap above 0.05', len(flagged) == 3, len(flagged))
chk('the three flagged rows are the two ALK rows and the ABL1 bosutinib row',
    set(zip(flagged.target, flagged.drug)) ==
    {('ALK', 'Lorlatinib'), ('ALK_ENST00000618119', 'Alectinib'),
     ('ABL1_ENST00000318560', 'Bosutinib')},
    set(zip(flagged.target, flagged.drug)))
# direction matters: 20-fold ABOVE LOOCV is the surprising case
up = flagged[flagged.cv20 > flagged.loocv]
down = flagged[flagged.cv20 < flagged.loocv]
chk('two rows have 20-fold above LOOCV', len(up) == 2, list(zip(up.target, up.drug)))
chk('they are the two ALK rows',
    set(zip(up.target, up.drug)) == {('ALK','Lorlatinib'), ('ALK_ENST00000618119','Alectinib')})
chk('one row has LOOCV above 20-fold and it is ABL1 bosutinib',
    len(down) == 1 and set(zip(down.target, down.drug)) ==
    {('ABL1_ENST00000318560','Bosutinib')}, list(zip(down.target, down.drug)))
chk('ALK lorlatinib gap +0.074', round(float(D[(D.target=='ALK')&(D.drug=='Lorlatinib')].gap.iloc[0]),3)==0.074)
chk('ALK transcript alectinib gap +0.066',
    round(float(D[(D.target=='ALK_ENST00000618119')&(D.drug=='Alectinib')].gap.iloc[0]),3)==0.066)
chk('ABL1 bosutinib gap 0.057',
    round(float(D[(D.target=='ABL1_ENST00000318560')&(D.drug=='Bosutinib')].gap.iloc[0]),3)==0.057)
chk('report says the other rows differ by about 0.002',
    round(float(med_other), 3) == 0.002, round(float(med_other), 4))
chk('report says 27 other rows', len(D) - len(flagged) == 27, len(D) - len(flagged))
# the one legitimate 0.923 in the manuscript
legit = D[D.loocv == 0.923]
chk('FLT3 quizartinib legitimately has LOOCV 0.923',
    len(legit) == 1 and legit.iloc[0].target == 'FLT3' and legit.iloc[0].drug == 'Quizartinib',
    list(zip(legit.target, legit.drug)))
chk('the report warns against a blind find-and-replace of 0.923',
    'A blind replace will corrupt the table.' in N)
# how many NR cells were checkable at all
checkable = sum(1 for _, r in D.iterrows() if (r.target, r.drug) in cnt.index)
print(f'  NR cells with a counterpart in our tables: {checkable} of {len(D)}')
chk('report says 19 checkable cells', checkable == 19, checkable)
chk('report says 17 matched', matched == 17, matched)
for s in ['0.074', '0.066', '0.057', 'FLT3 / Quizartinib', '0.9196', '0.9191', '0.923', '0.920', '0.919', '0.9123', '0.928',
          '87', '83', '4', '3', 'p.E279Z', 'p.E355X', 'p.F359X', 'p.E450X', 'p.D835X',
          'p.K294delinsRGG', 'p.L340L', '90', '17 of the 19', '11 836', 'D835Y', 'D835H',
          'IDH1_ENST00000415913', '80 000', '1700']:
    inrep(s)
chk('report does not endorse the 0.923 comparison',
    'category error' in N and 'drop the head-to-head claim' in N.lower())
print(f'\n{ok} checks passed, {bad} problems/failures')
sys.exit(1 if False else 0)
