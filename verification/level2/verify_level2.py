"""Verify every number in level2/REPORT_LEVEL2.md against its sources.

Sources of truth:
  level2_eval.csv / level2_eval.log      -- all AUC / CI / AUPRC / DeLong / Holm figures
  data/alphamissense/AF-*.csv (restaged) -- row counts, coverage, the pathogenic-share table
  esm2_650M_scores.csv                   -- ESM row count, coverage, canonical spot checks
  fold_*_s0 / fold_*_region_s0           -- our own variant sets

Written after the two pooled tables in that report were found to come from a pre-FASTA-fix run.
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

import glob, re, sys, numpy as np, pandas as pd

R = open(f'{_SAVDR_ROOT}/reports/REPORT_LEVEL2.md').read()
E = pd.read_csv(f'{_SAVDR_ROOT}/level2/level2_eval.csv')
LOG = open(f'{_SAVDR_ROOT}/level2/level2_eval.log').read()
AM = f'{_SAVDR_ROOT}/external/source-data/alphamissense'

ok = fail = 0
def chk(name, got, want):
    global ok, fail
    good = got == want
    ok, fail = ok + good, fail + (not good)
    print(f'  [{"PASS" if good else "FAIL"}] {name:56s} report={got!r:>24} source={want!r}')

ACC = {'ABL1_ENST00000318560': 'P00519', 'ALK': 'Q9UM73', 'EGFR': 'P00533',
       'KIT': 'P10721', 'NT5C2': 'P49902', 'SMO': 'Q99835'}
SHORT = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'EGFR': 'EGFR',
         'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}
# report label -> CSV model name
M = {'MultiPASS': 'MultiPASS', 'pos_only': 'pos_only', 'full_gbm': 'full_gbm',
     'subst_lr': 'subst_lr', 'REVEL': 'REVEL', 'AlphaMissense': 'AlphaMissense',
     'CADD': 'CADD', 'PolyPhen-2': 'PolyPhen2', 'SIFT (−)': 'SIFT_neg',
     'ESM-2 650M': 'ESM2_650M_neg', 'ESM-2 150M': 'ESM2_150M_neg'}

def cell(setting, tgt, model, col='auc'):
    r = E[(E.setting == setting) & (E.target == tgt) & (E.model == model)]
    assert len(r) == 1, (setting, tgt, model, len(r))
    return r.iloc[0][col]

def rows_of(block):
    out = []
    started = False
    for line in block.splitlines():
        if not line.strip().startswith('|'):
            if started: break
            continue
        f = [c.strip().replace('**', '') for c in line.strip('|').split('|')]
        if f[0] in ('model', 'target') or set(f[0]) <= {'-'}: continue
        started = True
        out.append(f)
    return out

# ---- 1. pooled table, full setting ---------------------------------------------------------
print('\n1. Pooled table, `full` setting')
blk = R.split('## Result 1')[1]
m = re.search(r'\((\d+) pos / (\d+) neg\)', blk)
chk('full pooled n_pos', int(m.group(1)), int(cell('full', 'POOLED', 'MultiPASS', 'n_pos')))
chk('full pooled n_neg', int(m.group(2)), int(cell('full', 'POOLED', 'MultiPASS', 'n_neg')))
seen = []
for f in rows_of(blk.split('|---|---|---|---|---|---|---|')[1]):
    mod = M[f[0]]; seen.append(mod)
    chk(f'full pooled {mod} auc', float(f[2]), round(float(cell('full', 'POOLED', mod)), 3))
    lo, hi = (float(x) for x in f[3].replace('–', '-').split('-'))
    chk(f'full pooled {mod} ci', (lo, hi), (round(float(cell('full', 'POOLED', mod, 'ci_lo')), 3),
                                            round(float(cell('full', 'POOLED', mod, 'ci_hi')), 3)))
    chk(f'full pooled {mod} auprc', float(f[4]),
        round(float(cell('full', 'POOLED', mod, 'auprc')), 3))
    if mod != 'MultiPASS':
        p = float(cell('full', 'POOLED', mod, 'delong_p_vs_multipass'))
        h = float(cell('full', 'POOLED', mod, 'delong_p_holm'))
        if f[5].startswith('<'):
            chk(f'full pooled {mod} p tiny', p < float(f[5][1:]), True)
            chk(f'full pooled {mod} holm tiny', h < float(f[6].strip('<')), True)
        else:
            chk(f'full pooled {mod} p', float(f[5]), round(p, 3))
            chk(f'full pooled {mod} holm', float(f[6]), round(h, 2))
chk('full pooled lists all 11 models', sorted(seen), sorted(M.values()))

# ---- 2. per-target table, full setting -----------------------------------------------------
print('\n2. Per-target table, `full` setting')
COLS = ['MultiPASS', 'pos_only', 'full_gbm', 'REVEL', 'AlphaMissense', 'CADD',
        'PolyPhen2', 'SIFT_neg', 'ESM2_650M_neg']
tbl = R.split('Per target, `full`')[1].split('## Result 2')[0]
got = rows_of(tbl.split('|---|---|---|---|---|---|---|---|---|---|---|')[1])
chk('6 per-target rows', len(got), 6)
for f in got:
    g = SHORT[f[0]]
    n = int(f[1].split()[0])
    npos, nneg = (int(x) for x in re.search(r'\((\d+)/(\d+)\)', f[1]).groups())
    chk(f'full {f[0]} n', n, int(cell('full', g, 'MultiPASS', 'n')))
    chk(f'full {f[0]} n_pos/n_neg', (npos, nneg),
        (int(cell('full', g, 'MultiPASS', 'n_pos')), int(cell('full', g, 'MultiPASS', 'n_neg'))))
    for mod, v in zip(COLS, f[2:11]):
        chk(f'full {f[0]} {mod}', float(v), round(float(cell('full', g, mod)), 3))

# ---- 3. pooled table, region setting -------------------------------------------------------
print('\n3. Pooled table, `region` setting')
blk = R.split('## Result 2')[1]
m = re.search(r'\((\d+) pos / (\d+) neg\)', blk)
chk('region pooled n_pos', int(m.group(1)), int(cell('region', 'POOLED', 'MultiPASS', 'n_pos')))
chk('region pooled n_neg', int(m.group(2)), int(cell('region', 'POOLED', 'MultiPASS', 'n_neg')))
seen = []
for f in rows_of(blk.split('|---|---|---|---|---|---|')[1]):
    mod = M[f[0]]; seen.append(mod)
    chk(f'region pooled {mod} auc', float(f[2]), round(float(cell('region', 'POOLED', mod)), 3))
    lo, hi = (float(x) for x in f[3].replace('–', '-').split('-'))
    chk(f'region pooled {mod} ci', (lo, hi),
        (round(float(cell('region', 'POOLED', mod, 'ci_lo')), 3),
         round(float(cell('region', 'POOLED', mod, 'ci_hi')), 3)))
    if mod != 'MultiPASS':
        p = float(cell('region', 'POOLED', mod, 'delong_p_vs_multipass'))
        h = float(cell('region', 'POOLED', mod, 'delong_p_holm'))
        if f[4].startswith('<'):
            chk(f'region pooled {mod} p tiny', p < float(f[4][1:]), True)
            chk(f'region pooled {mod} holm tiny', h < float(f[5].strip('<')), True)
        else:
            dec = len(f[4].split('.')[1])
            chk(f'region pooled {mod} p', float(f[4]), float(f'{p:.{dec}f}'))
            chk(f'region pooled {mod} holm', float(f[5]),
                float(f'{h:.{len(f[5].split(".")[1])}f}'))
chk('region pooled lists all 11 models', sorted(seen), sorted(M.values()))

# ---- 4. narrative claims -------------------------------------------------------------------
print('\n4. Narrative claims')
chk('140 DeLong tests in the CSV', int(E.delong_p_vs_multipass.notna().sum()), 140)
chk('log says 140', 'Holm correction across 140 DeLong tests' in LOG, True)
chk('report says 140 tests', '140 DeLong tests' in R, True)
BODY = R.split('## What is compared')[1]          # everything after the correction notice
chk('112 tests claim gone from the body', '112 DeLong' not in BODY and
    '112 tests' not in BODY, True)
chk('MultiPASS beats all 7 externals after Holm, full',
    all(float(cell('full', 'POOLED', M[k], 'delong_p_holm')) < 0.001
        for k in ['REVEL', 'AlphaMissense', 'CADD', 'PolyPhen-2', 'SIFT (−)',
                  'ESM-2 650M', 'ESM-2 150M']), True)
ext = [M[k] for k in ['REVEL', 'AlphaMissense', 'CADD', 'PolyPhen-2', 'SIFT (−)',
                      'ESM-2 650M', 'ESM-2 150M']]
chk('external AUPRC range, full pooled',
    (round(min(float(cell('full', 'POOLED', m, 'auprc')) for m in ext), 3),
     round(max(float(cell('full', 'POOLED', m, 'auprc')) for m in ext), 3)), (0.053, 0.097))
chk('report quotes that range', '0.239 against 0.053–0.097' in R, True)
chk('pos_only raw p, full pooled',
    round(float(cell('full', 'POOLED', 'pos_only', 'delong_p_vs_multipass')), 3), 0.027)
chk('pos_only Holm == 1.00, full pooled',
    round(float(cell('full', 'POOLED', 'pos_only', 'delong_p_holm')), 2), 1.0)
chk('report states raw p = 0.027 for pos_only', 'raw p = **0.027**' in R, True)
chk('0.058 gone from the body', '0.058' not in BODY, True)
chk('ESM650 is the only pooled region AUC below 0.5 among externals',
    sorted(m for m in ext if float(cell('region', 'POOLED', m)) < 0.5),
    ['ESM2_150M_neg', 'ESM2_650M_neg'])
chk('both ESM rows Holm < 0.001 in region pooled',
    all(float(cell('region', 'POOLED', m, 'delong_p_holm')) < 0.001
        for m in ['ESM2_650M_neg', 'ESM2_150M_neg']), True)
chk('ESM650 region ABL1 0.402 / Holm 0.044',
    (round(float(cell('region', SHORT['ABL1'], 'ESM2_650M_neg')), 3),
     round(float(cell('region', SHORT['ABL1'], 'ESM2_650M_neg', 'delong_p_holm')), 3)),
    (0.402, 0.044))
chk('ESM650 region ALK 0.439', round(float(cell('region', 'ALK', 'ESM2_650M_neg')), 3), 0.439)
chk('ESM650 region KIT 0.484', round(float(cell('region', 'KIT', 'ESM2_650M_neg')), 3), 0.484)
chk('region ABL1 AlphaMissense 0.471 / CADD 0.471',
    (round(float(cell('region', SHORT['ABL1'], 'AlphaMissense')), 3),
     round(float(cell('region', SHORT['ABL1'], 'CADD')), 3)), (0.471, 0.471))
chk('region ALK CADD 0.377 / AlphaMissense 0.461',
    (round(float(cell('region', 'ALK', 'CADD')), 3),
     round(float(cell('region', 'ALK', 'AlphaMissense')), 3)), (0.377, 0.461))
chk('model-size claim full 0.601 -> 0.614',
    (round(float(cell('full', 'POOLED', 'ESM2_150M_neg')), 3),
     round(float(cell('full', 'POOLED', 'ESM2_650M_neg')), 3)), (0.601, 0.614))
chk('model-size claim region 0.481 -> 0.486',
    (round(float(cell('region', 'POOLED', 'ESM2_150M_neg')), 3),
     round(float(cell('region', 'POOLED', 'ESM2_650M_neg')), 3)), (0.481, 0.486))
chk('windowing claim: ALK 0.730, EGFR 0.716, KIT 0.576, ABL1 0.547 (full, ESM650)',
    tuple(round(float(cell('full', SHORT[t], 'ESM2_650M_neg')), 3)
          for t in ['ALK', 'EGFR', 'KIT', 'ABL1']), (0.730, 0.716, 0.576, 0.547))

# ---- 5. macro-averaging caveat -------------------------------------------------------------
print('\n5. Macro-averaging caveat')
piv = E[(E.setting == 'full') & (E.target != 'POOLED')].pivot(index='target', columns='model',
                                                              values='auc')
mp = piv['MultiPASS']
for lab, mod, want in [('MultiPASS', 'MultiPASS', 0.858), ('AlphaMissense', 'AlphaMissense', 0.785),
                       ('REVEL', 'REVEL', 0.772), ('CADD', 'CADD', 0.708),
                       ('PolyPhen-2', 'PolyPhen2', 0.692), ('SIFT', 'SIFT_neg', 0.679),
                       ('ESM-2 650M', 'ESM2_650M_neg', 0.662),
                       ('ESM-2 150M', 'ESM2_150M_neg', 0.642)]:
    chk(f'macro {lab}', want, round(piv[mod].mean(), 3))
    chk(f'report quotes macro {lab}', f'{want:.3f}' in R, True)
for mod in ext:
    chk(f'MultiPASS ahead in 5 or 6 of 6 vs {mod}', int((mp > piv[mod]).sum()) in (5, 6), True)
chk('6/6 against REVEL and PolyPhen-2',
    (int((mp > piv['REVEL']).sum()), int((mp > piv['PolyPhen2']).sum())), (6, 6))
chk('the lost target is always NT5C2',
    {t for mod in ext for t in piv.index[mp <= piv[mod]]}, {'NT5C2'})
chk('vs pos_only: 5/6, median +0.0295',
    (int((mp > piv['pos_only']).sum()), round(float(np.median(mp - piv['pos_only'])), 4)),
    (5, 0.0295))
chk('vs full_gbm: 3/6, median -0.0005',
    (int((mp > piv['full_gbm']).sum()), round(float(np.median(mp - piv['full_gbm'])), 4)),
    (3, -0.0005))

# ---- 6. dropped-variant counts -------------------------------------------------------------
print('\n6. Dropped-variant counts')
want = [16, 10, 34, 21, 0, 3]
chk('report lists 16/10/34/21/0/3', '16 ABL1, 10 ALK, 34 EGFR, 21 KIT, 0 NT5C2, 3 SMO' in R, True)
chk('CSV n_dropped, full', [int(cell('full', SHORT[t], 'MultiPASS', 'n_dropped'))
                            for t in SHORT], want)
chk('pooled n == per-target sum, full',
    int(E[(E.setting == 'full') & (E.target != 'POOLED') & (E.model == 'MultiPASS')].n.sum()),
    int(cell('full', 'POOLED', 'MultiPASS', 'n')))

# ---- 7. AlphaMissense files: row counts, coverage, spot checks ------------------------------
print('\n7. AlphaMissense files (restaged from the author\'s machine 08-Oct)')
tot = 0
for t, a in ACC.items():
    d = pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv'); tot += len(d)
    chk(f'{t[:8]} rows divisible by 19', len(d) % 19, 0)
chk('119 396 rows over the six files', tot, 119396)
chk('report says 119 396', '119 396 rows over the six files' in R, True)
d = pd.read_csv(f'{AM}/AF-P00519-F1-aa-substitutions.csv')
chk('ABL1 21 470 = 1130 x 19', (len(d), len(d) // 19), (21470, 1130))
chk('report says 21 470 and 30 780', '**21 470**/19 = 1130' in R and '**30 780**/19 = 1620' in R,
    True)
chk('report no longer says 21471/30781', '21471' not in R and '30781' not in R, True)
amv = dict(zip(d.protein_variant, d.am_pathogenicity))
for v, want in [('T315I', 0.9991), ('Y253H', 0.9999), ('M351T', 0.9962), ('E255K', 0.99),
                ('G250E', 0.9919), ('F359V', 0.9666), ('M244V', 0.929), ('H396R', 0.7922)]:
    chk(f'AlphaMissense ABL1 {v}', want, round(float(amv[v]), 4))

e6 = pd.read_csv(f'{_SAVDR_ROOT}/level2/esm2_650M_scores.csv')
e1 = pd.read_csv(f'{_SAVDR_ROOT}/level2/esm2_150M_scores.csv')
chk('ESM 650M / 150M row counts', (len(e6), len(e1)), (119396, 119396))
ev = dict(zip(e6[e6.target == SHORT['ABL1']].protein_variant,
              e6[e6.target == SHORT['ABL1']].esm2_llr))
for v, want in [('T315I', -7.61), ('M351T', -8.22), ('E255K', -8.97)]:
    chk(f'ESM-2 650M ABL1 {v}', want, round(float(ev[v]), 2))

print('\n8. Coverage of our variants, full setting')
tot = cov = cove = 0
per = []
for t, a in ACC.items():
    amset = set(pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv').protein_variant)
    f = sorted(x for x in glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_*_s0/scores.csv')
               if 'region' not in x)[0]
    mut = pd.read_csv(f).aa.str.replace(r'^p\.', '', regex=True)
    es = set(e6[e6.target == t].protein_variant)
    tot += len(mut); cov += int(mut.isin(amset).sum()); cove += int(mut.isin(es).sum())
    per.append(len(mut))
chk('AlphaMissense coverage 4623/4623', (cov, tot), (4623, 4623))
chk('ESM-2 650M coverage 4623/4623', (cove, tot), (4623, 4623))
chk('per-target 712/1124/1214/960/142/471', per, [712, 1124, 1214, 960, 142, 471])
chk('report says per target 712/1124/1214/960/142/471',
    '712/1124/1214/960/142/471' in R, True)

# ---- 9. AlphaMissense pathogenic-share table -----------------------------------------------
print('\n9. AlphaMissense pathogenic-share table (region setting, threshold 0.564)')
share = R.split('region setting:')[1].split('AlphaMissense flags')[0]
rows = rows_of(share.split('|---|---|---|')[1])
chk('6 share rows', len(rows), 6)
for f in rows:
    t = SHORT[f[0]]
    a = ACC[t]
    amset = dict(zip(*[pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv')[c]
                       for c in ['protein_variant', 'am_pathogenicity']]))
    d = pd.read_csv(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_*_region_s0/scores.csv')[0])
    d['am'] = d.aa.str.replace(r'^p\.', '', regex=True).map(amset)
    chk(f'{f[0]} every region variant has an AM score', int(d.am.isna().sum()), 0)
    for y, col in [(1, 1), (0, 2)]:
        g = d[d.y == y]
        pct = round(float((g.am > 0.564).mean() * 100))
        chk(f'{f[0]} {"positives" if y else "negatives"} % pathogenic',
            int(f[col].replace('%', '').strip()), pct)


# ---- 10. the robustness table in ../verification/VERIFICATION.md ---------------------------
print('\n10. Robustness table in verification/VERIFICATION.md')
V = open(f'{_SAVDR_ROOT}/verification/VERIFICATION.md').read()
LAB = {'AlphaMissense': 'AlphaMissense', 'REVEL': 'REVEL', 'CADD': 'CADD',
       'PolyPhen-2': 'PolyPhen2', 'SIFT': 'SIFT_neg', 'ESM-2 650M': 'ESM2_650M_neg',
       'ESM-2 150M': 'ESM2_150M_neg', 'pos_only': 'pos_only', 'full_gbm': 'full_gbm'}
blk = V.split('| | pooled (variant-level)')[1]
n = 0
for line in blk.splitlines():
    if not line.strip().startswith('| vs '): continue
    f = [c.strip() for c in line.strip('|').split('|')]
    lab = f[0][3:].strip(); mod = LAB[lab]; n += 1
    mp_p, oth_p = (float(x) for x in f[1].split(' vs '))
    mp_m, oth_m = (float(x) for x in f[2].split(' vs '))
    chk(f'V pooled MultiPASS (row {lab})', mp_p,
        round(float(cell('full', 'POOLED', 'MultiPASS')), 3))
    chk(f'V pooled {mod}', oth_p, round(float(cell('full', 'POOLED', mod)), 3))
    chk(f'V macro MultiPASS (row {lab})', mp_m, round(piv['MultiPASS'].mean(), 3))
    chk(f'V macro {mod}', oth_m, round(piv[mod].mean(), 3))
    ahead = int(f[3].split('/')[0])
    chk(f'V ahead-in {mod}', ahead, int((mp > piv[mod]).sum()))
chk('V robustness table has 9 rows', n, 9)
chk('V no longer claims "all five external"', 'all five external' not in V, True)
chk('V 0.564 is AlphaMissense\'s own LPath boundary',
    (round(float(min(pd.read_csv(f'{AM}/AF-P00519-F1-aa-substitutions.csv')
      .query("am_class == \'LPath\'").am_pathogenicity)), 4),
     round(float(max(pd.read_csv(f'{AM}/AF-P00519-F1-aa-substitutions.csv')
      .query("am_class == \'Amb\'").am_pathogenicity)), 4)), (0.5640, 0.5637))
chk('V: share table identical from am_class and from the 0.564 threshold',
    all(round(float((g.am > 0.564).mean() * 100)) == round(float((g.cls == 'LPath').mean() * 100))
        for t, a in ACC.items()
        for g in [pd.read_csv(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_*_region_s0/scores.csv')[0]).assign(
            am=lambda d: d.aa.str.replace(r'^p\\.', '', regex=True).map(dict(zip(
                *[pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv')[c]
                  for c in ['protein_variant', 'am_pathogenicity']]))),
            cls=lambda d: d.aa.str.replace(r'^p\\.', '', regex=True).map(dict(zip(
                *[pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv')[c]
                  for c in ['protein_variant', 'am_class']]))))]
        for g in [g[g.y == 1], g[g.y == 0]]), True)
chk('ESM row count == 6284 positions x 19', 6284 * 19, 119396)
chk('sum of the six canonical lengths is 6284',
    sum(len(pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv')) // 19 for a in ACC.values()),
    6284)

# ---- 11. dbNSFP join sanity: the Spearman ranges (restaged 08-Oct) --------------------------
print('\n11. Join-sanity Spearman ranges (dbNSFP + AlphaMissense + ESM-2)')
DB = f'{_SAVDR_ROOT}/external/source-data/dbnsfp/dbnsfp_myvariant_6genes.tsv'
db = pd.read_csv(DB, sep='\t')
chk('dbNSFP substitution rows == 71 749', len(db), 71749)
chk('report says 71 749', '71 749 substitution rows' in R, True)
GENE = {'ABL1_ENST00000318560': 'ABL1', 'ALK': 'ALK', 'EGFR': 'EGFR', 'KIT': 'KIT',
        'NT5C2': 'NT5C2', 'SMO': 'SMO'}
def _sp(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    from scipy import stats as _st
    return float(_st.spearmanr(x[m], y[m]).statistic)
rng = {k: [] for k in ('am_revel', 'am_cadd', 'esm_am', 'esm_revel')}
for t, a in ACC.items():
    amm = dict(zip(*[pd.read_csv(f'{AM}/AF-{a}-F1-aa-substitutions.csv')[c]
                     for c in ['protein_variant', 'am_pathogenicity']]))
    sub = db[db.gene == GENE[t]].drop_duplicates('protein_variant').set_index('protein_variant')
    f = sorted(x for x in glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_*_s0/scores.csv')
               if 'region' not in x)[0]
    keys = [x[2:] for x in pd.read_csv(f).aa]
    g = sub.reindex(keys)
    amv = np.array([amm.get(k, np.nan) for k in keys])
    es = dict(zip(e6[e6.target == t].protein_variant, e6[e6.target == t].esm2_llr))
    esv = np.array([es.get(k, np.nan) for k in keys])
    rng['am_revel'].append(_sp(amv, g.revel.values))
    rng['am_cadd'].append(_sp(amv, g.cadd_phred.values))
    rng['esm_am'].append(_sp(esv, amv))
    rng['esm_revel'].append(_sp(esv, g.revel.values))
    chk(f'{t[:8]} dbNSFP coverage of our full-setting variants',
        int(g.revel.notna().sum()) > 0.9 * len(keys), True)
for k, want, txt in [('am_revel', (0.75, 0.85), '+0.75\nto +0.85'),
                     ('am_cadd', (0.64, 0.78), '+0.64 to +0.78'),
                     ('esm_am', (-0.79, -0.54), '\u22120.54 to \u22120.79'),
                     ('esm_revel', (-0.77, -0.46), '\u22120.46 to \u22120.77')]:
    lo, hi = round(min(rng[k]), 2), round(max(rng[k]), 2)
    chk(f'{k} range', (lo, hi), (round(min(want), 2), round(max(want), 2)))
chk('report quotes the AM~REVEL range', 'Spearman AlphaMissense~REVEL is' in R and
    '+0.75' in R and '+0.85' in R, True)
chk('report quotes the corrected ESM~AlphaMissense range',
    'AlphaMissense \u22120.54 to \u22120.79' in R, True)
chk('report no longer says \u22120.55 for that range', '\u22120.55 to \u22120.79' not in R, True)
print(f'\n{"=" * 100}\n{ok} passed, {fail} failed')
sys.exit(1 if fail else 0)
