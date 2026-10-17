"""Independent re-verification of every claim in the SAV-DrugResistance comparison.

Each check is written to fail loudly rather than to confirm. Run: python3 verify_all.py
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

import os, re, glob, collections, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats
RDLogger.DisableLog('rdApp.*')
U = f'{_SAVDR_ROOT}/external/source-data'
src = open('summarize.py').read(); exec(src.split('rows = []')[0])
CFG = {'ABL1_ENST00000318560': 17, 'ALK': 33, 'EGFR': 37, 'KIT': 23, 'NT5C2': 25, 'SMO': 31}
FAIL = []
def check(name, ok, detail=''):
    print(f'{"PASS" if ok else "FAIL"}  {name}' + (f'   {detail}' if detail else ''))
    if not ok: FAIL.append(name)

def extract(pep, mut, length):
    parts = mut.split('.'); pos = int(''.join(c for c in parts[1] if c.isdigit())) - 1
    if pos >= len(pep): return ''
    s = pos - length // 2; e = s + length
    if s < 0: return ''
    p = list(pep); p[pos] = parts[1][-1]; f = p[s:e]
    return ''.join(f) if len(f) == length else ''

def first_record(path):
    s, started = [], False
    for line in open(path):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)

def rebuild(name, L, region):
    seq = first_record(f'{U}/fasta/{name}.faa')
    g = pd.read_csv(f'{U}/csv/{name}.csv')
    g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
    res = set(pd.read_csv(f'{U}/mutations/{name}_mutations.csv')['AA Mutation'])
    rows = []
    for aa in g['AA Mutation'].dropna().unique():
        pep = extract(seq, aa, L)
        if not pep or Chem.MolFromSequence(pep) is None: continue
        rows.append(dict(aa=aa, pos=int(''.join(c for c in aa.split('.')[1] if c.isdigit())),
                         pep=pep, y=int(aa in res)))
    df = pd.DataFrame(rows)
    if region:
        lo, hi = df.loc[df.y == 1, 'pos'].min() - 10, df.loc[df.y == 1, 'pos'].max() + 10
        df = df[(df.pos >= lo) & (df.pos <= hi)].reset_index(drop=True)
    return df

def sdf_names(p):
    names, nxt = [], True
    for line in open(p):
        if nxt: names.append(line.strip()); nxt = False
        if line.startswith('$$$$'): nxt = True
    return [n for n in names if n]

print('\n== labels ==')
# the mutations file must contain only confirmed resistance rows
allyes = all(set(pd.read_csv(f)['Resistant Mutation'].astype(str)) == {'Yes'}
             for f in glob.glob(f'{U}/mutations/*_mutations.csv'))
check('every row of every *_mutations.csv is Resistant Mutation = Yes', allyes)

print('\n== fold integrity (the load-bearing claim) ==')
dirs = [d for d in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*_s*/')) if os.path.exists(d + 'scores.csv')]
leak = order = pos_leak = part = sdfmismatch = 0; nsdf = 0
for d in dirs:
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/$', d)
    name, L, region, seed = m.group(1), int(m.group(2)), bool(m.group(4)), int(m.group(5))
    sc = pd.read_csv(d + 'scores.csv'); df = rebuild(name, L, region)
    if list(sc.aa) != list(df.aa) or list(sc.y) != list(df.y): order += 1; continue
    y, grp = df.y.values, df.pos.values
    cover = np.zeros(len(df), int)
    for k, (tr, te) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(df, y, grp)):
        cover[te] += 1
        if set(grp[tr]) & set(grp[te]): pos_leak += 1
        f = d + f'train{k}.sdf'
        if not os.path.exists(f): continue
        nsdf += 1; tn = set(sdf_names(f))
        if tn & set(df.aa.values[te]): leak += 1
        if tn != set(df.aa.values[tr]): sdfmismatch += 1
    if not (cover == 1).all(): part += 1
check(f'scores.csv matches a fresh rebuild in all {len(dirs)} runs', order == 0, f'{order} mismatched')
check('no residue position on both sides of any split', pos_leak == 0, f'{pos_leak} leaks')
check(f'no test variant inside any train SDF ({nsdf} SDFs read)', leak == 0, f'{leak} leaks')
check('each train SDF equals its split train set exactly', sdfmismatch == 0, f'{sdfmismatch} differ')
check('test folds partition the variants exactly once', part == 0, f'{part} bad')

print('\n== no label leakage into the baselines ==')
d0 = pd.read_csv(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_ABL1_ENST00000318560_17_9_s0/scores.csv')
rng = np.random.default_rng(0); perm = collections.defaultdict(list)
for i in range(12):
    dd = d0.copy(); dd['y'] = rng.permutation(d0.y.values)
    for k, v in baselines(dd, 'ABL1_ENST00000318560', i).items():
        perm[k].append(roc_auc_score(dd.y.values, v))
for k, v in perm.items():
    check(f'permuted-label AUC near 0.5 for {k}', abs(np.mean(v) - .5) < .05, f'mean {np.mean(v):.3f}')

print('\n== statistics ==')
err = 0
for _ in range(200):
    n = rng.integers(40, 400); y = (rng.random(n) < rng.uniform(.1, .5)).astype(int)
    if y.sum() < 3 or (1 - y).sum() < 3: continue
    s1, s2 = rng.normal(y * 1.0, 1), rng.normal(y * .5, 1)
    a, _ = delong(y, s1, s2)
    err = max(err, abs(a[0] - roc_auc_score(y, s1)), abs(a[1] - roc_auc_score(y, s2)))
check('DeLong AUC equals sklearn AUC', err < 1e-12, f'max err {err:.1e}')
y = (rng.random(300) < .3).astype(int); s = rng.normal(y, 1)
check('DeLong p(s vs s) == 1', abs(delong(y, s, s)[1] - 1) < 1e-9)
ps = []
for _ in range(3000):
    yy = np.r_[np.ones(30), np.zeros(170)].astype(int)
    ps.append(delong(yy, rng.normal(yy, 1), rng.normal(yy, 1))[1])
ps = np.array(ps)
check('DeLong type-I error at 0.05 is nominal', abs((ps < .05).mean() - .05) < .015, f'{(ps<.05).mean():.3f}')
check('DeLong p-values uniform under the null', stats.kstest(ps, 'uniform').pvalue > .01,
      f'KS p={stats.kstest(ps,"uniform").pvalue:.3f}')
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); a = np.empty(m); r = 0.
    for rank, i in enumerate(o): r = max(r, (m - rank) * p[i]); a[i] = min(r, 1.)
    return a
try:
    from statsmodels.stats.multitest import multipletests
    e = max(np.abs(holm(p := rng.random(rng.integers(2, 60)) ** 2) - multipletests(p, method='holm')[1]).max()
            for _ in range(300))
    check('Holm matches statsmodels', e == 0, f'max diff {e:.1e}')
except ImportError:
    print('SKIP  statsmodels not installed')

print('\n== every CI brackets its own point estimate ==')
for f, c in (('summary.csv', 'auc_pooled'),
             (f'{_SAVDR_ROOT}/level2/level2_eval.csv', 'auc'),
             (f'{_SAVDR_ROOT}/civic/civic_external_eval.csv', 'auc')):
    x = pd.read_csv(f); v = ((x[c] < x.ci_lo - 1e-9) | (x[c] > x.ci_hi + 1e-9)).sum()
    check(f'{os.path.basename(f)}: CIs bracket the estimate', v == 0, f'{v}/{len(x)} violations')

print('\n== external joins ==')
ACC = {'ABL1_ENST00000318560': 'P00519', 'ALK': 'Q9UM73', 'EGFR': 'P00533',
       'KIT': 'P10721', 'NT5C2': 'P49902', 'SMO': 'Q99835'}
GENE = {k: k.split('_')[0] for k in ACC}
db = pd.read_csv(f'{U}/dbnsfp/dbnsfp_myvariant_6genes.tsv', sep='\t')
sc = pd.read_csv(f'{_SAVDR_ROOT}/multipass/scores_all.csv.gz')
for t, a in ACC.items():
    amd = pd.read_csv(f'{U}/alphamissense/AF-{a}-F1-aa-substitutions.csv')
    seq = first_record(f'{U}/fasta/{t}.faa')
    check(f'{t}: AlphaMissense CSV has len x 19 rows', len(amd) == len(seq) * 19,
          f'{len(amd)} vs {len(seq)*19}')
    m = dict(zip(amd.protein_variant, amd.am_pathogenicity))
    keys = [x[2:] for x in sc[sc.target == t].drop_duplicates('aa').aa]
    check(f'{t}: AlphaMissense covers every variant', all(k in m for k in keys),
          f'{sum(k in m for k in keys)}/{len(keys)}')
    g = db[db.gene == GENE[t]].drop_duplicates('protein_variant').set_index('protein_variant')
    A = np.array([m[k] for k in keys]); R = g.reindex(keys).revel.values
    ok = np.isfinite(A) & np.isfinite(R)
    rho = stats.spearmanr(A[ok], R[ok]).statistic
    check(f'{t}: AlphaMissense~REVEL correlated (a scrambled join would not be)', rho > .6, f'rho={rho:+.3f}')

print('\n== CIViC parsing, two independent routes ==')
cv = pd.read_csv(f'{U}/civic/01-Oct-2026-AcceptedClinicalEvidenceSummaries.tsv', sep='\t', low_memory=False)
r = cv[cv.significance == 'Resistance']
check('all Resistance evidence is evidence_type Predictive', set(r.evidence_type) == {'Predictive'})
sup = r[r.evidence_direction == 'Supports']
FUS = r'(?:BCR|EML4|NPM1|LMNA|TPM3|ATIC|CLTC|STRN|KIF5B)::'
GS = ['ABL1', 'ALK', 'EGFR', 'KIT', 'NT5C2', 'SMO']
A, B = set(), set()
for mp in sup.molecular_profile.astype(str):
    for g, wt, p, mu in re.findall(r'\b(' + '|'.join(GS) + r')\s+([A-Z])(\d+)([A-Z])\b',
                                   re.sub(r'\b' + FUS + r'(ABL1|ALK)\b', r'\1', mp)):
        A.add((g, f'p.{wt}{p}{mu}'))
    for cl in re.split(r'\s+AND\s+', mp):
        m2 = re.match(r'^(?:' + FUS + r')?(' + '|'.join(GS) + r')\s+([A-Z])(\d+)([A-Z])$', cl.strip())
        if m2: B.add((m2.group(1), f'p.{m2.group(2)}{m2.group(3)}{m2.group(4)}'))
check('two independent CIViC parsers agree', A == B, f'{len(A)} vs {len(B)}, symmetric diff {len(A^B)}')

print('\n' + ('ALL CHECKS PASSED' if not FAIL else f'{len(FAIL)} CHECK(S) FAILED: ' + '; '.join(FAIL)))
