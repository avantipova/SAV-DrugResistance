"""Re-derive the selection-optimism section of REPORT_SELECTION_OPTIMISM.md.

Six targets have a nested estimate; three are provably not computable. Checks the
optimism values, the not-computable verdicts, both Spearman correlations with exact
permutation p-values, and ALK's optimism at full precision from its own checkpoints --
the last because the earlier 6e-4 tolerance was wider than the 1e-4 rounding artefact it
was meant to catch.
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

import pandas as pd, numpy as np, json, glob, os, re, sys, itertools, math
from scipy import stats
REP=f'{_SAVDR_ROOT}/reports/REPORT_SELECTION_OPTIMISM.md'
R=open(REP).read(); N=re.sub(r'\s+',' ',R)
ok=bad=0
def chk(n,c,d=''):
    global ok,bad
    if c: ok+=1
    else: bad+=1; print(f'  FAIL {n} {d}')
def inrep(s): chk(f'text {s!r}', s in N, '-- missing')

J={}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/multipass/nested_*.json')):
    b=os.path.basename(f)
    if b.endswith('_flat.json') or 'probe' in b or 'summary' in b: continue
    d=json.load(open(f))
    if 'optimism' in d and d.get('nested_auc') is not None:
        J[d['target'].split('_')[0]]=d
chk('six targets have a nested result', len(J)==6, sorted(J))
EXP={'ALK':(24,0.9288,0.9277,0.0010),'EGFR':(17,0.9233,0.9214,0.0019),
     'KIT':(24,0.8977,0.8917,0.0060),'ABL1':(101,0.8436,0.8242,0.0194),
     'SMO':(18,0.8572,0.8209,0.0363),'NT5C2':(18,0.7099,0.5762,0.1337)}
for t,(np_,mx,ne,op) in EXP.items():
    d=J[t]
    chk(f'{t} positives',d['n_pos']==np_,d['n_pos'])
    chk(f'{t} flat max',round(d['flat_max_auc'],4)==mx,d['flat_max_auc'])
    chk(f'{t} nested',round(d['nested_auc'],4)==ne,d['nested_auc'])
    chk(f'{t} optimism',round(d['optimism'],4)==op,d['optimism'])
    # The stored optimism is round(unrounded_flat - unrounded_nested, 4); the two AUC columns are
    # rounded independently, so their difference can be off by up to 1e-4 in each. Anything beyond
    # 2e-4 is a real disagreement, not rounding.
    gap = abs(d['flat_max_auc']-d['nested_auc']-d['optimism'])
    chk(f'{t} optimism within rounding of max - nested', gap <= 2.0001e-4, f'gap {gap:.6f}')
opt=np.array([EXP[t][3] for t in EXP])
chk('range 0.0010 to 0.1337',round(opt.min(),4)==0.0010 and round(opt.max(),4)==0.1337,(opt.min(),opt.max()))
chk('factor of 130',round(opt.max()/opt.min())==134 or abs(opt.max()/opt.min()-134)<5,
    round(opt.max()/opt.min(),1))
# picks table
PK={'ALK':(3,0.9275,0.9288,0.0013),'KIT':(3,0.8647,0.8977,0.0330),
    'ABL1':(3,0.8123,0.8436,0.0313),'SMO':(2,0.8029,0.8572,0.0543),
    'NT5C2':(4,0.5706,0.7099,0.1393)}
for t,(nd,lo,hi,w) in PK.items():
    d=J[t]; picks=list(d['nested_picked'].values()); tr=[d['flat_auc'][p] for p in picks]
    chk(f'{t} distinct configs',len(set(picks))==nd,len(set(picks)))
    chk(f'{t} picks true min',round(min(tr),4)==lo,round(min(tr),4))
    chk(f'{t} picks true max',round(max(tr),4)==hi,round(max(tr),4))
    chk(f'{t} picks span',round(max(tr)-min(tr),4)==w,round(max(tr)-min(tr),4))
chk('SMO has the most stable choice', PK['SMO'][0]==min(v[0] for v in PK.values()))
chk('ALK spread 0.145 exceeds NT5C2 0.147 claim is stated correctly',
    '0.145 against 0.147' in N)
# predictors, exact p
D=pd.read_csv(f'{_SAVDR_ROOT}/results/data/nested_mechanism.csv'); n=len(D)
def exact_p(x,y):
    rx=stats.rankdata(x); ry=stats.rankdata(y)
    obs=abs(stats.spearmanr(rx,ry).statistic); cnt=tot=0
    for perm in itertools.permutations(range(n)):
        tot+=1
        if abs(stats.spearmanr(rx,[ry[i] for i in perm]).statistic)>=obs-1e-12: cnt+=1
    return cnt/tot
WANT={'penalty':(1.000,0.0167),'width':(0.900,0.0833),'grid_mean':(-1.000,0.0167),
      'n_pos':(-0.632,0.2667),'spread':(0.500,0.4500),'top3':(0.500,0.4500),
      'n_distinct':(0.224,1.0000)}
for a,(rho,pv) in WANT.items():
    r=float(stats.spearmanr(D[a],D.optimism).statistic); p=exact_p(D[a].values,D.optimism.values)
    chk(f'rho {a}',round(r,3)==rho,f'{r:+.4f} != {rho}')
    chk(f'exact p {a}',round(p,4)==pv,f'{p:.4f} != {pv}')
chk('floor p at n=5 is 0.0167',round(2/math.factorial(5),4)==0.0167)
chk('Holm threshold for 7 tests is below the floor',round(0.05/7,4)==0.0071)
chk('report states nothing survives correction','**Nothing survives correction.**' in R)
chk('report withdraws the stability claim',
    'was wrong and is withdrawn' in N and 'where the picks are stable the optimism is small' in N)
chk('report withdraws the positives claim',
    'that optimism grows as the number of positives falls' in N)
for s_ in ['0.0010','0.0060','0.0194','0.0363','0.1337','0.0167','0.0833','0.2667',
           '0.4500','1.0000','tautological','0.0024','exact permutation p-values',
           '0.9275–0.9288','0.5706–0.7099','+0.22']:
    inrep(s_)
# ---------------------------------------------------------------------------------------------
# Added 2026-10-16: EGFR completed, MET returned a not-computable verdict.
# ---------------------------------------------------------------------------------------------
print('\n--- addendum: completion, the not-computable three, and ALK at full precision')

NOTCOMP = {'MET': (8, 855, 6), 'ESR1_ENST00000206249': (7, 387, 6),
           'IDH1_ENST00000415913': (6, 226, 6)}
for t, (np_, nn_, nf) in NOTCOMP.items():
    d = json.load(open(f'{_SAVDR_ROOT}/multipass/nested_{t}.json'))
    chk(f'{t} flagged not computable', d.get('nested_computable') is False, d.get('nested_computable'))
    chk(f'{t} has no nested AUC', d.get('nested_auc') is None, d.get('nested_auc'))
    chk(f'{t} positives', d['n_pos'] == np_, d['n_pos'])
    chk(f'{t} negatives', d['n_neg'] == nn_, d['n_neg'])
    chk(f'{t} all six configurations failed on the inner folds', d['inner_failed'] == nf, d['inner_failed'])
    chk(f'{t} flat grid still trains', len(d['flat_auc']) >= 4, len(d['flat_auc']))
    chk(f'{t} note names the IAP failure', 'average IAP 0.0000' in d['note'], d['note'][:60])
chk('IDH1 is the only one with configurations failing on the flat grid too',
    json.load(open(f'{_SAVDR_ROOT}/multipass/nested_IDH1_ENST00000415913.json'))['flat_not_trainable']
    == ['L17_MNA4', 'L25_MNA4'])
for t in ('MET', 'ESR1_ENST00000206249'):
    chk(f'{t} flat grid fully trainable',
        json.load(open(f'{_SAVDR_ROOT}/multipass/nested_{t}.json'))['flat_not_trainable'] == [])

# Spearman correlations over the six, with exact permutation p over all 720 orderings
import itertools
rows = sorted(((EXP[t][0], EXP[t][3]) for t in EXP), key=lambda r: r[1])
npos = [r[0] for r in rows]
optv = [r[1] for r in rows]
SPREAD = {}
for t in EXP:
    d = J[t]; fa = d['flat_auc']
    SPREAD[t] = round(max(fa.values()) - min(fa.values()), 4)
spread = [SPREAD[t] for t in sorted(EXP, key=lambda t: EXP[t][3])]

def exact_p(x, y):
    r0 = stats.spearmanr(x, y).statistic
    perms = list(itertools.permutations(range(len(y))))
    c = sum(1 for q in perms if abs(stats.spearmanr(x, [y[i] for i in q]).statistic) >= abs(r0) - 1e-12)
    return r0, c / len(perms), len(perms)

r_np, p_np, nperm = exact_p(npos, optv)
r_sp, p_sp, _ = exact_p(spread, optv)
chk('720 orderings at n=6', nperm == 720, nperm)
chk('Spearman optimism~positives is -0.12', round(r_np, 2) == -0.12, round(r_np, 4))
chk('its exact p is 0.84', round(p_np, 2) == 0.84, round(p_np, 4))
chk('Spearman optimism~grid spread is +0.43', round(r_sp, 2) == 0.43, round(r_sp, 4))
chk('its exact p is 0.42', round(p_sp, 2) == 0.42, round(p_sp, 4))
chk('the positives correlation weakened from -0.63 to -0.12', abs(r_np) < 0.63)
chk('EGFR and NT5C2 are the matched pair (17 vs 18 positives)',
    EXP['EGFR'][0] == 17 and EXP['NT5C2'][0] == 18)
chk('their optimism differs by about seventy-fold',
    68 < EXP['NT5C2'][3] / EXP['EGFR'][3] < 72, EXP['NT5C2'][3] / EXP['EGFR'][3])

# ALK at full precision, recomputed from its own per-fold checkpoints
import pandas as pd
from sklearn.metrics import roc_auc_score
dA = J['ALK']
fl = {}
L, M = dA['flat_best_config'].replace('L', '').split('_MNA')
for k in range(5):
    fl.update(json.load(open(f'{_SAVDR_ROOT}/inputs/nested_checkpoints/nested_ALK/flat_L{L}_M{M}_f{k}.json')))
ne = {}
for k, cfg in dA['nested_picked'].items():
    l, m = cfg.replace('L', '').split('_MNA')
    ne.update(json.load(open(f'{_SAVDR_ROOT}/inputs/nested_checkpoints/nested_ALK/flat_L{l}_M{m}_f{k}.json')))
com = sorted(set(fl) & set(ne))
yv = pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/ALK.csv').set_index('aa').loc[com, 'y'].values
fa_full = roc_auc_score(yv, [fl[a] for a in com])
na_full = roc_auc_score(yv, [ne[a] for a in com])
chk('ALK flat max full precision 0.92875', round(fa_full, 5) == 0.92875, fa_full)
chk('ALK nested full precision 0.9277462', round(na_full, 7) == 0.9277462, na_full)
chk('ALK optimism full precision 0.0010038', round(fa_full - na_full, 7) == 0.0010038, fa_full - na_full)
# Compare with a tolerance, not ==: 0.9288 - 0.9277 is 0.0011000000000000232 in binary float.
chk('so 0.0010 is right and 0.0011 is the two-sided rounding artefact',
    abs(round(fa_full - na_full, 4) - 0.0010) < 1e-12
    and abs((round(fa_full, 4) - round(na_full, 4)) - 0.0011) < 1e-12,
    f'{round(fa_full - na_full, 4)} vs {round(fa_full, 4) - round(na_full, 4)}')
chk('the artefact is the sum of the two rounding moves',
    abs((round(fa_full, 4) - fa_full) + (na_full - round(na_full, 4)) - 0.0000962) < 2e-6,
    (round(fa_full, 4) - fa_full, na_full - round(na_full, 4)))

# report text
for s_ in ['0.0019', '0.0010 to 0.1337', '-0.12', '0.84', 'not computable',
           'declined to train on 459 variants (2 resistant)', '0.0010038',
           'average IAP 0.0000', 'factor of seventy', 'the procedure has a floor',
           'Superseded 2026-10-16']:
    inrep(s_)
chk('report withdraws the positives-dependence claim in the headline',
    'the dependence on the number of resistant variants is gone' in N)
chk('report states the floor is bracketed but not located',
    'bracketed but not located' in N)

print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
