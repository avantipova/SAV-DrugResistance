"""Leave-one-target-out: train on five proteins, test on the held-out sixth.
Answers 'did the model learn substitution chemistry, or memorise hot positions?'
Position features are meaningless across proteins, so pos_only is reported as a sanity floor.
Uses the same variant tables as the MultiPASS runs (scores.csv of seed 0)."""
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
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score
U = f'{_SAVDR_ROOT}/external/source-data'
AAS = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
VOL = dict(A=88.6,R=173.4,N=114.1,D=111.1,C=108.5,Q=143.8,E=138.4,G=60.1,H=153.2,I=166.7,L=166.7,K=168.6,M=162.9,F=189.9,P=112.7,S=89.0,T=116.1,W=227.8,Y=193.6,V=140.0)
CHG = dict(D=-1,E=-1,K=1,R=1,H=0.1)
TARGETS = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
           'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
frames = []
for t, (L, M) in TARGETS.items():
    cands = [f'fold_{t}_{L}_{M}_s0/scores.csv', f'fold_{t}_{L}_{M}_region_s0/scores.csv']
    f = next(c for c in cands if glob.glob(c))   # full-negative table when available
    d = pd.read_csv(f)
    print(f'{t}: {f}')
    seq = ''.join(l.strip() for l in open(f'{U}/fasta/{t}.faa') if not l.startswith('>'))
    d['target'] = t; d['seqlen'] = len(seq); d['seq'] = seq
    frames.append(d[['target', 'aa', 'pos', 'y', 'seqlen', 'seq']])
D = pd.concat(frames, ignore_index=True)
D['wt'] = [a.split('.')[1][0] for a in D.aa]; D['mut'] = [a.split('.')[1][-1] for a in D.aa]
def feats(d, with_pos, win=7):
    X = {}
    for a in AAS:
        X['wt_'+a] = [float(w == a) for w in d.wt]; X['mut_'+a] = [float(m == a) for m in d.mut]
    X['d_hyd'] = [KD[m]-KD[w] for w, m in zip(d.wt, d.mut)]
    X['d_vol'] = [VOL[m]-VOL[w] for w, m in zip(d.wt, d.mut)]
    X['d_chg'] = [CHG.get(m,0)-CHG.get(w,0) for w, m in zip(d.wt, d.mut)]
    if with_pos:
        X['pos_norm'] = list(d.pos / d.seqlen)
        for a in AAS: X['win_'+a] = [s[max(0,p-1-win):p+win].count(a) for p, s in zip(d.pos, d.seq)]
    return pd.DataFrame(X).values
M = {'pos_only': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.05, class_weight='balanced', random_state=0), False, True),
     'subst_lr': (lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight='balanced', max_iter=2000)), False, False),
     'full_gbm': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0), True, False)}
rows = []
for held in TARGETS:
    tr, te = D[D.target != held], D[D.target == held]
    for nm, (mk, wp, only_pos) in M.items():
        Xtr = (tr.pos.values/tr.seqlen.values).reshape(-1,1) if only_pos else feats(tr, wp)
        Xte = (te.pos.values/te.seqlen.values).reshape(-1,1) if only_pos else feats(te, wp)
        m = mk(); m.fit(Xtr, tr.y.values); s = m.predict_proba(Xte)[:, 1]
        rows.append(dict(held_out=held, model=nm, n_pos=int(te.y.sum()), n_neg=int((1-te.y).sum()),
                         auc=round(roc_auc_score(te.y, s), 3),
                         auprc=round(average_precision_score(te.y, s), 3),
                         prevalence=round(te.y.mean(), 3)))
r = pd.DataFrame(rows); r.to_csv('08_loto.csv', index=False)
print(r.pivot(index='held_out', columns='model', values='auc').to_string())
print('\nAUPRC'); print(r.pivot(index='held_out', columns='model', values='auprc').to_string())
print('\nprevalence'); print(r.drop_duplicates('held_out')[['held_out','n_pos','n_neg','prevalence']].to_string(index=False))
