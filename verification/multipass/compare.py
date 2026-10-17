"""Baselines on EXACTLY the variants MultiPASS scored, same single split (seed 0, 1 repeat)."""
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

import sys, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
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
d = pd.read_csv(sys.argv[1]); name = sys.argv[2]
seq = ''.join(l.strip() for l in open(f'{U}/fasta/{name}.faa') if not l.startswith('>'))
d['wt'] = [a.split('.')[1][0] for a in d.aa]; d['mut'] = [a.split('.')[1][-1] for a in d.aa]
y, g = d.y.values, d.pos.values
def feats(with_pos=True, win=7):
    X = {}
    for a in AAS: X['wt_'+a] = (d.wt == a).astype(float); X['mut_'+a] = (d.mut == a).astype(float)
    X['d_hyd'] = [KD[m]-KD[w] for w, m in zip(d.wt, d.mut)]
    X['d_vol'] = [VOL[m]-VOL[w] for w, m in zip(d.wt, d.mut)]
    X['d_chg'] = [CHG.get(m,0)-CHG.get(w,0) for w, m in zip(d.wt, d.mut)]
    if with_pos:
        X['pos_norm'] = d.pos / len(seq)
        for a in AAS: X['win_'+a] = [seq[max(0,p-1-win):p+win].count(a) for p in d.pos]
    return pd.DataFrame(X).values
M = {'pos_only': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.05, class_weight='balanced', random_state=0), (d.pos.values/len(seq)).reshape(-1,1)),
     'subst_lr': (lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight='balanced', max_iter=2000)), feats(False)),
     'full_gbm': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0), feats(True))}
out = []
for nm, (mk, X) in M.items():
    s = np.zeros(len(y))
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=0).split(X, y, g):
        m = mk(); m.fit(X[tr], y[tr]); s[te] = m.predict_proba(X[te])[:, 1]
    out.append((nm, round(roc_auc_score(y, s), 3), round(average_precision_score(y, s), 3)))
mp = d.multipass_pa.values
out.append(('MultiPASS', round(roc_auc_score(y, mp), 3), round(average_precision_score(y, mp), 3)))
print(f'{name}  n={len(y)}  pos={y.sum()}  neg={(1-y).sum()}  (single split, seed 0)')
for nm, a, p in out: print(f'  {nm:10s} AUC {a}  AUPRC {p}')
