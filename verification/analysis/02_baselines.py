"""Baseline classifiers on the same data MultiPASS saw (no MultiPASS needed).

Datasets
  pooled : per target, 'resistant to any drug' vs other COSMIC missense variants
  drug   : per target-drug pair (N>=15), one-vs-rest as in IAP (other-drug positives count as negatives)
Positives/negatives are built as in aas.ipynb cell 3: missense rows of csv/, labelled if the
AA change appears in mutations/. Only clean single-residue substitutions are used.

Models (fixed hyperparameters, no tuning -> no selection optimism)
  pos_only : gradient boosting on residue position only (hotspot memorisation)
  subst_lr : logistic regression on substitution features only (no position)
  full_gbm : gradient boosting on substitution features + position + local window composition
CV
  random  : repeated stratified 5-fold (same position may be in train and test)
  grouped : repeated stratified group 5-fold, groups = residue position
"""
import numpy as np, pandas as pd, warnings
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score
from common import *
warnings.filterwarnings('ignore')

AAS = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
VOL = dict(A=88.6,R=173.4,N=114.1,D=111.1,C=108.5,Q=143.8,E=138.4,G=60.1,H=153.2,I=166.7,L=166.7,K=168.6,M=162.9,F=189.9,P=112.7,S=89.0,T=116.1,W=227.8,Y=193.6,V=140.0)
CHG = dict(D=-1,E=-1,K=1,R=1,H=0.1)
_B = """A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
_rows = [l.split() for l in _B.splitlines()]; _ord = [r[0] for r in _rows]
BLOSUM = {(a, b): int(v) for r in _rows for b, v in zip(_ord, r[1:]) for a in [r[0]]}

def build(name):
    seq = read_fasta(name)[0][1]
    g = load_all_variants(name); g = g[g['Type'] == 'Substitution - Missense']
    r = load_resistance(name)
    rows = {}
    for m in g['AA Mutation'].dropna().unique():
        p = parse_missense(m)
        if p and p[1] <= len(seq) and seq[p[1]-1] == p[0]:
            rows[m] = dict(aa=m, wt=p[0], pos=p[1], mut=p[2], drugs=set())
    for m, d in zip(r['AA Mutation'], r['Drug']):
        if m in rows: rows[m]['drugs'].add(str(d))
    df = pd.DataFrame(rows.values())
    df['any'] = df['drugs'].apply(bool).astype(int)
    return df, seq

def features(df, seq, with_pos=True, with_window=True, win=7):
    X = {}
    for a in AAS:
        X['wt_'+a] = (df.wt == a).astype(float); X['mut_'+a] = (df.mut == a).astype(float)
    X['blosum'] = [BLOSUM[(w, m)] for w, m in zip(df.wt, df.mut)]
    X['d_hyd'] = [KD[m]-KD[w] for w, m in zip(df.wt, df.mut)]
    X['d_vol'] = [VOL[m]-VOL[w] for w, m in zip(df.wt, df.mut)]
    X['d_chg'] = [CHG.get(m,0)-CHG.get(w,0) for w, m in zip(df.wt, df.mut)]
    if with_pos: X['pos_norm'] = df.pos / len(seq)
    if with_window:
        for a in AAS:
            X['win_'+a] = [seq[max(0,p-1-win):p+win].count(a) for p in df.pos]
        X['win_hyd'] = [np.mean([KD[c] for c in seq[max(0,p-1-win):p+win] if c in KD]) for p in df.pos]
    return pd.DataFrame(X).values

MODELS = {
    'pos_only': lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.05, class_weight='balanced', random_state=0),
    'subst_lr': lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight='balanced', max_iter=2000)),
    'full_gbm': lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0),
}
def X_for(model, df, seq):
    if model == 'pos_only': return (df.pos.values / len(seq)).reshape(-1, 1)
    if model == 'subst_lr': return features(df, seq, with_pos=False, with_window=False)
    return features(df, seq)

def oof(model, X, y, groups, scheme, reps=5):
    P = np.zeros((reps, len(y)))
    for r in range(reps):
        if scheme == 'random':
            cv = StratifiedKFold(5, shuffle=True, random_state=r); splits = cv.split(X, y)
        else:
            cv = StratifiedGroupKFold(5, shuffle=True, random_state=r); splits = cv.split(X, y, groups)
        for tr, te in splits:
            m = MODELS[model](); m.fit(X[tr], y[tr]); P[r, te] = m.predict_proba(X[te])[:, 1]
    return P

def ci(y, s, B=1000, seed=0):
    rng = np.random.default_rng(seed); pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]; a = []; p = []
    for _ in range(B):
        i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        a.append(roc_auc_score(y[i], s[i])); p.append(average_precision_score(y[i], s[i]))
    return np.percentile(a, [2.5, 97.5]), np.percentile(p, [2.5, 97.5])

POOLED = ['ABL1_ENST00000318560', 'ALK', 'KIT', 'SMO', 'NT5C2', 'EGFR']
DRUG = [('ABL1','Imatinib'), ('ABL1_ENST00000318560','Imatinib'), ('ABL1_ENST00000318560','Tyrosine kinase inhibitor - NS'),
        ('ALK','Crizotinib'), ('EGFR','Osimertinib'), ('KIT','Imatinib'), ('NT5C2','Purine'), ('SMO','Vismodegib')]

if __name__ == '__main__':
    out = []
    import sys
    tasks = [(t, 'any') for t in POOLED] + DRUG
    sel = [int(i) for i in sys.argv[1:]] or range(len(tasks))
    tasks = [tasks[i] for i in sel]
    tag = '_'.join(map(str, sel)) if sys.argv[1:] else 'all'
    for t, drug in tasks:
        df, seq = build(t)
        y = df['any'].values if drug == 'any' else df['drugs'].apply(lambda s: drug in s).astype(int).values
        for model in MODELS:
            X = X_for(model, df, seq)
            for scheme in ('random', 'grouped'):
                P = oof(model, X, y, df.pos.values, scheme)
                s = P.mean(0)
                per_rep = [roc_auc_score(y, P[r]) for r in range(P.shape[0])]
                (alo, ahi), (plo, phi) = ci(y, s)
                out.append(dict(target=t, task=drug, n_pos=int(y.sum()), n_neg=int((1-y).sum()),
                                pos_positions=df.loc[y==1,'pos'].nunique(), model=model, cv=scheme,
                                auc=round(roc_auc_score(y, s),3), auc_lo=round(alo,3), auc_hi=round(ahi,3),
                                auc_rep_sd=round(float(np.std(per_rep)),3),
                                auprc=round(average_precision_score(y, s),3), auprc_lo=round(plo,3), auprc_hi=round(phi,3),
                                prevalence=round(y.mean(),4)))
                print(t, drug, model, scheme, out[-1]['auc'], flush=True)
    pd.DataFrame(out).to_csv(f'02_baselines_part{tag}.csv', index=False)
