"""Structure-descriptor models without MultiPASS.
Peptide = window of length L centred on the mutation with the MUTANT residue substituted
(same as extract_substring in aas.ipynb), built with RDKit Chem.MolFromSequence.
Descriptors = Morgan count fingerprint (radius 2, 2048 bits) — an open analogue of atom-neighbourhood
MNA descriptors, NOT MNA itself.
Lengths 11/17/25 are fixed in advance and all reported (no selection over lengths).
All models use the same variant subset: positions where a 25-residue window fits.
Models: morgan_nb (Bernoulli naive Bayes, closest in spirit to PASS), morgan_lr (L2 logistic),
pos_only (reference, same subset).
Usage: python3 06_morgan.py <task idx...>  (0-5 pooled full negatives, 6-11 pooled region-restricted)"""
import sys, numpy as np, pandas as pd, importlib
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.naive_bayes import BernoulliNB
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
RDLogger.DisableLog('rdApp.*')
B = importlib.import_module('02_baselines')
import os
LENS, LMAX = tuple(int(x) for x in os.environ.get('LENS','11,17,25').split(',')), 25
WITH_POS = os.environ.get('WITH_POS','1') == '1'
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
TASKS = [(t, False) for t in B.POOLED] + [(t, True) for t in B.POOLED]
MODELS = dict(B.MODELS)
MODELS['morgan_nb'] = lambda: BernoulliNB(alpha=1.0)
MODELS['morgan_lr'] = lambda: LogisticRegression(C=0.05, class_weight='balanced', max_iter=3000)
B.MODELS.update(MODELS)   # B.oof looks models up here

def peptide(seq, pos, mut, L):
    s = pos - 1 - L // 2; e = s + L
    if s < 0 or e > len(seq): return None
    p = list(seq[s:e]); p[L // 2] = mut; return ''.join(p)

def fp(pep):
    m = Chem.MolFromSequence(pep)
    return None if m is None else GEN.GetCountFingerprintAsNumPy(m).astype(np.float32)

def dataset(t, region):
    df, seq = B.build(t)
    ok = [(p - 1 - LMAX // 2 >= 0) and (p - 1 - LMAX // 2 + LMAX <= len(seq)) for p in df.pos]
    df = df[ok].reset_index(drop=True); y = df['any'].values
    if region:
        lo, hi = df.loc[y==1,'pos'].min()-10, df.loc[y==1,'pos'].max()+10
        df = df[(df.pos>=lo)&(df.pos<=hi)].reset_index(drop=True); y = df['any'].values
    return df, seq, y

out = []
for i in map(int, sys.argv[1:]):
    t, region = TASKS[i]
    df, seq, y = dataset(t, region); g = df.pos.values
    feats = {'pos_only': ('-', B.X_for('pos_only', df, seq))} if WITH_POS else {}
    for L in LENS:
        X = np.vstack([fp(peptide(seq, p, m, L)) for p, m in zip(df.pos, df.mut)])
        feats[f'L{L}'] = (L, X)
    for key, (L, X) in feats.items():
        models = ['pos_only'] if key == 'pos_only' else ['morgan_nb', 'morgan_lr']
        for model in models:
            Xm = (X > 0).astype(np.float32) if model == 'morgan_nb' else (np.log1p(X) if model == 'morgan_lr' else X)
            for scheme in ('random', 'grouped'):
                s = B.oof(model, Xm, y, g, scheme).mean(0)
                (alo, ahi), _ = B.ci(y, s, B=500)
                out.append(dict(target=t, region=region, n_pos=int(y.sum()), n_neg=int((1-y).sum()),
                                model=model, length=L, cv=scheme, auc=round(roc_auc_score(y, s), 3),
                                auc_lo=round(alo, 3), auc_hi=round(ahi, 3), auprc=round(average_precision_score(y, s), 3)))
                print(t, 'region' if region else 'full', model, L, scheme, out[-1]['auc'], flush=True)
pd.DataFrame(out).to_csv('06_morgan_part' + '_'.join(sys.argv[1:]) + '_L' + '-'.join(map(str,LENS)) + ('' if WITH_POS else '_nopos') + '.csv', index=False)
