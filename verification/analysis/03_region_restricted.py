"""Robustness: restrict negatives to the protein segment spanned by resistance positions (+-10 aa).
Asks whether models beat chance once 'is it in the drug-binding region' is taken away."""
import sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
import importlib; B = importlib.import_module('02_baselines')
out = []
for t in [B.POOLED[int(i)] for i in sys.argv[1:]]:
    df, seq = B.build(t); y0 = df['any'].values
    lo, hi = df.loc[y0==1,'pos'].min()-10, df.loc[y0==1,'pos'].max()+10
    df = df[(df.pos>=lo)&(df.pos<=hi)].reset_index(drop=True); y = df['any'].values
    for model in ['pos_only','subst_lr','full_gbm']:
        X = B.X_for(model, df, seq)
        for scheme in ('random','grouped'):
            s = B.oof(model, X, y, df.pos.values, scheme).mean(0)
            (alo,ahi),(plo,phi) = B.ci(y, s)
            out.append(dict(target=t, span=f'{lo}-{hi}', n_pos=int(y.sum()), n_neg=int((1-y).sum()), model=model, cv=scheme,
                            auc=round(roc_auc_score(y,s),3), auc_lo=round(alo,3), auc_hi=round(ahi,3),
                            auprc=round(average_precision_score(y,s),3), prevalence=round(y.mean(),3)))
            print(out[-1], flush=True)
pd.DataFrame(out).to_csv('03_region_part'+'_'.join(sys.argv[1:])+'.csv', index=False)
