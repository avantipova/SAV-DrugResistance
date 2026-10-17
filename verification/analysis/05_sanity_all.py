"""04_sanity extended to every task. Usage: python3 05_sanity_all.py <idx...>
Tasks 0-13 = 02_baselines tasks (full negatives); 14-19 = region-restricted pooled targets (03).
Per task:
 leak        positions shared between train/test across 5 grouped repeats (must be 0)
 shuf_*      mean/max AUC over 5 label permutations, grouped CV, 2 repeats each (expect ~0.5)
 mw_diff     |sklearn AUC - manual Mann-Whitney AUC| for pos_only grouped (must be ~0)
 repro_*     rerun AUC minus saved AUC for pos_only/full_gbm grouped (must be 0)"""
import sys, numpy as np, pandas as pd, importlib
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
B = importlib.import_module('02_baselines')
TASKS = [(t, 'any', False) for t in B.POOLED] + [(t, d, False) for t, d in B.DRUG] + [(t, 'any', True) for t in B.POOLED]
S2 = pd.read_csv('02_baselines.csv'); S3 = pd.read_csv('03_region_restricted.csv')

def data(t, drug, region):
    df, seq = B.build(t)
    y = df['any'].values if drug == 'any' else df['drugs'].apply(lambda s: drug in s).astype(int).values
    if region:
        lo, hi = df.loc[y==1,'pos'].min()-10, df.loc[y==1,'pos'].max()+10
        keep = ((df.pos>=lo)&(df.pos<=hi)).values
        df = df[keep].reset_index(drop=True); y = y[keep]
    return df, seq, y

def saved(t, drug, region, model):
    S = S3 if region else S2.query('task==@drug')
    return S.query("target==@t and model==@model and cv=='grouped'").auc.item()

out = []
for i in map(int, sys.argv[1:]):
    t, drug, region = TASKS[i]
    df, seq, y = data(t, drug, region); g = df.pos.values
    leak = sum(len(set(g[tr]) & set(g[te])) for r in range(5)
               for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=r).split(df, y, g))
    row = dict(idx=i, target=t, task=drug, region=region, n_pos=int(y.sum()), n_neg=int((1-y).sum()), leak=leak)
    rng = np.random.default_rng(1)
    perms = [rng.permutation(y) for _ in range(5)]
    for model in ('pos_only', 'full_gbm'):
        X = B.X_for(model, df, seq)
        a = [roc_auc_score(ys, B.oof(model, X, ys, g, 'grouped', reps=2).mean(0)) for ys in perms]
        row[f'shuf_{model}_mean'] = round(float(np.mean(a)), 3); row[f'shuf_{model}_max'] = round(float(np.max(a)), 3)
        s = B.oof(model, X, y, g, 'grouped').mean(0)
        row[f'repro_{model}'] = round(round(roc_auc_score(y, s), 3) - saved(t, drug, region, model), 4)
        if model == 'pos_only':
            p, n = s[y==1], s[y==0]
            mw = ((p[:,None] > n[None,:]).sum() + 0.5*(p[:,None] == n[None,:]).sum()) / (len(p)*len(n))
            row['mw_diff'] = float(abs(mw - roc_auc_score(y, s)))
    out.append(row); print(row, flush=True)
pd.DataFrame(out).to_csv('05_sanity_part' + '_'.join(sys.argv[1:]) + '.csv', index=False)
