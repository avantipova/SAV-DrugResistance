"""Shuffled-label null for morgan_lr, region-restricted, grouped CV. Usage: 07_morgan_null.py <task idx 6-11> <L> [n_perm]"""
import sys, numpy as np, pandas as pd, importlib
from sklearn.metrics import roc_auc_score
M = importlib.import_module('06_morgan'); B = M.B
i, L = int(sys.argv[1]), int(sys.argv[2]); n = int(sys.argv[3]) if len(sys.argv) > 3 else 20
t, region = M.TASKS[i]; df, seq, y = M.dataset(t, region); g = df.pos.values
X = np.log1p(np.vstack([M.fp(M.peptide(seq, p, m, L)) for p, m in zip(df.pos, df.mut)]))
obs = roc_auc_score(y, B.oof('morgan_lr', X, y, g, 'grouped').mean(0))
rng = np.random.default_rng(7)
null = [roc_auc_score(ys, B.oof('morgan_lr', X, ys, g, 'grouped', reps=2).mean(0)) for ys in (rng.permutation(y) for _ in range(n))]
p = (1 + sum(a >= obs for a in null)) / (n + 1)
print(dict(target=t, L=L, obs=round(obs,3), null_mean=round(np.mean(null),3), null_95=round(np.percentile(null,95),3), null_max=round(max(null),3), p=round(p,3), n_perm=n))
