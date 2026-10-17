"""Sanity checks for 02/03: (a) grouped folds never share a position between train and test;
(b) shuffled labels give AUC ~0.5 through the same pipeline; (c) sklearn AUC == manual Mann-Whitney AUC;
(d) rerun of one cell reproduces the saved number."""
import numpy as np, pandas as pd, importlib
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
B = importlib.import_module('02_baselines')
df, seq = B.build('ALK'); y = df['any'].values; g = df.pos.values
# (a)
leak = 0
for r in range(5):
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=r).split(df, y, g):
        leak += len(set(g[tr]) & set(g[te]))
print('(a) positions shared train/test in grouped CV:', leak)
# (b)
rng = np.random.default_rng(1); res = []
for k in range(5):
    ys = rng.permutation(y)
    for model in ('pos_only', 'full_gbm'):
        s = B.oof(model, B.X_for(model, df, seq), ys, g, 'grouped', reps=2).mean(0)
        res.append((model, roc_auc_score(ys, s)))
r = pd.DataFrame(res, columns=['model', 'auc']).groupby('model').auc.agg(['mean', 'min', 'max']).round(3)
print('(b) shuffled-label AUC (ALK, grouped, 5 permutations):'); print(r)
# (c)+(d)
s = B.oof('pos_only', B.X_for('pos_only', df, seq), y, g, 'grouped').mean(0)
p, n = s[y == 1], s[y == 0]
mw = ((p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum()) / (len(p) * len(n))
saved = pd.read_csv('02_baselines.csv').query("target=='ALK' and task=='any' and model=='pos_only' and cv=='grouped'").auc.item()
print('(c) sklearn', round(roc_auc_score(y, s), 4), 'manual MW', round(mw, 4), '| (d) saved', saved)
