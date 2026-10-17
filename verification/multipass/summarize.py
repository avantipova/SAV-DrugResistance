"""Pool MultiPASS + baseline scores over seeds; bootstrap CIs and a paired DeLong test.

For each target/setting and each seed, baselines are refit on exactly the variants MultiPASS
scored, with the same StratifiedGroupKFold(5, seed) split. AUC is reported as the mean over
seeds; the CI is a bootstrap over variants of the seed-averaged out-of-fold score; DeLong
compares MultiPASS with each baseline on those same averaged scores (paired, same variants).
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

import glob, re, sys, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

def first_record(path):
    """FIX: KIT (8 records), SMO (3) and FLT3 (11) are multi-record FASTA files. Reading them
    with a plain "skip the > lines" concatenates every isoform into one string, which lets a
    peptide window run past the end of the canonical sequence and into the next record's
    N-terminus instead of being rejected. Only the first record is the canonical sequence."""
    s, started = [], False
    for line in open(path):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)
U = f'{_SAVDR_ROOT}/external/source-data'
AAS = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
VOL = dict(A=88.6,R=173.4,N=114.1,D=111.1,C=108.5,Q=143.8,E=138.4,G=60.1,H=153.2,I=166.7,L=166.7,K=168.6,M=162.9,F=189.9,P=112.7,S=89.0,T=116.1,W=227.8,Y=193.6,V=140.0)
CHG = dict(D=-1,E=-1,K=1,R=1,H=0.1)

def midrank(x):
    s = np.argsort(x); z = x[s]; n = len(x); r = np.empty(n); i = 0
    while i < n:
        j = i
        while j < n - 1 and z[j+1] == z[i]: j += 1
        r[i:j+1] = 0.5 * (i + j) + 1; i = j + 1
    out = np.empty(n); out[s] = r; return out

def delong(y, s1, s2):
    """Paired DeLong test for two correlated ROC curves (DeLong 1988; Sun & Xu 2014)."""
    p = np.c_[s1, s2][y == 1].T; n = np.c_[s1, s2][y == 0].T
    m, nn, k = p.shape[1], n.shape[1], 2
    tx = np.array([midrank(p[r]) for r in range(k)])
    ty = np.array([midrank(n[r]) for r in range(k)])
    tz = np.array([midrank(np.r_[p[r], n[r]]) for r in range(k)])
    auc = (tz[:, :m].sum(1) / m - (m + 1) / 2) / nn
    v01 = (tz[:, :m] - tx) / nn; v10 = 1 - (tz[:, m:] - ty) / m
    S = np.cov(v01) / m + np.cov(v10) / nn
    L = np.array([[1, -1]]); d = L @ auc
    var = float(np.asarray(L @ S @ L.T).ravel()[0])
    z = float(np.asarray(d).ravel()[0]) / np.sqrt(var) if var > 0 else 0.0
    return auc, 2 * stats.norm.sf(abs(z))

def baselines(d, name, seed):
    seq = first_record(f'{U}/fasta/{name}.faa')
    wt = [a.split('.')[1][0] for a in d.aa]; mu = [a.split('.')[1][-1] for a in d.aa]
    y, g = d.y.values, d.pos.values
    def F(with_pos, win=7):
        X = {}
        for a in AAS:
            X['wt_'+a] = [float(w == a) for w in wt]; X['mut_'+a] = [float(m == a) for m in mu]
        X['d_hyd'] = [KD[m]-KD[w] for w, m in zip(wt, mu)]
        X['d_vol'] = [VOL[m]-VOL[w] for w, m in zip(wt, mu)]
        X['d_chg'] = [CHG.get(m,0)-CHG.get(w,0) for w, m in zip(wt, mu)]
        if with_pos:
            X['pos_norm'] = list(d.pos / len(seq))
            for a in AAS: X['win_'+a] = [seq[max(0,p-1-win):p+win].count(a) for p in d.pos]
        return pd.DataFrame(X).values
    M = {'pos_only': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.05, class_weight='balanced', random_state=0), (d.pos.values/len(seq)).reshape(-1,1)),
         'subst_lr': (lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight='balanced', max_iter=2000)), F(False)),
         'full_gbm': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0), F(True))}
    out = {}
    for nm, (mk, X) in M.items():
        s = np.zeros(len(y))
        for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(X, y, g):
            m = mk(); m.fit(X[tr], y[tr]); s[te] = m.predict_proba(X[te])[:, 1]
        out[nm] = s
    return out

def ci(y, s, B=600, seed=0):
    rng = np.random.default_rng(seed); pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]; a = []
    for _ in range(B):
        i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        a.append(roc_auc_score(y[i], s[i]))
    return np.percentile(a, [2.5, 97.5])

rows = []
groups = {}
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m: continue
    key = (m.group(1), int(m.group(2)), int(m.group(3)), bool(m.group(4)))
    groups.setdefault(key, []).append((int(m.group(5)), f))

for (name, L, MNA, region), runs in sorted(groups.items()):
    acc = {}
    per_seed = {}
    declined = []
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        mp = d.multipass_pa.values
        nmiss = int(np.isnan(mp).sum())
        if nmiss: declined.append((seed, nmiss, len(mp)))
        mp = np.nan_to_num(mp, nan=0.0)   # MultiPASS declared the class not predictable in that fold
        sc = baselines(d, name, seed); sc['MultiPASS'] = mp
        y = d.y.values
        for k, v in sc.items():
            acc.setdefault(k, []).append(stats.rankdata(v) / len(v))   # rank-average across seeds
            per_seed.setdefault(k, []).append(roc_auc_score(y, v))
    avg = {k: np.mean(v, 0) for k, v in acc.items()}
    tag = f'{name} L{L} MNA{MNA} {"region" if region else "full"}'
    print(f'\n{tag}: {y.sum()} pos / {(1-y).sum()} neg, {len(runs)} seeds'
          + (f'  [MultiPASS declined to predict the class in some folds: ' +
             ', '.join(f'seed {s} {n}/{t}' for s, n, t in declined) + ' -> scored 0]' if declined else ''))
    for k in ['pos_only', 'subject', 'subst_lr', 'full_gbm', 'MultiPASS']:
        if k not in avg: continue
        lo, hi = ci(y, avg[k])
        line = dict(target=name, setting='region' if region else 'full', model=k,
                    n_pos=int(y.sum()), n_neg=int((1-y).sum()), seeds=len(runs),
                    auc_mean_over_seeds=round(float(np.mean(per_seed[k])), 3),
                    auc_sd_over_seeds=round(float(np.std(per_seed[k])), 3),
                    auc_pooled=round(roc_auc_score(y, avg[k]), 3),
                    ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                    auprc=round(average_precision_score(y, avg[k]), 3))
        if k != 'MultiPASS':
            _, p = delong(y, avg['MultiPASS'], avg[k]); line['delong_p_vs_multipass'] = round(p, 4)
        rows.append(line)
        print('  ', {kk: vv for kk, vv in line.items() if kk not in ('target', 'setting', 'n_pos', 'n_neg')})
df = pd.DataFrame(rows)
# Holm-Bonferroni across every DeLong test in this table
mask = df.delong_p_vs_multipass.notna()
p = df.loc[mask, 'delong_p_vs_multipass'].values
order = np.argsort(p); m = len(p); adj = np.empty(m)
run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)
df.to_csv('summary.csv', index=False)
print('\nHolm-corrected DeLong p-values (', m, 'tests )')
print(df.loc[mask, ['target','setting','model','auc_pooled','delong_p_vs_multipass','delong_p_holm']].to_string(index=False))
