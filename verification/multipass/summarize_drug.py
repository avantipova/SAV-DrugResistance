"""Per-drug comparison, in the shape of the paper's own Table 1.

Table 1 reports one figure per drug; the binary run (summarize.py) collapses every drug into
"Resistant". This rebuilds the comparison per drug, on the same position-grouped folds.

For a drug d:
    positives = variants annotated resistant to d
    negatives = variants with no resistance annotation at all
    excluded  = variants resistant only to a DIFFERENT drug (their status for d is unknown,
                so counting them as negatives would be wrong and counting them as positives
                would be worse)
MultiPASS's score is Pa(d) from the multiclass base of the same fold (one base predicts every
drug at once, exactly as the project's own SDF builder sets it up). Baselines are refit on the
same restricted variant set with the same seed and the same position grouping.
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

import glob, re, collections, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats
src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # AAS/KD/VOL/CHG, midrank, delong, baselines, ci, first_record

MINPOS = 10
runs = collections.defaultdict(list)
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/drug_*/scores_drug.csv')):
    m = re.match(r'drug_(.+)_(\d+)_(\d+)_s(\d+)/scores_drug\.csv$', f)
    runs[(m.group(1), int(m.group(2)), int(m.group(3)))].append((int(m.group(4)), f))

rows = []
for (name, L, MNA), rr in sorted(runs.items()):
    d0 = pd.read_csv(sorted(rr)[0][1])
    d0['drugs'] = d0.drugs.fillna('')
    drugs = sorted({x for s in d0.drugs for x in (s.split(';') if s else [])})
    for drug in drugs:
        y_d = np.array([int(drug in (s.split(';') if s else [])) for s in d0.drugs])
        neg = (d0.y.values == 0)
        keep = (y_d == 1) | neg
        if y_d[keep].sum() < MINPOS: continue
        sub = d0[keep].reset_index(drop=True)
        yy = y_d[keep]
        acc = collections.defaultdict(list); declined = []
        for seed, f in sorted(rr):
            d = pd.read_csv(f); d['drugs'] = d.drugs.fillna('')
            s = d.loc[keep, 'pa::' + drug].values if 'pa::' + drug in d.columns else np.full(keep.sum(), np.nan)
            nmiss = int(np.isnan(s).sum())
            if nmiss: declined.append((seed, nmiss))
            acc['MultiPASS'].append(stats.rankdata(np.nan_to_num(s, nan=0.0)) / len(s))
            b = baselines(sub.assign(y=yy), name, seed)
            for k, v in b.items(): acc[k].append(stats.rankdata(v) / len(v))
        avg = {k: np.mean(v, 0) for k, v in acc.items()}
        for k in ['MultiPASS', 'pos_only', 'subst_lr', 'full_gbm']:
            lo, hi = ci(yy, avg[k])
            line = dict(target=name, drug=drug, model=k, n_pos=int(yy.sum()), n_neg=int((1 - yy).sum()),
                        n_excluded=int(((y_d == 0) & ~neg).sum()), seeds=len(rr),
                        auc=round(roc_auc_score(yy, avg[k]), 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                        auprc=round(average_precision_score(yy, avg[k]), 3),
                        multipass_declined='; '.join(f'seed {s} {n}' for s, n in declined) if k == 'MultiPASS' else '')
            if k != 'MultiPASS':
                _, p = delong(yy, avg['MultiPASS'], avg[k]); line['delong_p_vs_multipass'] = round(p, 4)
            rows.append(line)

df = pd.DataFrame(rows)
mask = df.delong_p_vs_multipass.notna()
p = df.loc[mask, 'delong_p_vs_multipass'].values
o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, i in enumerate(o):
    run = max(run, (m - rank) * p[i]); adj[i] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)
df.to_csv(f'{_SAVDR_ROOT}/multipass/per_drug.csv', index=False)

for (t, d), g in df.groupby(['target', 'drug']):
    h = g.iloc[0]
    print(f'\n{t}  /  {d}   n+={h.n_pos} n-={h.n_neg} (excluded {h.n_excluded} resistant only to other drugs)')
    if h.multipass_declined: print(f'   MultiPASS gave no Pa for this drug in some folds: {h.multipass_declined}')
    for _, r in g.sort_values('auc', ascending=False).iterrows():
        ph = '' if pd.isna(r.get('delong_p_holm')) else f'  p={r.delong_p_vs_multipass:<7.4f} Holm={r.delong_p_holm:.3f}'
        print(f'   {"*" if r.model=="MultiPASS" else " "}{r.model:10s} AUC={r.auc:.3f} CI[{r.ci_lo:.3f}-{r.ci_hi:.3f}] AUPRC={r.auprc:.3f}{ph}')
print(f'\nHolm across {m} DeLong tests')
