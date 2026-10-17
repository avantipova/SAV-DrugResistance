"""External-label evaluation: can each model rank CIViC-annotated resistance substitutions
above variants that neither COSMIC nor CIViC flags?

Scores are the out-of-fold scores already computed for the 'full' setting (all protein-wide
negatives), rank-averaged over the 5 fold-assignment seeds exactly as summarize.py does, so
no model is refitted here and nothing from a test fold ever entered its own score.

Two label sets, one shared negative pool:
  negatives (both)  variants in our set with COSMIC label 0 AND no CIViC resistance evidence
  A  all CIViC       positives = every CIViC resistance SAV (Supports) present in our set
  B  discordant      positives = only those CIViC SAVs our COSMIC labels call negative
B is the harder test: those variants were negatives in training, so a model that ranks them
highly is generalising rather than recalling its labels.
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

import glob, re, csv, collections, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # AAS/KD/VOL/CHG, midrank, delong, baselines, ci

CIVIC_F = U + '/civic/01-Oct-2026-AcceptedClinicalEvidenceSummaries.tsv'
GENE2TARGET = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'EGFR': 'EGFR',
               'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}
FUS = re.compile(r'\b(?:BCR|EML4|NPM1|LMNA|TPM3|ATIC|CLTC|STRN|KIF5B)::(ABL1|ALK)\b')
PAT = re.compile(r'\b(ABL1|ALK|EGFR|KIT|NT5C2|SMO)\s+([A-Z])(\d+)([A-Z])\b')

civic = collections.defaultdict(set)
for r in csv.DictReader(open(CIVIC_F, newline='', encoding='utf-8'), delimiter='\t'):
    if r['significance'] != 'Resistance' or r['evidence_type'] != 'Predictive': continue
    if r['evidence_direction'] != 'Supports': continue
    for g, wt, pos, mu in PAT.findall(FUS.sub(r'\1', r['molecular_profile'])):
        civic[GENE2TARGET[g]].add(f'p.{wt}{pos}{mu}')

# locate the 'full' runs per target
groups = collections.defaultdict(list)
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m or m.group(4): continue
    groups[(m.group(1), int(m.group(2)), int(m.group(3)))].append((int(m.group(5)), f))

MODELS = ['MultiPASS', 'pos_only', 'full_gbm', 'subst_lr']
rows, pooled = [], collections.defaultdict(lambda: collections.defaultdict(list))

for (name, L, MNA), runs in sorted(groups.items()):
    if name not in civic: continue
    acc = collections.defaultdict(list)
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        sc = baselines(d, name, seed)
        sc['MultiPASS'] = np.nan_to_num(d.multipass_pa.values, nan=0.0)
        for k, v in sc.items():
            acc[k].append(stats.rankdata(v) / len(v))
    avg = {k: np.mean(v, 0) for k, v in acc.items()}
    aa, y0 = d.aa.values, d.y.values
    C = np.array([a in civic[name] for a in aa])
    neg = (y0 == 0) & ~C                        # shared clean-negative pool
    sets = {'A_all_civic': C, 'B_discordant': C & (y0 == 0)}
    for tag, posmask in sets.items():
        keep = posmask | neg
        yy = posmask[keep].astype(int)
        if yy.sum() < 1: continue
        for k in MODELS:
            s = avg[k][keep]
            pooled[tag][k].append((yy, s))
            lo, hi = ci(yy, s)
            line = dict(label_set=tag, target=name, model=k, n_pos=int(yy.sum()),
                        n_neg=int((1 - yy).sum()),
                        auc=round(roc_auc_score(yy, s), 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                        auprc=round(average_precision_score(yy, s), 3))
            if k != 'MultiPASS':
                _, p = delong(yy, avg['MultiPASS'][keep], s)
                line['delong_p_vs_multipass'] = round(p, 4)
            rows.append(line)

# pooled over targets: within-target ranks are comparable, so concatenate them
for tag, per in pooled.items():
    ys = np.concatenate([y for y, _ in per['MultiPASS']])
    mps = np.concatenate([s for _, s in per['MultiPASS']])
    for k in MODELS:
        s = np.concatenate([s for _, s in per[k]])
        lo, hi = ci(ys, s)
        line = dict(label_set=tag, target='POOLED', model=k, n_pos=int(ys.sum()),
                    n_neg=int((1 - ys).sum()), auc=round(roc_auc_score(ys, s), 3),
                    ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                    auprc=round(average_precision_score(ys, s), 3))
        if k != 'MultiPASS':
            _, p = delong(ys, mps, s); line['delong_p_vs_multipass'] = round(p, 4)
        rows.append(line)

df = pd.DataFrame(rows)
mask = df.delong_p_vs_multipass.notna()
p = df.loc[mask, 'delong_p_vs_multipass'].values
order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)
df.to_csv(f'{_SAVDR_ROOT}/civic/civic_external_eval.csv', index=False)

for tag in ('A_all_civic', 'B_discordant'):
    print(f'\n===== {tag} =====')
    t = df[df.label_set == tag]
    for tgt in list(dict.fromkeys(t.target)):
        g = t[t.target == tgt]
        head = g.iloc[0]
        print(f'{tgt:24s} {head.n_pos:4d} pos / {head.n_neg:5d} neg')
        for _, r in g.iterrows():
            ph = '' if pd.isna(r.get('delong_p_holm')) else f'  p={r.delong_p_vs_multipass:<7.4f} pHolm={r.delong_p_holm:.3f}'
            print(f'    {r.model:10s} AUC={r.auc:.3f} CI[{r.ci_lo:.3f}-{r.ci_hi:.3f}] AUPRC={r.auprc:.3f}{ph}')
print(f'\nHolm correction applied across {m} DeLong tests in this file')
