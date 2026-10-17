"""MultiPASS and the trained baselines against five external variant-effect predictors,
on exactly the variants of the main table.

Trained-here models (out-of-fold, rank-averaged over the 5 fold-assignment seeds, as in
summarize.py): MultiPASS, pos_only, subst_lr, full_gbm.
External, not refitted and never trained on drug-resistance labels:
  AlphaMissense  am_pathogenicity, per-protein CSV from AlphaFold DB (keyed to the canonical
                 UniProt we verified by CRC64), higher = more pathogenic
  REVEL          dbNSFP via MyVariant.info, higher = more deleterious
  CADD           phred, higher = more deleterious
  PolyPhen2      HDIV score, higher = more damaging
  SIFT           lower = more deleterious, so the score used is -SIFT
  ESM-2 650M     wt-marginal log-ratio log p(mut) - log p(wt) from esm2_t33_650M_UR50D over the
                 FULL canonical protein (windowed at 1022 residues with 256 overlap for the
                 three proteins longer than that); lower = more deleterious, so -LLR is used
  ESM-2 150M     the same from esm2_t30_150M_UR50D, kept alongside to show the size effect
dbNSFP reports per-nucleotide variants, so several genomic variants can give one amino-acid
substitution; scores are aggregated per substitution as max REVEL / max CADD / max PolyPhen2 /
min SIFT. dbNSFP's aa.pos is an array over isoforms and a row was emitted for every element,
so the join lands on our own numbering without an offset assumption.

Every model is scored on the SAME variants: the intersection where all eight have a value, so
each paired DeLong test uses identical variants. The number dropped is reported.
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
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # AAS/KD/VOL/CHG, midrank, delong, baselines, ci

ACC = {'ABL1_ENST00000318560': 'P00519', 'ALK': 'Q9UM73', 'EGFR': 'P00533',
       'KIT': 'P10721', 'NT5C2': 'P49902', 'SMO': 'Q99835'}
GENE = {'ABL1_ENST00000318560': 'ABL1', 'ALK': 'ALK', 'EGFR': 'EGFR',
        'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}

am = {}
for t, a in ACC.items():
    d = pd.read_csv(f'{U}/alphamissense/AF-{a}-F1-aa-substitutions.csv')
    am[t] = dict(zip(d.protein_variant, d.am_pathogenicity))
db = pd.read_csv(f'{U}/dbnsfp/dbnsfp_myvariant_6genes.tsv', sep='\t')
esm650 = pd.read_csv(f'{_SAVDR_ROOT}/level2/esm2_650M_scores.csv')
esm150 = pd.read_csv(f'{_SAVDR_ROOT}/level2/esm2_150M_scores.csv')
esmg650 = {t: dict(zip(g.protein_variant, g.esm2_llr)) for t, g in esm650.groupby('target')}
esmg150 = {t: dict(zip(g.protein_variant, g.esm2_llr)) for t, g in esm150.groupby('target')}
dbg = {}
for t, g in GENE.items():
    s = db[db.gene == g].drop_duplicates('protein_variant').set_index('protein_variant')
    dbg[t] = s

TRAINED = ['MultiPASS', 'pos_only', 'subst_lr', 'full_gbm']
EXTERNAL = ['AlphaMissense', 'REVEL', 'CADD', 'PolyPhen2', 'SIFT_neg', 'ESM2_650M_neg', 'ESM2_150M_neg']
ALL = TRAINED + EXTERNAL

groups = collections.defaultdict(list)
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv')):
    m = re.match(r'fold_(.+)_(\d+)_(\d+)(_region)?_s(\d+)/scores\.csv$', f)
    if not m: continue
    groups[(m.group(1), int(m.group(2)), int(m.group(3)), bool(m.group(4)))].append((int(m.group(5)), f))

rows = []
pooled = collections.defaultdict(list)
for (name, L, MNA, region), runs in sorted(groups.items()):
    if name not in ACC: continue
    acc = collections.defaultdict(list)
    for seed, f in sorted(runs):
        d = pd.read_csv(f)
        sc = baselines(d, name, seed)
        sc['MultiPASS'] = np.nan_to_num(d.multipass_pa.values, nan=0.0)
        for k, v in sc.items(): acc[k].append(stats.rankdata(v) / len(v))
    S = {k: np.mean(v, 0) for k, v in acc.items()}
    keys = [a[2:] for a in d.aa]                      # 'p.T315I' -> 'T315I'
    g = dbg[name].reindex(keys)
    S['AlphaMissense'] = np.array([am[name].get(k, np.nan) for k in keys])
    S['REVEL'] = g.revel.values
    S['CADD'] = g.cadd_phred.values
    S['PolyPhen2'] = g.polyphen2_hdiv.values
    S['SIFT_neg'] = -g.sift.values
    S['ESM2_650M_neg'] = np.array([-esmg650[name].get(k, np.nan) for k in keys])
    S['ESM2_150M_neg'] = np.array([-esmg150[name].get(k, np.nan) for k in keys])
    y = d.y.values
    keep = np.all([np.isfinite(S[k]) for k in ALL], axis=0)
    yy = y[keep]
    setting = 'region' if region else 'full'
    for k in ALL:
        s = stats.rankdata(S[k][keep]) / keep.sum()
        pooled[(setting, k)].append((yy, s))
        lo, hi = ci(yy, s)
        line = dict(setting=setting, target=name, model=k,
                    kind='trained_here' if k in TRAINED else 'external',
                    n=int(keep.sum()), n_dropped=int((~keep).sum()),
                    n_pos=int(yy.sum()), n_neg=int((1 - yy).sum()),
                    auc=round(roc_auc_score(yy, s), 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                    auprc=round(average_precision_score(yy, s), 3))
        if k != 'MultiPASS':
            ref = stats.rankdata(S['MultiPASS'][keep]) / keep.sum()
            _, p = delong(yy, ref, s); line['delong_p_vs_multipass'] = round(p, 4)
        rows.append(line)

for (setting, k), lst in pooled.items():
    ys = np.concatenate([a for a, _ in lst]); ss = np.concatenate([b for _, b in lst])
    ref = np.concatenate([b for _, b in pooled[(setting, 'MultiPASS')]])
    lo, hi = ci(ys, ss)
    line = dict(setting=setting, target='POOLED', model=k,
                kind='trained_here' if k in TRAINED else 'external',
                n=len(ys), n_dropped=0, n_pos=int(ys.sum()), n_neg=int((1 - ys).sum()),
                auc=round(roc_auc_score(ys, ss), 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                auprc=round(average_precision_score(ys, ss), 3))
    if k != 'MultiPASS':
        _, p = delong(ys, ref, ss); line['delong_p_vs_multipass'] = round(p, 4)
    rows.append(line)

df = pd.DataFrame(rows)
mask = df.delong_p_vs_multipass.notna()
p = df.loc[mask, 'delong_p_vs_multipass'].values
order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)
df.to_csv(f'{_SAVDR_ROOT}/level2/level2_eval.csv', index=False)

for setting in ('full', 'region'):
    print(f'\n################  {setting}  ################')
    t = df[df.setting == setting]
    for tgt in list(dict.fromkeys(t.target)):
        gg = t[t.target == tgt]; h = gg.iloc[0]
        print(f'\n{tgt}   n={h.n} ({h.n_pos} pos / {h.n_neg} neg), dropped {h.n_dropped} for missing external scores')
        for _, r in gg.sort_values('auc', ascending=False).iterrows():
            ph = '' if pd.isna(r.get('delong_p_holm')) else f'  p={r.delong_p_vs_multipass:<7.4f} Holm={r.delong_p_holm:.3f}'
            star = '*' if r.model == 'MultiPASS' else ' '
            print(f'  {star}{r.model:14s} {r.kind:12s} AUC={r.auc:.3f} CI[{r.ci_lo:.3f}-{r.ci_hi:.3f}] AUPRC={r.auprc:.3f}{ph}')
print(f'\nHolm correction across {m} DeLong tests')
