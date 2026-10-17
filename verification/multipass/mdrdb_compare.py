"""MultiPASS against MdrDB -- the first QUANTITATIVE external benchmark in this project.

MdrDB (Yang, Ye, Qiu, Feng, Li, Hsieh, Allcock & Zhang, *Commun Chem* 6, 123, 2023;
quantum.tencent.com/mdrdb) collects binding-affinity changes upon mutation, ddG, for
protein-drug pairs. It is a database, not a predictor, so this is not a method comparison: it is
the first chance to ask whether our resistance score tracks a CONTINUOUS, externally measured
quantity rather than a curated binary label.

Why it matters that it is not COSMIC-derived. Every external check so far shares our labels'
ancestry: DRSP agreed with us on 439 of 440 variants, and CIViC overlaps COSMIC heavily. MdrDB's
ddG comes from GDSC/DepMap IC50 ratios, the Platinum affinity database, AIMMS, KinaseMD and a
tyrosine-kinase-inhibitor set -- none of them a resistance-annotation catalogue.

The data: `MdrDB_CoreSet_release_v1.0.2022.tsv`, taken from the project's own GitHub repository
(raw.githubusercontent.com/tencent-quantum-lab/MdrDB), which is ungated. The website's download
page requires registering a name, e-mail, telephone and institution, so the repository copy is
what is used here.

MdrDB's own convention: ddG > 1.36 kcal/mol is a ten-fold affinity drop and is called resistant.
Positive ddG therefore means weaker binding, i.e. resistance-ward.

Two readings are reported, and they answer different questions:
  * Spearman correlation between our score and ddG -- does the ranking track the measurement?
  * AUC against the ddG > 1.36 call -- does it separate the resistant tail?
Per source as well as pooled, because GDSC's ddG is derived from averaged IC50 values and is
much noisier than the Platinum/TKI measurements.
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

import glob, sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # midrank, delong, ci, baselines, first_record, ...

TSV = f'{_SAVDR_ROOT}/external/mdrdb/core.tsv'
TAB = f'{_SAVDR_ROOT}/inputs/variants'
FASTA = f'{_SAVDR_ROOT}/inputs/fasta/{{}}.faa'
LOTO = '/home/claude/loto11_L25_mna9_{}'
THR = 1.36                                # MdrDB's own resistance threshold

OURS = {'P00519': 'ABL1_ENST00000318560', 'Q9UM73': 'ALK', 'P00533': 'EGFR',
        'P10721': 'KIT', 'Q06187': 'BTK_ENST00000308731', 'P36888': 'FLT3',
        'P08581': 'MET', 'P03372': 'ESR1_ENST00000206249', 'Q99835': 'SMO',
        'P49902': 'NT5C2', 'O75874': 'IDH1_ENST00000415913'}
SHORT = {v: (k if False else v.split('_')[0]) for k, v in OURS.items()}
CVDIR = {'ABL1_ENST00000318560': 'fold_ABL1_ENST00000318560_17_9_s*', 'ALK': 'fold_ALK_33_4_s*',
         'EGFR': 'fold_EGFR_37_30_s*', 'KIT': 'fold_KIT_23_9_s*', 'SMO': 'fold_SMO_31_7_s*'}

d = pd.read_csv(TSV, sep='\t', low_memory=False)
print(f'MdrDB core set: {len(d)} samples, {d.UNIPROT_ID.nunique()} proteins, '
      f'{d.DRUG.nunique()} drugs, {d.MUTATION.nunique()} mutations')
s = d[(d.UNIPROT_ID.isin(OURS)) & (d.TYPE == 'Single Substitution')].copy()
s['target'] = s.UNIPROT_ID.map(OURS)
print(f'single substitutions on our targets: {len(s)} rows over '
      f'{s.target.nunique()} proteins, {s.DRUG.nunique()} drugs')
print(f'ddG present for all rows: {s["DDG.EXP"].notna().all()}; '
      f'{int((s["DDG.EXP"] > THR).sum())} rows above MdrDB\'s {THR} threshold')


def detect_offset(savs, wt_at):
    best = (0, -1, 0)
    for o in range(-60, 61):
        hit = tot = 0
        for v in savs:
            p = int(v[1:-1]) - o
            if p in wt_at:
                tot += 1
                hit += (wt_at[p] == v[0])
        if tot >= 5 and hit > best[1]:
            best = (o, hit, tot)
    return best


# MdrDB turns out to be internally inconsistent for ABL1: 107 of its 113 distinct ABL1
# substitutions are numbered on the 1130-residue canonical sequence (our ENST00000318560, i.e.
# offset 0) and six are numbered on the 1149-residue isoform b (+19). Rather than force a single
# global offset, each row is resolved against a SHORT, justified candidate list -- the isoforms
# that actually exist for that gene -- and accepted only when the wild-type residue matches.
# Allowing an arbitrary per-row offset would let one match in twenty through by chance.
CAND = {'ABL1_ENST00000318560': (0, 19)}


# ---- map onto our variant tables, with the usual wild-type / offset discipline --------------
rows, offs, lost = [], {}, []
for tgt, sub in s.groupby('target'):
    ours = pd.read_csv(f'{TAB}/{tgt}.csv')
    ours['mut'] = ours.aa.str.replace(r'^p\.', '', regex=True)
    tab = ours.set_index('mut')
    # the wild-type residue comes from the reference SEQUENCE, not from our variant table, so
    # that "numbering disagrees" and "we have no such variant" stay separate diagnoses
    seq = first_record(FASTA.format(tgt))
    wt_at = {i + 1: a for i, a in enumerate(seq)}
    savs = [str(x).strip() for x in sub.MUTATION]
    o, hit, tot = detect_offset(savs, wt_at)
    cands = CAND.get(tgt, (o,))
    offs[tgt.split('_')[0]] = dict(best_offset=o, wt_agree=f'{hit}/{tot}',
                                   frac=round(hit / tot, 4) if tot else None,
                                   candidates=str(cands))
    lo = pd.read_csv(LOTO.format(tgt) + '/scores.csv')
    lo['mut'] = lo.aa.str.replace(r'^p\.', '', regex=True)
    lopa = dict(zip(lo['mut'], lo.multipass_pa))
    cv = None
    if tgt in CVDIR:
        per = []
        for f in sorted(glob.glob(CVDIR[tgt] + '/scores.csv')):
            t = pd.read_csv(f)
            per.append(pd.Series(stats.rankdata(np.nan_to_num(t.multipass_pa.values, nan=0.0))
                                 / len(t), index=t.aa.str.replace(r'^p\.', '', regex=True)))
        cv = pd.concat(per, axis=1).mean(axis=1) if per else None
    for _, r in sub.iterrows():
        v = str(r.MUTATION).strip()
        use = None
        for c in cands:                      # first candidate whose wild-type residue matches
            q = int(v[1:-1]) - c
            if wt_at.get(q) == v[0]:
                use = c; break
        if use is None:
            lost.append((tgt, v, 'wt matches no isoform numbering')); continue
        p = int(v[1:-1]) - use
        ov = f'{v[0]}{p}{v[-1]}'
        if ov not in tab.index:
            lost.append((tgt, v, 'not in our table')); continue
        rows.append(dict(gene=tgt.split('_')[0], target=tgt, mdrdb_mut=v, mut=ov, pos=p,
                         offset_used=use,
                         drug=r.DRUG, dataset=r.DATASET, ddg=float(r['DDG.EXP']),
                         y_mdrdb=int(float(r['DDG.EXP']) > THR), y_ours=int(tab.loc[ov].y),
                         our_drugs=('' if pd.isna(tab.loc[ov].drugs) else tab.loc[ov].drugs),
                         loto_pa=lopa.get(ov, np.nan),
                         cv_rank=(float(cv.loc[ov]) if cv is not None and ov in cv.index
                                  else np.nan)))
M = pd.DataFrame(rows)
print('\nnumbering offsets (MdrDB position = our position + offset):')
print(pd.DataFrame(offs).T.to_string())
print(f'\nmatched {len(M)} of {len(s)} rows; dropped {len(lost)}')
if lost:
    ld = pd.DataFrame(lost, columns=['target', 'mut', 'why'])
    print(ld.groupby(['target', 'why']).size().to_string())
M.to_csv(f'{_SAVDR_ROOT}/multipass/mdrdb_overlap.csv', index=False)

print('\nmatched rows per gene:')
print(M.groupby('gene').agg(rows=('mut', 'size'), muts=('mut', 'nunique'),
                            drugs=('drug', 'nunique'),
                            res=('y_mdrdb', 'sum')).to_string())
print('\nby source:')
print(M.groupby('dataset').agg(rows=('mut', 'size'), res=('y_mdrdb', 'sum'),
                               ddg_med=('ddg', 'median')).round(3).to_string())

# one score per (gene, mutation): resistance to ANY drug = the largest ddG seen
agg = (M.groupby(['gene', 'target', 'mut'])
        .agg(ddg_max=('ddg', 'max'), ddg_med=('ddg', 'median'), n_drugs=('drug', 'nunique'),
             y_mdrdb=('y_mdrdb', 'max'), y_ours=('y_ours', 'first'),
             loto_pa=('loto_pa', 'first'), cv_rank=('cv_rank', 'first'))
        .reset_index())
print(f'\ncollapsed to {len(agg)} distinct (gene, mutation) pairs, '
      f'{int(agg.y_mdrdb.sum())} resistant by MdrDB\'s threshold')


def boot_spearman(x, y, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    n = len(x)
    out = [stats.spearmanr(x[i], y[i]).statistic
           for i in (rng.integers(0, n, n) for _ in range(B))]
    out = [v for v in out if not np.isnan(v)]
    return np.percentile(out, [2.5, 97.5]) if out else (np.nan, np.nan)


def panel(df, score, name, target_col='ddg_max'):
    df = df[df[score].notna()]
    if len(df) < 8: return None
    x = np.nan_to_num(df[score].values, nan=0.0)
    y = df[target_col].values
    rho = stats.spearmanr(x, y)
    lo, hi = boot_spearman(x, y)
    out = dict(panel=name, n=len(df), rho=round(float(rho.statistic), 3),
               rho_p=float(f'{rho.pvalue:.3g}'), rho_ci=f'{lo:.3f}..{hi:.3f}')
    yb = (y > THR).astype(int)
    if yb.sum() >= 3 and (1 - yb).sum() >= 3:
        a, b = ci(yb, x)
        out |= dict(n_res=int(yb.sum()), auc=round(roc_auc_score(yb, x), 3),
                    auc_ci=f'{a:.3f}..{b:.3f}',
                    auprc=round(average_precision_score(yb, x), 3))
    return out


print('\n' + '=' * 108)
print('Does our score track MdrDB\'s measured ddG?   (Spearman; AUC against ddG > 1.36)')
print('=' * 108)
res = []
for score, lab in [('cv_rank', 'within-protein CV'), ('loto_pa', 'LOTO')]:
    r = panel(agg, score, f'pooled, {lab}')
    if r: res.append(r)
    for dsname in ['TKI', 'AIMMS', 'GDSC']:
        keep = M[M.dataset == dsname][['gene', 'target', 'mut']].drop_duplicates()
        sub = agg.merge(keep, on=['gene', 'target', 'mut'])
        r = panel(sub, score, f'{dsname} only, {lab}')
        if r: res.append(r)
R = pd.DataFrame(res)
pd.set_option('display.width', 260)
print(R.to_string(index=False))
R.to_csv(f'{_SAVDR_ROOT}/multipass/mdrdb_compare.csv', index=False)

print('\nper gene (LOTO, pooled over drugs):')
pg = [panel(sub, 'loto_pa', g) for g, sub in agg.groupby('gene')]
PG = pd.DataFrame([x for x in pg if x])
print(PG.to_string(index=False))
PG.to_csv(f'{_SAVDR_ROOT}/multipass/mdrdb_pergene.csv', index=False)

print('\nagreement between MdrDB\'s threshold call and our COSMIC label:')
print(pd.crosstab(agg.y_mdrdb, agg.y_ours).to_string())
print(f'  our positives among these {len(agg)} mutations: {int(agg.y_ours.sum())}')
print(f'  MdrDB-resistant among them:                     {int(agg.y_mdrdb.sum())}')
print(f'  both:                                           {int(((agg.y_ours==1)&(agg.y_mdrdb==1)).sum())}')
ct = pd.crosstab(agg.y_mdrdb, agg.y_ours).reindex(index=[0,1], columns=[0,1]).fillna(0).values
n = ct.sum(); po = (ct[0,0] + ct[1,1]) / n
pe = ((ct[0].sum()*ct[:,0].sum()) + (ct[1].sum()*ct[:,1].sum())) / n**2
print(f'  raw agreement {po:.3f}; Cohen kappa {(po-pe)/(1-pe):.3f}')
print('  (contrast DRSP, whose labels agreed with ours on 439 of 440, kappa 0.995 --')
print('   MdrDB is a genuinely independent source, and it disagrees a great deal)')

# does the continuous score add anything over the binary label it was trained on?
yb = (agg.ddg_max > THR).astype(int).values
print('\nas predictors of MdrDB\'s ddG > 1.36 call, on the same 168 mutations:')
for nm, v in [('our COSMIC label (binary)', agg.y_ours.values.astype(float)),
              ('MultiPASS LOTO Pa', np.nan_to_num(agg.loto_pa.values, nan=0.0))]:
    a, b = ci(yb, v)
    print(f'  {nm:28s} AUC {roc_auc_score(yb, v):.3f}  CI [{a:.3f}, {b:.3f}]')
sub5 = agg[agg.cv_rank.notna()]
yb5 = (sub5.ddg_max > THR).astype(int).values
for nm, v in [('our COSMIC label (binary)', sub5.y_ours.values.astype(float)),
              ('MultiPASS within-protein CV', np.nan_to_num(sub5.cv_rank.values, nan=0.0))]:
    a, b = ci(yb5, v)
    print(f'  [5 CV genes, n={len(sub5)}] {nm:26s} AUC {roc_auc_score(yb5, v):.3f}  CI [{a:.3f}, {b:.3f}]')

print('\nmean ddG by our own label:')
print(agg.groupby('y_ours').ddg_max.agg(['count', 'mean', 'median']).round(3).to_string())
u = stats.mannwhitneyu(agg[agg.y_ours == 1].ddg_max, agg[agg.y_ours == 0].ddg_max,
                       alternative='greater')
print(f'  Mann-Whitney, our positives have higher ddG: U={u.statistic:.0f}, p={u.pvalue:.4g}')
