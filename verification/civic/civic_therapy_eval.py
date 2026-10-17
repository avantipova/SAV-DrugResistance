"""Therapy-matched CIViC evaluation -- the "Not done" caveat in REPORT_CIVIC.md.

That report matched CIViC's therapies against our drug coverage descriptively, but its scored
evaluation collapsed drug identity: a CIViC cetuximab call was tested against a score that
pooled every EGFR drug. Now that per-drug MultiPASS Pa exists for all 8 target x drug cells,
each CIViC call can be scored against the Pa for the drug CIViC actually names.

For a drug d:
  positives  CIViC resistance SAVs (Predictive / Resistance / Supports) whose therapy list
             contains d, present in our variant set
  negatives  our variants with COSMIC label 0 AND no CIViC resistance evidence at all
  score      Pa(d) from the per-drug multiclass run, rank-averaged over the 5 seeds
             (out-of-fold; nothing refit here)

Two positive sets, as in civic_eval.py:
  A  every drug-matched CIViC call
  B  only those our COSMIC labels call negative for d -- these were negatives during training,
     so ranking them highly is generalisation, not label recall

The comparison that matters: Pa(d), the drug-specific score, against Pa(Resistant), the binary
score, on the SAME drug-matched positives. If per-drug modelling adds anything, it shows here.
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

import glob, re, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # midrank, delong, ci, baselines, first_record

CFG = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
       'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
SHORT = {'ABL1_ENST00000318560': 'ABL1'}
# our COSMIC export labels NT5C2 with the drug CLASS; CIViC names the members
CLASS_MEMBERS = {'purine': {'mercaptopurine', 'thioguanine', '6-mercaptopurine', '6-thioguanine',
                            'purine'}}
MINPOS = 5

tm = pd.read_csv(f'{_SAVDR_ROOT}/civic/therapy_match.csv')
tm['civic_drug_set'] = tm.civic_drugs.fillna('').apply(lambda s: set(x for x in s.split('|') if x))
civic_variants = {(g, v) for g, v in zip(tm.gene, tm.variant)}


def matches(drug, civic_set):
    d = drug.strip().lower()
    if d in CLASS_MEMBERS:
        return bool(civic_set & CLASS_MEMBERS[d])
    return d in civic_set


def load_perdrug(target, L, MNA):
    """Rank-averaged Pa per activity, plus the binary Pa, keyed by variant name."""
    fs = sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/drug_{target}_{L}_{MNA}_s*/scores_drug.csv'))
    if not fs: return None, None, None
    per_act, base = {}, []
    d0 = None
    for f in fs:
        d = pd.read_csv(f)
        d0 = d if d0 is None else d0
        for c in [c for c in d.columns if c.startswith('pa::')]:
            act = c[4:]
            if act == 'Not resistant': continue
            v = stats.rankdata(np.nan_to_num(d[c].values, nan=0.0)) / len(d)
            per_act.setdefault(act, []).append(pd.Series(v, index=d.aa.values))
    bf = sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{target}_{L}_{MNA}_s*/scores.csv'))
    for f in bf:
        d = pd.read_csv(f)
        base.append(pd.Series(stats.rankdata(np.nan_to_num(d.multipass_pa.values, nan=0.0)) / len(d),
                              index=d.aa.values))
    act_avg = {a: pd.concat(v, axis=1).mean(1) for a, v in per_act.items()}
    bin_avg = pd.concat(base, axis=1).mean(1) if base else None
    return d0, act_avg, bin_avg


rows = []
for target, (L, MNA) in CFG.items():
    d0, act_avg, bin_avg = load_perdrug(target, L, MNA)
    if d0 is None:
        print(f'{target}: no per-drug run, skipped'); continue
    gene = SHORT.get(target, target)
    d0 = d0.copy(); d0['drugs'] = d0.drugs.fillna('')
    in_civic = np.array([(gene, a) in civic_variants for a in d0.aa])
    neg_pool = (d0.y.values == 0) & ~in_civic          # COSMIC-negative and CIViC-silent
    tmg = tm[tm.gene == gene].set_index('variant')

    for act, sc in sorted(act_avg.items()):
        if act.strip().lower() == 'tyrosine kinase inhibitor - ns':
            print(f'{gene}/{act}: class with no CIViC member list, cannot be therapy-matched')
            continue
        matched = np.array([a in tmg.index and matches(act, tmg.loc[a, 'civic_drug_set'])
                            for a in d0.aa])
        ours_d = np.array([act.lower() in [x.lower() for x in s.split(';')] if s else False
                           for s in d0.drugs])
        for setname, posmask in (('A_all_civic', matched),
                                 ('B_discordant_only', matched & ~ours_d)):
            keep = posmask | neg_pool
            y = posmask[keep].astype(int)
            if y.sum() < MINPOS:
                continue
            sub_aa = d0.aa.values[keep]
            s_drug = stats.rankdata(sc.loc[sub_aa].values) / keep.sum()
            s_bin = stats.rankdata(bin_avg.loc[sub_aa].values) / keep.sum()
            lo, hi = ci(y, s_drug)
            _, p = delong(y, s_drug, s_bin)
            rows.append(dict(target=gene, drug=act, label_set=setname,
                             n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                             auc_pa_drug=round(roc_auc_score(y, s_drug), 3),
                             ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                             auprc_pa_drug=round(average_precision_score(y, s_drug), 3),
                             auc_pa_binary=round(roc_auc_score(y, s_bin), 3),
                             delta_drug_minus_binary=round(roc_auc_score(y, s_drug)
                                                           - roc_auc_score(y, s_bin), 3),
                             delong_p_drug_vs_binary=round(p, 4)))

df = pd.DataFrame(rows)
# Holm across every test in this file
p = df.delong_p_drug_vs_binary.values
order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df['delong_p_holm'] = np.round(adj, 4)
df.to_csv(f'{_SAVDR_ROOT}/civic/civic_therapy_eval.csv', index=False)

print('\n' + '=' * 110)
print('THERAPY-MATCHED CIViC EVALUATION: drug-specific Pa vs the binary Pa, same positives')
print('=' * 110)
print(df.to_string(index=False))
for s in ('A_all_civic', 'B_discordant_only'):
    g = df[df.label_set == s]
    if g.empty: continue
    print(f'\n{s}: {len(g)} cells, mean delta (drug-specific minus binary) = '
          f'{g.delta_drug_minus_binary.mean():+.3f}, drug-specific ahead in '
          f'{int((g.delta_drug_minus_binary > 0).sum())} of {len(g)}')
print(f'\nHolm across {m} tests: significant = {int((df.delong_p_holm < 0.05).sum())}')
