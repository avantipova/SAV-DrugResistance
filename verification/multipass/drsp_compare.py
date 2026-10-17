"""MultiPASS vs DRSP (Lin, Y.-F. et al., Pharmaceuticals 2022, 15, 136) on DRSP's own variant set.

DRSP is the second resistance-specific comparator the manuscript cites that has no web server
and no code -- it is reference 22, the one dismissed in the Comparison section with "limited by
the lack of publicly accessible web interfaces". It cannot be run. But its **labelled dataset is
published** as supplementary material (`Dataset.xlsx` inside `pharmaceuticals-15-00136-s001.zip`),
so our model can be scored on exactly their variants and the result put beside their published
figures. That is a real head-to-head without their server.

What DRSP is: an ML classifier of single amino acid variants as drug-resistant or not, built on
45 features (drug-SAV interaction energies from iGEMDOCK, weighted contact number, B-factor and
DSSP structure terms, BLOSUM62/PAM250/PSSM conservation) computed on a drug-protein complex,
either a crystal structure or a docking model.

Their data design, read from the file rather than assumed:
  * Training set 610 SAVs over 11 genes; positives carry the drug name, negatives are labelled
    `NA_deleterious` -- deleterious-but-not-resistance variants.
  * Testing set 115 SAVs over BRAF, MAP2K2, ROS1 -- none of which we model, so their test set is
    unusable for us and only the training set overlaps our panel.

**The asymmetry runs in DRSP's favour and must be stated.** These 610 variants are the set they
fitted and cross-validated on; our scores are out-of-fold under position-grouped CV and were
produced by a model that never saw a DRSP label. A tie therefore favours us.

**Their negative class is not ours.** `NA_deleterious` means "damaging but not resistance", while
our negatives mean "no resistance annotation in COSMIC". Using DRSP's labels on DRSP's variants
makes this an independent label source, in the same role CIViC plays, and the label disagreement
with our own COSMIC labels is reported separately below.

Two score panels, because our coverage differs by target:
  A. LOTO (leave-one-target-out, L=25 / MNA 9): the held-out protein was never in training. All
     nine shared genes have this. The most conservative reading.
  B. Within-protein CV (`full` setting, per-target best configuration, rank-averaged over the 5
     fold-assignment seeds): only ABL1, ALK, EGFR, KIT, SMO have it.
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

import glob, json, sys, numpy as np, pandas as pd
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef,
                             f1_score, accuracy_score, precision_score, recall_score)
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # midrank, delong, ci, baselines, first_record, ...

XLSX = f'{_SAVDR_ROOT}/external/drsp/x/pharmaceuticals-1557619-supplementary/Dataset.xlsx'
TAB = f'{_SAVDR_ROOT}/inputs/variants'
LOTO = '/home/claude/loto11_L25_mna9_{}'

# DRSP gene -> our target key in tab11_L25 / loto11
GENE2TGT = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'BTK': 'BTK_ENST00000308731',
            'EGFR': 'EGFR', 'ESR1': 'ESR1_ENST00000206249', 'FLT3': 'FLT3', 'KIT': 'KIT',
            'MET': 'MET', 'SMO': 'SMO'}
# our own within-protein CV runs (full setting), target -> fold dir glob
CVDIR = {'ABL1_ENST00000318560': 'fold_ABL1_ENST00000318560_17_9_s*',
         'ALK': 'fold_ALK_33_4_s*', 'EGFR': 'fold_EGFR_37_30_s*',
         'KIT': 'fold_KIT_23_9_s*', 'SMO': 'fold_SMO_31_7_s*'}

tr = pd.read_excel(XLSX, sheet_name='Training set')
te = pd.read_excel(XLSX, sheet_name='Testing set')
print(f'DRSP training set {len(tr)} SAVs over {tr.Gene.nunique()} genes; '
      f'testing set {len(te)} over {sorted(te.Gene.unique())}')
print(f'DRSP genes we do not model: {sorted(set(tr.Gene) - set(GENE2TGT))}')

tr['y_drsp'] = (tr.Drug != 'NA_deleterious').astype(int)
print(f'DRSP training labels: {int(tr.y_drsp.sum())} resistant / {int((1-tr.y_drsp).sum())} '
      f'NA_deleterious  ({tr.y_drsp.mean():.1%} positive)')

def detect_offset(savs, wt_at):
    """Find the integer offset o with DRSP position = our position + o, by wild-type agreement.

    DRSP does not state which isoform it numbers on, and for ABL1 it is demonstrably not ours:
    its M263V is our M244V (the +19 isoform-b offset). Rather than assume, the offset that
    maximises wild-type-residue agreement is chosen from a window, and both the agreement at
    that offset and the agreement at 0 are reported so the choice is visible.
    """
    best = (None, -1, 0)
    for o in range(-60, 61):
        hit = tot = 0
        for s in savs:
            p = int(s[1:-1]) - o
            if p in wt_at:
                tot += 1
                hit += (wt_at[p] == s[0])
        if tot >= 5 and hit > best[1]:
            best = (o, hit, tot)
    return best


rows, missing, wtbad, offsets = [], [], [], {}
for g, tgt in GENE2TGT.items():
    sub = tr[tr.Gene == g]
    ours = pd.read_csv(f'{TAB}/{tgt}.csv')
    ours['mut'] = ours.aa.str.replace(r'^p\.', '', regex=True)
    byname = ours.set_index('mut')
    # numbering check: our table's wild-type residue at that position must equal DRSP's
    wt_at = {}
    for m in ours['mut']:
        wt_at.setdefault(int(m[1:-1]), m[0])
    savs = [str(x).strip() for x in sub.SAV]
    o, hit, tot = detect_offset(savs, wt_at)
    hit0 = sum(wt_at.get(int(s[1:-1]), '') == s[0] for s in savs)
    tot0 = sum(int(s[1:-1]) in wt_at for s in savs)
    offsets[g] = dict(offset=o, wt_agree=f'{hit}/{tot}', wt_agree_at_0=f'{hit0}/{tot0}')
    if tot and hit / tot < 0.9:
        print(f'  !! {g}: best offset {o} still only agrees {hit}/{tot} -- not trusted')
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
        sav = str(r.SAV).strip()
        pos = int(sav[1:-1]) - o                       # translate into our numbering
        ours_sav = f'{sav[0]}{pos}{sav[-1]}'
        if pos in wt_at and wt_at[pos] != sav[0]:
            wtbad.append((g, sav, ours_sav, wt_at[pos])); continue
        if ours_sav not in byname.index:
            missing.append((g, sav, ours_sav, r.Drug)); continue
        row = byname.loc[ours_sav]
        rows.append(dict(gene=g, target=tgt, drsp_sav=sav, sav=ours_sav, pos=pos,
                         drsp_drug=r.Drug, y_drsp=int(r.y_drsp), y_ours=int(row.y),
                         our_drugs=('' if pd.isna(row.drugs) else row.drugs),
                         loto_pa=lopa.get(ours_sav, np.nan),
                         cv_rank=(float(cv.loc[ours_sav]) if cv is not None
                                  and ours_sav in cv.index else np.nan)))
D = pd.DataFrame(rows)
print('\nnumbering offsets detected (DRSP position = our position + offset):')
print(pd.DataFrame(offsets).T.to_string())
print(f'\nmatched {len(D)} of {len(tr[tr.Gene.isin(GENE2TGT)])} DRSP variants in the 9 shared genes')
print(f'  wild-type residue mismatches (numbering): {len(wtbad)}'
      + (f'  e.g. {wtbad[:5]}' if wtbad else ''))
print(f'  absent from our COSMIC-derived variant table: {len(missing)}'
      + (f'  e.g. {missing[:5]}' if missing else ''))
print(f'  LOTO score present for {int(D.loto_pa.notna().sum())} / {len(D)}')
print(f'  within-protein CV score present for {int(D.cv_rank.notna().sum())} / {len(D)}')
D.to_csv(f'{_SAVDR_ROOT}/multipass/drsp_overlap.csv', index=False)

print('\nper gene (DRSP labels on the matched set):')
print(D.groupby('gene').agg(n=('sav', 'size'), drsp_pos=('y_drsp', 'sum'),
                            our_pos=('y_ours', 'sum')).to_string())

# ---- label agreement between DRSP and our COSMIC labels ------------------------------------
tab = pd.crosstab(D.y_drsp, D.y_ours)
print('\nDRSP label (rows) vs our COSMIC label (cols):')
print(tab.to_string())
agree = (D.y_drsp == D.y_ours).mean()
print(f'  agreement {agree:.1%}; Cohen kappa '
      f'{(agree - ((D.y_drsp.mean()*D.y_ours.mean())+((1-D.y_drsp.mean())*(1-D.y_ours.mean())))) / (1 - ((D.y_drsp.mean()*D.y_ours.mean())+((1-D.y_drsp.mean())*(1-D.y_ours.mean())))):.3f}')


def panel(d, score, name, y='y_drsp'):
    d = d[d[score].notna()].copy()
    yv = d[y].values
    s = np.nan_to_num(d[score].values, nan=0.0)
    if len(set(yv)) < 2:
        print(f'  {name}: only one class, skipped'); return None
    lo, hi = ci(yv, s)
    auc = roc_auc_score(yv, s)
    out = dict(panel=name, n=len(d), n_pos=int(yv.sum()), n_neg=int((1 - yv).sum()),
               auc=round(auc, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
               auprc=round(average_precision_score(yv, s), 3))
    # threshold-based metrics at the prevalence-matched cut (top-k = number of positives),
    # which needs no calibration, plus the best achievable cut as an upper bound
    k = int(yv.sum())
    thr = np.sort(s)[::-1][k - 1]
    pred = (s >= thr).astype(int)
    out |= dict(acc_prev=round(accuracy_score(yv, pred), 3),
                mcc_prev=round(matthews_corrcoef(yv, pred), 3),
                f1_prev=round(f1_score(yv, pred), 3),
                prec_prev=round(precision_score(yv, pred, zero_division=0), 3),
                rec_prev=round(recall_score(yv, pred), 3))
    best = max(((matthews_corrcoef(yv, (s >= t).astype(int)), t) for t in np.unique(s)),
               key=lambda x: x[0])
    pb = (s >= best[1]).astype(int)
    out |= dict(mcc_best=round(best[0], 3), acc_at_best=round(accuracy_score(yv, pb), 3),
                f1_at_best=round(f1_score(yv, pb), 3))
    out['acc_all_negative'] = round(1 - yv.mean(), 3)
    return out


print('\n' + '=' * 104)
print('MultiPASS on DRSP\'s own variants, scored against DRSP\'s own labels')
print('=' * 104)
res = []
for score, name in [('loto_pa', 'A. LOTO, 9 shared genes'),
                    ('cv_rank', 'B. within-protein CV, 5 genes')]:
    r = panel(D, score, name)
    if r: res.append(r)
# panel A restricted to the same 5 genes as B, so A and B are comparable
r = panel(D[D.target.isin(CVDIR)], 'loto_pa', 'A5. LOTO, the same 5 genes as B')
if r: res.append(r)
R = pd.DataFrame(res)
pd.set_option('display.width', 250)
print(R.to_string(index=False))
R.to_csv(f'{_SAVDR_ROOT}/multipass/drsp_compare.csv', index=False)

# ---- DRSP's own published numbers, and an arithmetic check on them -------------------------
# Table 1 of Lin et al. 2022, DRSPce+DRSPci (crystal, combined), 5-fold CV on the 610-SAV
# training set: accuracy 0.8557, sensitivity 0.6541, specificity 0.9119, MCC 0.5724,
# precision 0.6744, F1 0.6641. Abstract rounds these to 86% / 0.57 / 0.66.
NP, NN = 133, 477
SENS, SPEC = 0.6541, 0.9119
TP = round(SENS * NP); FN = NP - TP
TN = round(SPEC * NN); FP = NN - TN
acc = (TP + TN) / (NP + NN)
prec = TP / (TP + FP)
f1 = 2 * prec * SENS / (prec + SENS)
mcc = ((TP * TN - FP * FN)
       / np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
print('\nDRSP\'s own published figures (Lin et al. 2022 Table 1, crystal combined, 5-fold CV')
print('on their 610-SAV training set) and an arithmetic re-derivation from sens/spec:')
print(f'  reported   accuracy 0.8557  MCC 0.5724  F1 0.6641  precision 0.6744')
print(f'  re-derived accuracy {acc:.4f}  MCC {mcc:.4f}  F1 {f1:.4f}  precision {prec:.4f}'
      f'   (TP {TP} FP {FP} FN {FN} TN {TN})')
print(f'  all-negative baseline on their set: {NN/(NP+NN):.3f}; '
      f'on our matched subset: {1 - D.y_drsp.mean():.3f}')
print('  Their figures are cross-validated on the set they fitted, and their feature set was')
print('  chosen by a genetic algorithm on that same data, best of 8 runs by optimised MCC --')
print('  the same selection optimism we measured for our own grid. Ours are out-of-fold.')

# ---- matched operating point: our sensitivity at DRSP's specificity ------------------------
print('\n' + '=' * 104)
print(f'MATCHED OPERATING POINT -- our sensitivity at DRSP\'s specificity ({SPEC:.4f})')
print('=' * 104)
print('Accuracy is prevalence-dependent and the two sets differ in prevalence, so the')
print('comparable quantity is sensitivity at a fixed specificity.')


def sens_at_spec(y, s, spec):
    neg = np.sort(s[y == 0])[::-1]
    k = int(np.floor((1 - spec) * len(neg)))
    thr = neg[k] + 1e-12 if k < len(neg) else -np.inf
    pred = (s > thr).astype(int)
    return (pred[y == 1].mean(), float((1 - pred[y == 0]).mean()), thr)


rows2 = []
for score, name, sub in [('cv_rank', 'MultiPASS, within-protein CV (5 genes)', D),
                         ('loto_pa', 'MultiPASS, LOTO (9 genes)', D),
                         ('loto_pa', 'MultiPASS, LOTO (same 5 genes)', D[D.target.isin(CVDIR)])]:
    d = sub[sub[score].notna()]
    y = d.y_drsp.values
    s = np.nan_to_num(d[score].values, nan=0.0)
    se, sp, _ = sens_at_spec(y, s, SPEC)
    bs = []
    rng = np.random.default_rng(0)
    ip, ineg = np.where(y == 1)[0], np.where(y == 0)[0]
    for _ in range(2000):
        idx = np.concatenate([rng.choice(ip, len(ip), True), rng.choice(ineg, len(ineg), True)])
        bs.append(sens_at_spec(y[idx], s[idx], SPEC)[0])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows2.append(dict(score=name, n=len(d), n_pos=int(y.sum()),
                      spec_achieved=round(sp, 4), sens=round(se, 3),
                      sens_ci=f'{lo:.3f}-{hi:.3f}'))
    print(f'  {name:42s} n={len(d):4d}  specificity {sp:.4f}  '
          f'sensitivity {se:.3f} [{lo:.3f}, {hi:.3f}]')
print(f'  {"DRSP crystal combined (published)":42s} n={NP+NN:4d}  specificity {SPEC:.4f}  '
      f'sensitivity {SENS:.3f}  (their own CV, their own fit)')
pd.DataFrame(rows2).to_csv(f'{_SAVDR_ROOT}/multipass/drsp_operating_point.csv',
                           index=False)

print('\nper gene, LOTO panel:')
pg = []
for g, sub in D.groupby('gene'):
    r = panel(sub, 'loto_pa', g)
    if r: pg.append(r)
PG = pd.DataFrame(pg)
print(PG.to_string(index=False))
PG.to_csv(f'{_SAVDR_ROOT}/multipass/drsp_pergene.csv', index=False)
