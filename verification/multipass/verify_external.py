"""Independent re-verification of the SUSPECT-ABL and Dr. Kinase comparisons.

Every number is recomputed by a different route than the script that produced it, or an
invariant is tested that would fail if that script were wrong. AUCs are recomputed from the
Mann-Whitney U statistic rather than sklearn's roc_auc_score; the DeLong p-value is checked
against an independent paired bootstrap of the AUC difference.
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

import glob, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, matthews_corrcoef
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])

ok, bad = [], []
def chk(name, cond, detail=''):
    (ok if cond else bad).append(name)
    print(('  PASS  ' if cond else '  FAIL  ') + name + ('   ' + detail if detail else ''))

FOUR = ['imatinib', 'dasatinib', 'nilotinib', 'bosutinib']
S = pd.read_csv(f'{_SAVDR_ROOT}/inputs/suspect/suspect_abl.csv')
SET = pd.read_csv(f'{_SAVDR_ROOT}/inputs/suspect/abl1_suspect_set.csv')
TAB = pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/ABL1_ENST00000318560.csv')

print('=' * 100); print('1. THE VARIANT SET'); print('=' * 100)
ref = TAB[(TAB.pos >= 242) & (TAB.pos <= 493)].reset_index(drop=True)
chk('set is exactly our ABL1 variants in 242-493',
    sorted(SET.aa) == sorted(ref.aa) and len(SET) == 280, f'{len(SET)} variants')
chk('positives in the window match the source table',
    int(SET.y.sum()) == int(ref.y.sum()) == 97, f'{int(SET.y.sum())}')
chk('ABL1 total positives unchanged from the main analysis', int(TAB.y.sum()) == 101)
chk('SUSPECT returned every variant exactly once',
    len(S) == 280 and S.mut.nunique() == 280)
chk('SUSPECT mutations == our variant names minus the p. prefix',
    sorted(S.mut) == sorted(SET.aa.str.replace(r'^p\.', '', regex=True)))

print('=' * 100); print('2. LABEL CONSTRUCTION, recomputed from scratch'); print('=' * 100)
d = SET.copy()
d['mut'] = d.aa.str.replace(r'^p\.', '', regex=True)
d = d.merge(S, on='mut', validate='one_to_one')
def drugset(s):
    return {x.strip().lower() for x in str(s).split(';') if x and str(x) != 'nan'}
y4_ind = np.array([int(bool(drugset(s) & set(FOUR))) for s in d.drugs.fillna('')])
excl_ind = (d.y.values == 1) & (y4_ind == 0)
chk('10 positives are TKI-NS-only and get excluded', int(excl_ind.sum()) == 10)
for v in d.loc[excl_ind, 'drugs']:
    if drugset(v) != {'tyrosine kinase inhibitor - ns'}:
        chk('an excluded variant is TKI-NS-only', False, str(v)); break
else:
    chk('every excluded variant is TKI-NS-only and nothing else', True)
k = ~excl_ind
y = y4_ind[k]
chk('final set 270 = 87 positives + 183 negatives',
    len(y) == 270 and int(y.sum()) == 87 and int((1 - y).sum()) == 183)
chk('no negative carries any resistance annotation',
    all(drugset(s) == set() for s in d.loc[k, 'drugs'].fillna('')[y == 0]))

print('=' * 100); print('3. ddG PARSING spot-checks against the raw server output'); print('=' * 100)
spot = {('T315I', 'imatinib_ddg'): 2.32, ('T315I', 'dasatinib_ddg'): 5.08,
        ('T315I', 'nilotinib_ddg'): 3.75, ('T315I', 'bosutinib_ddg'): 2.45,
        ('M244V', 'imatinib_ddg'): 0.21, ('E459K', 'imatinib_ddg'): 0.52,
        ('I242M', 'imatinib_ddg'): 0.12, ('F493L', 'bosutinib_ddg'): 0.24,
        ('H246D', 'imatinib_ddg'): -0.27}
allgood = True
for (mu, col), want in spot.items():
    got = float(S.loc[S.mut == mu, col].iloc[0])
    if abs(got - want) > 1e-9: chk(f'{mu}/{col}', False, f'{got} vs {want}'); allgood = False
chk(f'{len(spot)} ddG cells match the raw server text', allgood)
chk('T315I is called Resistant for all four inhibitors',
    all(S.loc[S.mut == 'T315I', f'{x}_cls'].iloc[0] == 'R' for x in FOUR))
chk('no ddG cell is missing', S[[f'{x}_ddg' for x in FOUR]].notna().all().all())
chk('every class cell is R or S',
    set(np.ravel(S[[f'{x}_cls' for x in FOUR]].values)) <= {'R', 'S'})

print('=' * 100); print('4. AUCs RECOMPUTED from Mann-Whitney U (not sklearn)'); print('=' * 100)
ddg = d.loc[k, [f'{x}_ddg' for x in FOUR]].astype(float)
sc = {'SUSPECT-ABL (max ddG)': ddg.max(axis=1).values,
      'SUSPECT-ABL (mean ddG)': ddg.mean(axis=1).values}
per = []
for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_ABL1_ENST00000318560_17_9_s*/scores.csv')):
    t = pd.read_csv(f)
    per.append(pd.Series(stats.rankdata(np.nan_to_num(t.multipass_pa.values, nan=0.0)) / len(t),
                         index=t.aa.values))
chk('MultiPASS seeds found', len(per) == 5, f'{len(per)} seeds')
mp = pd.concat(per, axis=1).mean(axis=1)
sc['MultiPASS'] = mp.loc[d.loc[k, 'aa']].values

def auc_mw(y_, s_):
    u = stats.mannwhitneyu(s_[y_ == 1], s_[y_ == 0], alternative='two-sided', method='asymptotic')
    return u.statistic / ((y_ == 1).sum() * (y_ == 0).sum())

T = pd.read_csv(f'{_SAVDR_ROOT}/multipass/suspect_abl_compare.csv').set_index('model')
for name, v in sc.items():
    a_mw, a_sk = auc_mw(y, v), roc_auc_score(y, v)
    want = float(T.loc[name, 'auc'])
    chk(f'{name}: U-statistic AUC == sklearn == reported',
        abs(a_mw - a_sk) < 1e-9 and abs(round(a_sk, 3) - want) < 6e-4,
        f'U {a_mw:.4f} / sklearn {a_sk:.4f} / reported {want}')

print('=' * 100); print('5. DeLong p-value vs an INDEPENDENT paired bootstrap'); print('=' * 100)
for other in ['SUSPECT-ABL (max ddG)', 'SUSPECT-ABL (mean ddG)']:
    _, p_dl = delong(y, sc['MultiPASS'], sc[other])
    rng = np.random.default_rng(11)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    diffs = []
    for _ in range(3000):
        i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        yy = y[i]
        if yy.sum() in (0, len(yy)): continue
        diffs.append(roc_auc_score(yy, sc['MultiPASS'][i]) - roc_auc_score(yy, sc[other][i]))
    diffs = np.array(diffs)
    p_bs = 2 * min((diffs > 0).mean(), (diffs < 0).mean())
    chk(f'DeLong vs bootstrap, MultiPASS vs {other}',
        abs(p_dl - p_bs) < 0.15,
        f'DeLong {p_dl:.3f}, bootstrap {p_bs:.3f}, mean diff {diffs.mean():+.4f}')

print('=' * 100); print('6. SUSPECT-ABL CLASSIFIER confusion matrix, recomputed'); print('=' * 100)
cls_any = np.array([int(any(d.loc[k].iloc[i][f'{x}_cls'] == 'R' for x in FOUR)) for i in range(len(y))])
tp = int(((cls_any == 1) & (y == 1)).sum()); fp = int(((cls_any == 1) & (y == 0)).sum())
fn = int(((cls_any == 0) & (y == 1)).sum()); tn = int(((cls_any == 0) & (y == 0)).sum())
chk('TP/FP/FN/TN = 14/17/73/166', (tp, fp, fn, tn) == (14, 17, 73, 166), f'{tp}/{fp}/{fn}/{tn}')
chk('counts sum to 270', tp + fp + fn + tn == 270)
chk('sensitivity 0.161', abs(tp / (tp + fn) - 0.161) < 6e-4, f'{tp/(tp+fn):.4f}')
chk('specificity 0.907', abs(tn / (tn + fp) - 0.907) < 6e-4, f'{tn/(tn+fp):.4f}')
chk('MCC 0.100', abs(matthews_corrcoef(y, cls_any) - 0.100) < 6e-4,
    f'{matthews_corrcoef(y, cls_any):.4f}')

print('=' * 100); print('7. DRUG-MATCHED RECALL, recomputed'); print('=' * 100)
hits = tot = 0; pd_ = {x: [0, 0] for x in FOUR}
for _, r in d.iterrows():
    for dr in drugset(r.drugs) & set(FOUR):
        tot += 1; pd_[dr][1] += 1
        if r[f'{dr}_cls'] == 'R': hits += 1; pd_[dr][0] += 1
chk('overall drug-matched recall 16/106', (hits, tot) == (16, 106), f'{hits}/{tot}')
for dr, (h, t) in pd_.items():
    want = {'imatinib': (5, 82), 'dasatinib': (5, 10), 'nilotinib': (4, 9), 'bosutinib': (2, 5)}[dr]
    chk(f'  {dr} {want[0]}/{want[1]}', (h, t) == want, f'{h}/{t}')
anyR = int(sum(1 for _, r in d.iterrows() if any(r[f'{x}_cls'] == 'R' for x in FOUR)))
chk('SUSPECT calls R for 32 of 280 (11.4%)', anyR == 32, f'{anyR}')
chk('our prevalence in the window 0.346', abs(d.y.mean() - 0.346) < 6e-4, f'{d.y.mean():.4f}')

print('=' * 100); print('8. IMATINIB per-drug numbers'); print('=' * 100)
# The per-drug analysis runs on the 270-variant set, i.e. AFTER the 10 TKI-NS-only positives
# are removed. Counting the exclusions on the 280 base instead gives 5 + 10 = 15, which is what
# an earlier version of this check did; the base, not the analysis, was wrong.
d270 = d[~excl_ind].reset_index(drop=True)
im = np.array([int('imatinib' in drugset(s)) for s in d270.drugs.fillna('')])
other_only = (d270.y.values == 1) & (im == 0)
kk = ~other_only
yi = im[kk]
chk('82 imatinib positives vs 183 negatives, 5 excluded (on the 270-variant base)',
    int(yi.sum()) == 82 and int((1 - yi).sum()) == 183 and int(other_only.sum()) == 5,
    f'{int(yi.sum())}/{int((1-yi).sum())}, {int(other_only.sum())} excluded')
chk('on the 280 base the same exclusion is 15 (= 5 + the 10 TKI-NS-only)',
    int((((d.y.values == 1) & (np.array([int("imatinib" in drugset(s)) for s in d.drugs.fillna("")]) == 0))).sum()) == 15)
a_sus = roc_auc_score(yi, d270.imatinib_ddg.values[kk])
chk('SUSPECT ddG(imatinib) AUC 0.560', abs(round(a_sus, 3) - 0.560) < 6e-4, f'{a_sus:.4f}')
pdr = sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/drug_ABL1_ENST00000318560_17_9_s*/scores_drug.csv'))
pers = []
for f in pdr:
    t = pd.read_csv(f)
    c = [x for x in t.columns if x.lower() == 'pa::imatinib']
    if c: pers.append(pd.Series(stats.rankdata(np.nan_to_num(t[c[0]].values, nan=0.0)) / len(t),
                                index=t.aa.values))
chk('per-drug imatinib seeds found', len(pers) == 5, f'{len(pers)}')
if pers:
    mpi = pd.concat(pers, axis=1).mean(axis=1).loc[d270.loc[kk, "aa"]].values
    a_mp = roc_auc_score(yi, mpi)
    chk('MultiPASS Pa(Imatinib) AUC 0.571', abs(round(a_mp, 3) - 0.571) < 6e-4, f'{a_mp:.4f}')
    _, pv = delong(yi, mpi, d270.imatinib_ddg.values[kk])
    chk('imatinib DeLong p 0.824', abs(pv - 0.8242) < 5e-3, f'{pv:.4f}')

print('=' * 100); print('9. Dr. KINASE'); print('=' * 100)
DK = pd.read_csv(f'{_SAVDR_ROOT}/multipass/drkinase_abl1.csv')
chk('280 rows, 9 scored', len(DK) == 280 and int(DK.drkinase_gk.notna().sum()) == 9)
chk('6 of the 9 scored are our positives', int(DK[DK.drkinase_gk.notna()].y.sum()) == 6)
chk('Dr. Kinase scored 6 of our 97 positives in the window',
    int(DK[(DK.y == 1) & DK.drkinase_gk.notna()].shape[0]) == 6)
chk('T315I scores 0.5495, below the 0.90 reference line',
    abs(float(DK.loc[DK.mut == 'T315I', 'drkinase_gk'].iloc[0]) - 0.5495) < 1e-9)
chk('only M351K is at or above 0.90',
    list(DK[DK.drkinase_gk >= 0.90].mut) == ['M351K'], str(list(DK[DK.drkinase_gk >= 0.90].mut)))
# the substitution-awareness claim: distinct scores at a shared position
for p_, n_ in [(315, 3), (351, 3)]:
    vals = sorted(DK[(DK.pos == p_) & DK.drkinase_gk.notna()].drkinase_gk.round(4).tolist())
    chk(f'position {p_}: {n_} substitutions with distinct scores',
        len(vals) == n_ and len(set(vals)) == n_, str(vals))
chk('unscored share is 271/280 = 96.8%',
    int(DK.drkinase_gk.isna().sum()) == 271 and abs(271/280 - 0.968) < 6e-4)

print('=' * 100); print('10. NARRATIVE NUMBERS quoted from the other reports'); print('=' * 100)
ms = pd.read_csv(f'{_SAVDR_ROOT}/multipass/margin_sweep.csv')
m = ms[(ms.target == 'ABL1_ENST00000318560') & (ms.model == 'multipass_full_trained')].set_index('margin')
for mg, want in [('5', 0.589), ('10', 0.599), ('full', 0.860)]:
    chk(f'ABL1 margin {mg} = {want}', abs(float(m.loc[mg, 'auc']) - want) < 6e-4,
        f'{float(m.loc[mg,"auc"]):.3f}')
chk('positive span inside 242-493', SET[SET.y == 1].pos.min() == 244 and SET[SET.y == 1].pos.max() == 486,
    f'{SET[SET.y==1].pos.min()}-{SET[SET.y==1].pos.max()}')

print('=' * 100)
print(f'{len(ok)} checks passed, {len(bad)} failed')
if bad:
    print('FAILURES:'); [print('  ' + b) for b in bad]
print("""
NOT verified here:
  * whether SUSPECT-ABL's or Dr. Kinase's training data actually overlaps COSMIC v100 -- both
    are inferred from their papers, not established, and neither can be retrained.
  * that the four inhibitor names map to the same co-crystals our drug labels implicitly assume.
  * that every one of the 280 variants was genuinely submitted to Dr. Kinase: 280 results came
    back and the typed list matches the source exactly, but with 271 unscored a silently
    dropped submission would be indistinguishable from an unscored one.
  * the three Dr. Kinase channels other than gatekeeper, for which the web endpoint returns only
    top-loci lists and no per-residue track.
  * MultiPASS's internals (closed Win32 binary).
""")
