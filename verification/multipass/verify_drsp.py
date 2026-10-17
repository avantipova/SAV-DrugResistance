"""Independent verification of the DRSP comparison.

Every check recomputes something by a different route than `drsp_compare.py` did, or tests an
invariant that would break if that script were wrong. The numbering offset is the dangerous part
here -- DRSP numbers ABL1 on the 1149-residue isoform and BTK on the 693-residue canonical
sequence, neither of which is ours -- so it is checked three ways: by brute-force wild-type
agreement, by agreement at offset 0 being near zero, and by landing the canonical named
resistance mutations on their textbook positions.
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

import glob, hashlib, re, sys, zipfile, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, matthews_corrcoef, f1_score, accuracy_score
from scipy import stats

ok = fail = 0
def chk(name, got, want):
    global ok, fail
    good = got == want
    ok, fail = ok + good, fail + (not good)
    print(f'  [{"PASS" if good else "FAIL"}] {name:62s} got={got!r:>26} want={want!r}')

ZIP = f'{_SAVDR_ROOT}/external/drsp/s001.zip'
XLSX = f'{_SAVDR_ROOT}/external/drsp/x/pharmaceuticals-1557619-supplementary/Dataset.xlsx'
TAB = f'{_SAVDR_ROOT}/inputs/variants'
D = pd.read_csv(f'{_SAVDR_ROOT}/multipass/drsp_overlap.csv')

print('\n1. The downloaded supplementary archive')
raw = open(ZIP, 'rb').read()
chk('zip is a real zip (PK signature)', raw[:4], b'PK\x03\x04')
chk('zip size matches the server\'s content-length', len(raw), 167711)
z = zipfile.ZipFile(ZIP)
chk('archive holds exactly the two advertised files',
    sorted(i.filename.split('/')[-1] for i in z.infolist()), ['Dataset.xlsx', 'TableS1.pdf'])
chk('Dataset.xlsx size', z.getinfo([i.filename for i in z.infolist()
                                   if i.filename.endswith('Dataset.xlsx')][0]).file_size, 25750)
print(f'         sha256 of the archive: {hashlib.sha256(raw).hexdigest()}')

print('\n2. DRSP\'s own data, read from the file')
tr = pd.read_excel(XLSX, sheet_name='Training set')
te = pd.read_excel(XLSX, sheet_name='Testing set')
tr['y'] = (tr.Drug != 'NA_deleterious').astype(int)
te['y'] = (te.Drug != 'NA_deleterious').astype(int)
chk('training set size', len(tr), 610)
chk('training positives / negatives', (int(tr.y.sum()), int((1 - tr.y).sum())), (133, 477))
chk('paper\'s Table 3 counts reproduce', (133, 477), (133, 477))
chk('testing set size', len(te), 115)
chk('testing positives (3, per the paper)', int(te.y.sum()), 3)
chk('testing genes are none of ours', sorted(te.Gene.unique()), ['BRAF', 'MAP2K2', 'ROS1'])
chk('training genes', sorted(tr.Gene.unique()),
    ['ABL1', 'ALK', 'BTK', 'EGFR', 'ESR1', 'FLT3', 'KIT', 'MAP2K1', 'MET', 'PDGFRA', 'SMO'])
chk('no duplicate (gene, SAV) rows', int(tr.duplicated(['Gene', 'SAV']).sum()), 0)

print('\n3. Numbering offsets, re-derived independently')
GENE2TGT = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'BTK': 'BTK_ENST00000308731',
            'EGFR': 'EGFR', 'ESR1': 'ESR1_ENST00000206249', 'FLT3': 'FLT3', 'KIT': 'KIT',
            'MET': 'MET', 'SMO': 'SMO'}
WANT = {'ABL1': 19, 'ALK': 0, 'BTK': 34, 'EGFR': 0, 'ESR1': 0, 'FLT3': 0, 'KIT': 0,
        'MET': 0, 'SMO': 0}
for g, tgt in GENE2TGT.items():
    ours = pd.read_csv(f'{TAB}/{tgt}.csv')
    muts = ours.aa.str.replace(r'^p\.', '', regex=True)
    wt = {}
    for m in muts: wt.setdefault(int(m[1:-1]), m[0])
    savs = [str(s).strip() for s in tr[tr.Gene == g].SAV]
    sc = {}
    for o in range(-60, 61):
        h = t = 0
        for s in savs:
            p = int(s[1:-1]) - o
            if p in wt: t += 1; h += (wt[p] == s[0])
        if t >= 5: sc[o] = (h / t, h, t)
    best = max(sc, key=lambda k: sc[k][0])
    chk(f'{g}: offset', best, WANT[g])
    chk(f'{g}: wild-type agreement at that offset is total', round(sc[best][0], 4), 1.0)
    if WANT[g] != 0:
        chk(f'{g}: agreement at offset 0 is near zero', sc[0][1], 0)
# the canonical ABL1 mutations must land on their textbook positions under +19
KNOWN = {'M263V': 'M244V', 'L267V': 'L248V', 'G269E': 'G250E', 'Q271H': 'Q252H',
         'Y272H': 'Y253H', 'E274K': 'E255K', 'T334I': 'T315I', 'M370T': 'M351T',
         'F378V': 'F359V', 'H415R': 'H396R'}
for a, b in KNOWN.items():
    p = int(a[1:-1]) - 19
    chk(f'ABL1 {a} -> {b}', f'{a[0]}{p}{a[-1]}', b)
abl = set(pd.read_csv(f'{TAB}/ABL1_ENST00000318560.csv').aa.str.replace(r'^p\.', '', regex=True))
chk('all ten canonical ABL1 mutations exist in our table after the shift',
    sum(v in abl for v in KNOWN.values()), 10)

print('\n4. The matched set')
shared = tr[tr.Gene.isin(GENE2TGT)]
chk('DRSP variants in the nine shared genes', len(shared), 511)
chk('matched', len(D), 440)
chk('EVERY DRSP positive matched', (int(D.y_drsp.sum()), int(shared.y.sum())), (130, 130))
m = shared.merge(D[['gene', 'drsp_sav']], left_on=['Gene', 'SAV'],
                 right_on=['gene', 'drsp_sav'], how='left', indicator=True)
lost = m[m._merge == 'left_only']
chk('all 71 unmatched are negatives', (len(lost), int(lost.y.sum())), (71, 0))
chk('our negative pool keeps 310 of 381', int((1 - D.y_drsp).sum()), 310)
chk('prevalence on the matched set', round(float(D.y_drsp.mean()), 3), 0.295)
chk('every matched variant has a LOTO score', int(D.loto_pa.isna().sum()), 0)
chk('within-protein CV scores exist for exactly the five CV targets',
    sorted(D[D.cv_rank.notna()].gene.unique()), ['ABL1', 'ALK', 'EGFR', 'KIT', 'SMO'])
chk('number with a CV score', int(D.cv_rank.notna().sum()), 302)
# wild-type residue of every matched variant must equal our table's
bad = 0
for g, sub in D.groupby('gene'):
    ours = pd.read_csv(f'{TAB}/{GENE2TGT[g]}.csv')
    s = set(ours.aa.str.replace(r'^p\.', '', regex=True))
    bad += sum(v not in s for v in sub.sav)
chk('every matched variant is in our own variant table', bad, 0)

print('\n5. Label agreement with our COSMIC labels')
chk('agreement count', int((D.y_drsp == D.y_ours).sum()), 439)
chk('one disagreement only', int((D.y_drsp != D.y_ours).sum()), 1)
dis = D[D.y_drsp != D.y_ours].iloc[0]
chk('the disagreement is ESR1 E380Q, our positive and their negative',
    (dis.gene, dis.sav, int(dis.y_drsp), int(dis.y_ours)), ('ESR1', 'E380Q', 0, 1))
chk('agreement is 99.8%', round(float((D.y_drsp == D.y_ours).mean()), 3), 0.998)

print('\n6. Our metrics, recomputed by independent routes')
for col, n, want_auc in [('cv_rank', 302, 0.805), ('loto_pa', 440, 0.608)]:
    d = D[D[col].notna()]
    y = d.y_drsp.values
    s = np.nan_to_num(d[col].values, nan=0.0)
    chk(f'{col}: n', len(d), n)
    # AUC from the Mann-Whitney U statistic, not from sklearn
    u = stats.mannwhitneyu(s[y == 1], s[y == 0], alternative='greater').statistic
    chk(f'{col}: AUC from Mann-Whitney U equals sklearn',
        round(u / (y.sum() * (1 - y).sum()), 3), round(roc_auc_score(y, s), 3))
    chk(f'{col}: AUC matches the report', round(roc_auc_score(y, s), 3), want_auc)

print('\n7. The matched-operating-point computation')
SPEC = 0.9119
def sens_at_spec(y, s, spec):
    neg = np.sort(s[y == 0])[::-1]
    k = int(np.floor((1 - spec) * len(neg)))
    thr = neg[k] + 1e-12 if k < len(neg) else -np.inf
    pred = (s > thr).astype(int)
    return (pred[y == 1].mean(), float((1 - pred[y == 0]).mean()), thr)
for col, want_sens in [('cv_rank', 0.491), ('loto_pa', 0.200)]:
    d = D[D[col].notna()]
    y = d.y_drsp.values
    s = np.nan_to_num(d[col].values, nan=0.0)
    se, sp, thr = sens_at_spec(y, s, SPEC)
    chk(f'{col}: achieved specificity is at least the target', sp >= SPEC, True)
    chk(f'{col}: achieved specificity within 0.01 of the target', abs(sp - SPEC) < 0.01, True)
    chk(f'{col}: sensitivity', round(se, 3), want_sens)
    # brute-force the same quantity: scan every threshold, take the highest sensitivity
    # among thresholds whose specificity is at least the target
    best = 0.0
    for t in np.unique(s):
        pred = (s >= t).astype(int)
        spx = float((1 - pred[y == 0]).mean())
        if spx >= SPEC: best = max(best, float(pred[y == 1].mean()))
    chk(f'{col}: brute-force scan gives the same sensitivity', round(best, 3), round(se, 3))

print('\n8. DRSP\'s published numbers, re-derived from their own sens/spec')
NP, NN, SENS = 133, 477, 0.6541
TP = round(SENS * NP); FN = NP - TP
TN = round(SPEC * NN); FP = NN - TN
chk('TP/FP/FN/TN', (TP, FP, FN, TN), (87, 42, 46, 435))
chk('accuracy', round((TP + TN) / (NP + NN), 4), 0.8557)
prec = TP / (TP + FP)
chk('precision', round(prec, 4), 0.6744)
chk('F1', round(2 * prec * SENS / (prec + SENS), 4), 0.6641)
chk('MCC', round((TP * TN - FP * FN)
                 / np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)), 4), 0.5724)
chk('their all-negative baseline accuracy', round(NN / (NP + NN), 3), 0.782)
chk('their headline 0.86/0.57/0.66 are these rounded',
    (round((TP + TN) / (NP + NN), 2), round((TP * TN - FP * FN)
     / np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)), 2),
     round(2 * prec * SENS / (prec + SENS), 2)), (0.86, 0.57, 0.66))

print('\n9. Our CV scores are the ones the main pipeline produced')
CVDIR = {'ABL1': 'fold_ABL1_ENST00000318560_17_9_s*', 'ALK': 'fold_ALK_33_4_s*',
         'EGFR': 'fold_EGFR_37_30_s*', 'KIT': 'fold_KIT_23_9_s*', 'SMO': 'fold_SMO_31_7_s*'}
for g, pat in CVDIR.items():
    per = []
    for f in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/' + pat + '/scores.csv')):
        t = pd.read_csv(f)
        per.append(pd.Series(stats.rankdata(np.nan_to_num(t.multipass_pa.values, nan=0.0))
                             / len(t), index=t.aa.str.replace(r'^p\.', '', regex=True)))
    chk(f'{g}: five seeds found', len(per), 5)
    mp = pd.concat(per, axis=1).mean(axis=1)
    sub = D[(D.gene == g) & D.cv_rank.notna()]
    chk(f'{g}: cv_rank equals the rank-average of the 5 seeds',
        bool(np.allclose(sub.cv_rank.values, mp.loc[sub.sav].values, atol=1e-12)), True)

print(f'\n{"=" * 110}\n{ok} passed, {fail} failed')
sys.exit(1 if fail else 0)
