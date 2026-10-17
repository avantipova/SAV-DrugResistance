"""Verify every number in civic/REPORT_CIVIC.md against its source artefacts.

Sources of truth:
  civic_external_eval.csv   -- the evaluation tables (A/B x 6 targets x 4 models + POOLED)
  civic_overlap.log         -- the 90/70/14/6 overlap table and the discordance list
  therapy_match.log         -- the therapy-matched analysis
  civic_eval.py             -- whether an `interpretable` column is actually written

This was written after the report was found to be 19 hours older than the CSV it quotes.
Each check prints PASS/FAIL with the report value and the source value side by side.
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

import re, sys, pandas as pd

R = open(f'{_SAVDR_ROOT}/reports/REPORT_CIVIC.md').read()
E = pd.read_csv(f'{_SAVDR_ROOT}/civic/civic_external_eval.csv')
OV = open(f'{_SAVDR_ROOT}/civic/civic_overlap.log').read()
TH = open(f'{_SAVDR_ROOT}/civic/therapy_match.log').read()
EVPY = open(f'{_SAVDR_ROOT}/civic/civic_eval.py').read()

ok = fail = 0
def chk(name, got, want):
    global ok, fail
    good = got == want
    ok, fail = ok + good, fail + (not good)
    print(f'  [{"PASS" if good else "FAIL"}] {name:58s} report={got!r:>22} source={want!r}')

def cell(ls, tgt, model, col='auc'):
    r = E[(E.label_set == ls) & (E.target == tgt) & (E.model == model)]
    assert len(r) == 1, (ls, tgt, model, len(r))
    return r.iloc[0][col]

ABL = 'ABL1_ENST00000318560'
TGT = {'ABL1': ABL, 'ALK': 'ALK', 'EGFR': 'EGFR', 'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}
MODELS = ['MultiPASS', 'pos_only', 'full_gbm', 'subst_lr']

# ---- 1. per-target table of set A ----------------------------------------------------------
print('\n1. Set A per-target AUC table')
sec = R.split('Per target (AUC):')[1].split('**This reproduces')[0]
rows = []
for l in sec.splitlines():
    if not l.strip().startswith('|'): continue
    f = [c.strip().replace('*', '').replace('†', '').strip() for c in l.strip('|').split('|')]
    if f[0] in TGT: rows.append(f)
assert len(rows) == 6, len(rows)
for f in rows:
    g = f[0]
    npos, nneg = (int(x) for x in f[1].split('/'))
    chk(f'A {g} n_pos', npos, int(cell('A_all_civic', TGT[g], 'MultiPASS', 'n_pos')))
    chk(f'A {g} n_neg', nneg, int(cell('A_all_civic', TGT[g], 'MultiPASS', 'n_neg')))
    for m, v in zip(MODELS, f[2:6]):
        chk(f'A {g} {m}', float(v), round(float(cell('A_all_civic', TGT[g], m)), 3))

# ---- 2. pooled tables ----------------------------------------------------------------------
for tag, ls, head in [('A', 'A_all_civic', '### A — all CIViC positives'),
                      ('B', 'B_discordant', '### B — discordant positives only')]:
    print(f'\n2{tag}. Set {tag} pooled table')
    blk = R.split(head)[1]
    m = re.search(r'pooled:\s*(\d+)\s*pos\s*/\s*(\d+)\s*neg', blk)
    chk(f'{tag} pooled n_pos', int(m.group(1)), int(cell(ls, 'POOLED', 'MultiPASS', 'n_pos')))
    chk(f'{tag} pooled n_neg', int(m.group(2)), int(cell(ls, 'POOLED', 'MultiPASS', 'n_neg')))
    body = blk.split('|---|---|---|---|---|---|')[1]
    seen = set()
    started = False
    for line in body.splitlines():
        if not line.strip().startswith('|'):
            if started: break
            continue
        started = True
        f = [c.strip().replace('*', '') for c in line.strip('|').split('|')]
        mod = f[0]
        if mod not in MODELS: continue
        seen.add(mod)
        chk(f'{tag} pooled {mod} auc', float(f[1]), round(float(cell(ls, 'POOLED', mod)), 3))
        lo, hi = (float(x) for x in f[2].replace('–', '-').split('-'))
        chk(f'{tag} pooled {mod} ci_lo', lo, round(float(cell(ls, 'POOLED', mod, 'ci_lo')), 3))
        chk(f'{tag} pooled {mod} ci_hi', hi, round(float(cell(ls, 'POOLED', mod, 'ci_hi')), 3))
        chk(f'{tag} pooled {mod} auprc', float(f[3]), round(float(cell(ls, 'POOLED', mod, 'auprc')), 3))
        if mod != 'MultiPASS':
            p = float(cell(ls, 'POOLED', mod, 'delong_p_vs_multipass'))
            h = float(cell(ls, 'POOLED', mod, 'delong_p_holm'))
            rp, rh = f[4], f[5]
            if rp.startswith('<'):
                chk(f'{tag} pooled {mod} p is tiny', p < float(rp[1:]), True)
                chk(f'{tag} pooled {mod} holm is tiny', h < float(rh.lstrip('<')), True)
            else:
                chk(f'{tag} pooled {mod} delong p', float(rp), round(p, 3))
                chk(f'{tag} pooled {mod} holm', float(rh), round(h, 2))
    chk(f'{tag} pooled table lists all 4 models', sorted(seen), sorted(MODELS))

# ---- 3. narrative claims in the A section --------------------------------------------------
print('\n3. Narrative claims, set A')
chk('A pooled: full_gbm > MultiPASS', cell('A_all_civic', 'POOLED', 'full_gbm') >
    cell('A_all_civic', 'POOLED', 'MultiPASS'), True)
chk('A pooled subst_lr Holm < 0.001',
    float(cell('A_all_civic', 'POOLED', 'subst_lr', 'delong_p_holm')) < 0.001, True)
chk('ABL1 A full_gbm raw p 0.0044 (text)', '0.0044' in R, True)
chk('ABL1 A full_gbm raw p value', round(float(cell('A_all_civic', ABL, 'full_gbm',
    'delong_p_vs_multipass')), 4), 0.0044)
chk('ABL1 A full_gbm Holm 0.150 (text)', '(Holm 0.150)' in R, True)
chk('ABL1 A full_gbm Holm value', round(float(cell('A_all_civic', ABL, 'full_gbm',
    'delong_p_holm')), 3), 0.15)
chk('ALK A pos_only beats MultiPASS', cell('A_all_civic', 'ALK', 'pos_only') >
    cell('A_all_civic', 'ALK', 'MultiPASS'), True)
chk('A pooled pos_only raw p = 0.101 (text)', 'p = 0.101 raw' in R, True)
chk('A pooled pos_only raw p value',
    round(float(cell('A_all_civic', 'POOLED', 'pos_only', 'delong_p_vs_multipass')), 3), 0.101)
chk('MultiPASS top point estimate on EGFR/KIT/SMO in A',
    all(cell('A_all_civic', t, 'MultiPASS') >= max(cell('A_all_civic', t, m) for m in MODELS)
        for t in ['EGFR', 'KIT', 'SMO']), True)
chk('SMO A spread of 3 non-chemistry models <= 0.006',
    round(max(cell('A_all_civic', 'SMO', m) for m in MODELS[:3]) -
          min(cell('A_all_civic', 'SMO', m) for m in MODELS[:3]), 3) <= 0.006, True)
chk('KIT A full_gbm raw p 0.0045',
    round(float(cell('A_all_civic', 'KIT', 'full_gbm', 'delong_p_vs_multipass')), 4), 0.0045)
chk('KIT A subst_lr raw p 0.0028',
    round(float(cell('A_all_civic', 'KIT', 'subst_lr', 'delong_p_vs_multipass')), 4), 0.0028)
chk('KIT A full_gbm Holm 0.150',
    round(float(cell('A_all_civic', 'KIT', 'full_gbm', 'delong_p_holm')), 3), 0.15)
chk('KIT A subst_lr Holm 0.098',
    round(float(cell('A_all_civic', 'KIT', 'subst_lr', 'delong_p_holm')), 3), 0.098)
# Superseded values may appear ONLY as an explicit before->after correction, never as a claim.
# The correction block and the two "X -> Y" annotations are the only licensed occurrences.
_nb = R.index('> **Numbers corrected 2026-10-08.**')
_ne = R.index('## How CIViC was read')
def only_as_correction(val):
    """val appears nowhere except inside the correction notice or beside a before->after arrow."""
    for m in re.finditer(re.escape(val), R):
        if _nb <= m.start() < _ne:                       # inside the correction notice
            continue
        ctx = R[max(0, m.start() - 300):m.end() + 60]
        if '→' in ctx[-120:] or 'earlier numbers' in ctx:   # annotated as superseded
            continue
        return False
    return True
chk('KIT 0.920 only as a correction', only_as_correction('0.920'), True)
chk('EGFR 0.686 only as a correction', only_as_correction('0.686'), True)
chk('stale neg-pool 4425 only as a correction', only_as_correction('4425'), True)
chk('stale neg-pool 4301 only as a correction', only_as_correction('4301'), True)
chk('stale B pooled 0.674 only as a correction', only_as_correction('0.674'), True)
chk('stale pooled A full_gbm 0.812 gone entirely', R.count('0.812') == 0, True)
chk('stale pooled B full_gbm 0.577 gone entirely', R.count('0.577') == 0, True)
chk('stale pooled B subst_lr 0.588 gone entirely', R.count('0.588') == 0, True)
chk('stale A per-target EGFR pos_only 0.567 gone', R.count('0.567') == 0, True)
chk('stale A per-target KIT subst_lr 0.669 gone', R.count('0.669') == 0, True)
chk('stale A per-target SMO subst_lr 0.075 gone', R.count('0.075') == 0, True)

# ---- 4. narrative claims in the B section --------------------------------------------------
print('\n4. Narrative claims, set B')
chk('B pooled MultiPASS CI includes 0.5',
    float(cell('B_discordant', 'POOLED', 'MultiPASS', 'ci_lo')) < 0.5 <
    float(cell('B_discordant', 'POOLED', 'MultiPASS', 'ci_hi')), True)
chk('B pooled: all Holm == 1.00', sorted({round(float(cell("B_discordant", "POOLED", m,
    "delong_p_holm")), 2) for m in MODELS[1:]}), [1.0])
chk('B pooled: MultiPASS is the top AUC', cell('B_discordant', 'POOLED', 'MultiPASS') ==
    max(cell('B_discordant', 'POOLED', m) for m in MODELS), True)
lo = min(cell('B_discordant', 'POOLED', m) for m in MODELS)
hi = max(cell('B_discordant', 'POOLED', m) for m in MODELS)
chk('B pooled range 0.566-0.626', (round(lo, 3), round(hi, 3)), (0.566, 0.626))
chk('B KIT pos_only raw p 0.0011',
    round(float(cell('B_discordant', 'KIT', 'pos_only', 'delong_p_vs_multipass')), 4), 0.0011)
chk('B KIT pos_only Holm 0.0396',
    round(float(cell('B_discordant', 'KIT', 'pos_only', 'delong_p_holm')), 4), 0.0396)
chk('B KIT MultiPASS 0.657 / pos_only 0.583',
    (round(float(cell('B_discordant', 'KIT', 'MultiPASS')), 3),
     round(float(cell('B_discordant', 'KIT', 'pos_only')), 3)), (0.657, 0.583))
chk('B KIT n_pos == 3', int(cell('B_discordant', 'KIT', 'MultiPASS', 'n_pos')), 3)
chk('B ABL1 full_gbm 0.970 on 1 positive',
    (round(float(cell('B_discordant', ABL, 'full_gbm')), 3),
     round(float(cell('B_discordant', ABL, 'full_gbm', 'ci_lo')), 3),
     round(float(cell('B_discordant', ABL, 'full_gbm', 'ci_hi')), 3),
     int(cell('B_discordant', ABL, 'full_gbm', 'n_pos'))), (0.97, 0.957, 0.982, 1))
chk('B per-target n_pos range is 1-7',
    (int(E[(E.label_set == 'B_discordant') & (E.target != 'POOLED')].n_pos.min()),
     int(E[(E.label_set == 'B_discordant') & (E.target != 'POOLED')].n_pos.max())), (1, 7))
chk('NT5C2 absent from set B (0 discordant)',
    len(E[(E.label_set == 'B_discordant') & (E.target == 'NT5C2')]), 0)

# ---- 5. internal consistency of the CSV ----------------------------------------------------
print('\n5. CSV internal consistency')
for ls in ['A_all_civic', 'B_discordant']:
    per = E[(E.label_set == ls) & (E.target != 'POOLED') & (E.model == 'MultiPASS')]
    chk(f'{ls}: pooled n_pos == sum of targets', int(per.n_pos.sum()),
        int(cell(ls, 'POOLED', 'MultiPASS', 'n_pos')))
    chk(f'{ls}: pooled n_neg == sum of targets', int(per.n_neg.sum()),
        int(cell(ls, 'POOLED', 'MultiPASS', 'n_neg')))
chk('A n_pos == 84 == CIViC SAVs present in our set', int(cell('A_all_civic', 'POOLED',
    'MultiPASS', 'n_pos')), 84)
chk('B n_pos == 14 == discordant calls', int(cell('B_discordant', 'POOLED',
    'MultiPASS', 'n_pos')), 14)
chk('A neg pool - B neg pool == NT5C2 negatives',
    int(cell('A_all_civic', 'POOLED', 'MultiPASS', 'n_neg')) -
    int(cell('B_discordant', 'POOLED', 'MultiPASS', 'n_neg')),
    int(cell('A_all_civic', 'NT5C2', 'MultiPASS', 'n_neg')))
chk('39 DeLong tests in the file', int(E.delong_p_vs_multipass.notna().sum()), 39)
chk('report says 39 tests', '39 tests' in ' '.join(R.split()), True)

# ---- 6. overlap table and discordance list (should be unaffected) --------------------------
print('\n6. Overlap table vs civic_overlap.log')
chk('90 distinct CIViC resistance SAVs', '90' in OV.splitlines()[0], True)
chk('0 wild-type mismatches', 'mismatches against our reference FASTA: 0' in OV, True)
for g, tot, pos, neg, absent, cont in [('ABL1', 45, 39, 1, 5, 33), ('ALK', 13, 11, 2, 0, 0),
                                       ('EGFR', 11, 3, 7, 1, 0), ('KIT', 14, 11, 3, 0, 5),
                                       ('NT5C2', 5, 5, 0, 0, 5), ('SMO', 2, 1, 1, 0, 0)]:
    line = [l for l in OV.splitlines() if l.strip().startswith(g)][0]
    nums = [int(x) for x in line.split()[1:6]]
    chk(f'overlap {g}', [tot, pos, neg, absent, cont], nums)
    rl = [l for l in R.splitlines() if l.startswith(f'| {g} |')][0]
    rn = [int(c.strip()) for c in rl.strip('|').split('|')[1:6]]
    chk(f'report row {g}', rn, nums)
tot = [l for l in R.splitlines() if l.startswith('| **total**')][0]
chk('report totals row', [int(c.strip().strip('*')) for c in tot.strip('|').split('|')[1:6]],
    [90, 70, 14, 6, 43])
disc = [l.split()[0] + ' ' + l.split()[1] for l in OV.split('call NOT resistant:')[1]
        .split('In CIViC but absent')[0].splitlines() if l.strip()]
chk('14 discordant calls in the log', len(disc), 14)
chk('all 14 appear in the report', all(f'{g} {v}' in R.replace('   ', ' ').replace('  ', ' ')
    or v in R for g, v in (d.split() for d in disc)), True)

print('\n7. Therapy-matched section vs therapy_match.log')
chk('9 of 14 name a drug we label',
    'of 14 discordant calls, 9 name a drug we do label and 5 name a drug outside' in TH, True)
chk('report says 9 and 5', '**9 name a drug we do label** and only 5 fall outside' in R, True)
chk('3 antibody-only calls', 'antibody-only calls among them: 3' in TH, True)
chk('70 agreeing, 67 share a drug', '70 such variants; drug sets overlap for 67 of them' in TH, True)
chk('report says 67 share a drug', '**67 share at least one drug**' in R, True)
chk('axitinib trio in log', all(v in TH for v in ['p.F311I', 'p.F317R', 'p.M343T']), True)
chk('axitinib trio in report', 'F311I, F317R and\nM343T' in R or
    'F311I, F317R and M343T' in R, True)

print('\n8. The `interpretable` column claim')
chk('civic_eval.py does not write an interpretable column', 'interpret' in EVPY, False)
chk('CSV has no interpretable column', 'interpretable' in E.columns, False)
chk('report now says there is no such column', 'There is **no** `interpretable` column' in R, True)

print(f'\n{"=" * 100}\n{ok} passed, {fail} failed')
sys.exit(1 if fail else 0)
