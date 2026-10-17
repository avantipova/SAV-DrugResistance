"""Independent re-verification of everything computed in this round.

Each check recomputes a published number by a different route than the script that produced it,
or tests an invariant that would be violated if the producing script were wrong. Anything that
is NOT checked here is listed at the end, so the scope of this verification is explicit.
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

import json, glob, os, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats

O = f'{_SAVDR_ROOT}/multipass'
ok, bad = [], []


def chk(name, cond, detail=''):
    (ok if cond else bad).append(f'{name}: {detail}')
    print(('  PASS  ' if cond else '  FAIL  ') + name + ('   ' + detail if detail else ''))


print('=' * 100)
print('1. POWER / MINIMUM DETECTABLE EFFECT')
print('=' * 100)
p = pd.read_csv(f'{O}/power.csv')
se = pd.read_csv(f'{O}/power_secheck.csv')
# the MDE must equal 2.8016 * se, and excl_95 must equal delta + 1.6449 * se
z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)
chk('mde_80 == (z_0.975+z_0.80)*se', np.allclose(p.mde_80, np.round(z * p.se, 3), atol=1.5e-3),
    f'max dev {np.abs(p.mde_80 - z * p.se).max():.4f}')
chk('excl_95 == delta + 1.6449*se',
    np.allclose(p.excl_95, np.round(p.delta_obs + stats.norm.ppf(0.95) * p.se, 3), atol=1.5e-3),
    f'max dev {np.abs(p.excl_95 - (p.delta_obs + stats.norm.ppf(0.95) * p.se)).max():.4f}')
chk('delta_obs == auc_multipass - auc_baseline',
    np.allclose(p.delta_obs, np.round(p.auc_multipass - p.auc_baseline, 3), atol=1.5e-3))
chk('bootstrap SE agrees with DeLong SE within 10%',
    bool(((se.ratio > 0.9) & (se.ratio < 1.1)).all()),
    f'ratio range {se.ratio.min():.3f}-{se.ratio.max():.3f}')
chk('no cell has |delta_obs| > mde_80 for pos_only',
    int((p[p.baseline == "pos_only"].delta_obs.abs() > p[p.baseline == "pos_only"].mde_80).sum()) == 0)

_FOLDS = pd.read_csv(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_assignments.csv')

print('=' * 100)
print('2. nearest_pos BASELINE -- recomputed from scratch, independent implementation')
print('=' * 100)
sm = pd.read_csv(f'{O}/summary_with_nearestpos.csv')
CFG = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
       'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
for (name, (L, M)), region in [((n, c), r) for n, c in CFG.items() for r in (False, True)]:
    tag = 'fold_%s_%d_%d%s_s' % (name, L, M, '_region' if region else '')
    fs = sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/' + tag + '*/scores.csv'))
    if not fs: continue
    per = []
    for f in fs:
        d = pd.read_csv(f)
        y, pos = d.y.values, d.pos.values
        s = np.zeros(len(y))
        # The fold assignment is READ, not rebuilt: scikit-learn's StratifiedGroupKFold does
        # not give the same split across versions, so rebuilding it here made this block
        # disagree with the published numbers on a different scikit-learn (observed: 12 of
        # these 12 checks failed on 1.7.2 while passing on 1.8.0). The assignment that
        # produced the reported values is stored in inputs/cv_folds/fold_assignments.csv.
        fa = _FOLDS[_FOLDS.dir == _os.path.basename(_os.path.dirname(f))]
        assert list(fa.sort_values('row').aa) == list(d.aa), \
            f'stored fold assignment does not line up with {f}'
        for k in sorted(fa.fold.unique()):
            te = fa[fa.fold == k].row.values
            tr = fa[fa.fold != k].row.values
            # independent route: loop instead of the broadcast min used in extras.py
            hot = sorted(set(pos[tr][y[tr] == 1]))
            s[te] = [-min(abs(pp - h) for h in hot) if hot else 0.0 for pp in pos[te]]
        per.append(stats.rankdata(s) / len(s))
    v = np.mean(per, 0)
    want = sm[(sm.target == name) & (sm.setting == ('region' if region else 'full')) &
              (sm.model == 'nearest_pos')].auc_pooled.iloc[0]
    got = round(roc_auc_score(y, v), 3)
    chk(f'nearest_pos {name} {"region" if region else "full"}', abs(got - want) < 1e-9,
        f'reported {want}  recomputed {got}')

print('=' * 100)
print('3. CALIBRATION / RECALIBRATION')
print('=' * 100)
cal = pd.read_csv(f'{O}/calibration.csv')
rec = pd.read_csv(f'{O}/recalibration.csv')
chk('brier_baserate == prevalence*(1-prevalence)',
    np.allclose(cal.brier_baserate, np.round(cal.prevalence * (1 - cal.prevalence), 4), atol=1.5e-4),
    'the Brier of a constant base-rate forecast is p(1-p)')
chk('raw Pa worse than base rate in 10 of 12', int((cal.brier > cal.brier_baserate).sum()) == 10)
chk('Platt beats base rate in 12 of 12', int((rec.brier_platt < rec.brier_baserate).sum()) == 12)
chk('both recalibration maps monotone within every fold',
    bool(rec.monotone_within_fold_iso.all() and rec.monotone_within_fold_platt.all()))
chk('Platt leaves within-fold AUC essentially unchanged',
    float((rec.auc_withinfold_platt - rec.auc_withinfold_raw).abs().max()) < 0.002,
    f'max |delta| {float((rec.auc_withinfold_platt - rec.auc_withinfold_raw).abs().max()):.5f}')
chk('reliability bins sum to the full variant count',
    True, 'checked per cell below')
rel = pd.read_csv(f'{O}/calibration_reliability.csv')
for (t, s), g in rel.groupby(['target', 'setting']):
    row = cal[(cal.target == t) & (cal.setting == s)].iloc[0]
    chk(f'  bins sum {t} {s}', int(g.n.sum()) == int(row.n_pos + row.n_neg),
        f'{int(g.n.sum())} vs {int(row.n_pos + row.n_neg)}')

print('=' * 100)
print('4. MARGIN SWEEP')
print('=' * 100)
ms = pd.read_csv(f'{O}/margin_sweep.csv')
chk('n_pos constant across margins within a target',
    bool(ms.groupby('target').n_pos.nunique().eq(1).all()),
    'the margin only admits negatives, never changes the positive set')
mono = []
for t, g in ms[ms.model == 'multipass_full_trained'].groupby('target'):
    g = g[g.margin != 'full']
    g = g.assign(mg=g.margin.astype(int)).sort_values('mg')
    mono.append(bool((g.n_neg.diff().dropna() > 0).all()))
chk('n_neg strictly increases with margin', all(mono))
chk('"full" row has the largest n_neg in every target',
    bool(all(ms[(ms.target == t) & (ms.margin == 'full')].n_neg.iloc[0] >
             ms[(ms.target == t) & (ms.margin != 'full')].n_neg.max() for t in ms.target.unique())))
sp = {}
for t, g in ms[ms.model == 'multipass_full_trained'].groupby('target'):
    ORD = {'5': 0, '10': 1, '20': 2, '30': 3, '50': 4, 'full': 5}
    g = g.assign(o=g.margin.astype(str).map(ORD)).sort_values('o')
    sp[t] = stats.spearmanr(g.o, g.auc).statistic
chk('MultiPASS AUC rises with margin in all six targets', all(v > 0 for v in sp.values()),
    'rho ' + ', '.join(f'{t.split("_")[0]}={v:+.2f}' for t, v in sp.items()))
sl = {}
for t, g in ms[ms.model == 'subst_lr'].groupby('target'):
    ORD = {'5': 0, '10': 1, '20': 2, '30': 3, '50': 4, 'full': 5}
    g = g.assign(o=g.margin.astype(str).map(ORD)).sort_values('o')
    sl[t] = stats.spearmanr(g.o, g.auc).statistic
chk('the position-free control (subst_lr) shows no margin trend',
    abs(float(np.mean(list(sl.values())))) < 0.2, f'mean rho {float(np.mean(list(sl.values()))):+.3f}')

print('=' * 100)
print('5. LEAVE-ONE-TARGET-OUT')
print('=' * 100)
lo = pd.read_csv(f'{O}/loto.csv')
for lvl in (9, 4):
    for held, (L, M) in CFG.items():
        f = f'{_SAVDR_ROOT}/inputs/loto_bases/loto_{held}_mna{lvl}/scores.csv'
        if not os.path.exists(f): continue
        d = pd.read_csv(f)
        main = pd.read_csv(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{held}_{L}_{M}_s0/scores.csv')
        chk(f'  LOTO variant set identical to main table, {held} MNA{lvl}',
            list(d.aa) == list(main.aa) and list(d.y) == list(main.y))
        want = lo[(lo.held_out == held) & (lo.model == f'MultiPASS_loto_mna{lvl}')].auc.iloc[0]
        got = round(roc_auc_score(d.y, np.nan_to_num(d.multipass_pa.values, nan=0.0)), 3)
        chk(f'  LOTO AUC {held} MNA{lvl}', abs(got - want) < 1e-9, f'reported {want} got {got}')
# the decisive invariant: the held-out protein contributed nothing to its own base
import re
for held in CFG:
    for lvl in (9, 4):
        w = f'{_SAVDR_ROOT}/inputs/loto_bases/loto_{held}_mna{lvl}'
        if not os.path.exists(f'{w}/B.HST'): continue
        txt = open(f'{w}/B.HST', errors='ignore').read()
        n = int(re.search(r'(\d+)\s+Substances', txt).group(1))
        exp = sum(len(pd.read_csv(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_{c[0]}_{c[1]}_s0/scores.csv'))
                  for t, c in CFG.items() if t != held)
        chk(f'  base for held-out {held} MNA{lvl} trained on exactly the other five',
            n == exp, f'{n} substances, expected {exp}')

print('=' * 100)
print('6. NESTED CONFIGURATION SEARCH')
print('=' * 100)
for t in ('NT5C2', 'ABL1_ENST00000318560', 'SMO'):
    f = f'{O}/nested_{t}.json'
    if not os.path.exists(f): continue
    j = json.load(open(f))
    if not j.get('nested_computable', True): continue
    chk(f'  {t}: optimism == flat_max - nested',
        abs(j['optimism'] - (j['flat_max_auc'] - j['nested_auc'])) < 1e-9,
        f'{j["flat_max_auc"]} - {j["nested_auc"]} = {j["optimism"]}')
    chk(f'  {t}: flat_max is really the max over the grid',
        abs(j['flat_max_auc'] - max(j['flat_auc'].values())) < 1e-9)
    chk(f'  {t}: flat_best_config matches that max',
        max(j['flat_auc'], key=j['flat_auc'].get) == j['flat_best_config'])
    # nested AUC must be reproducible from the per-fold checkpoints of the picked configs
    root = f'{_SAVDR_ROOT}/inputs/nested_checkpoints/nested_{t}'
    outer_sets = {}
    allsc = {}
    for k, cfg in j['nested_picked'].items():
        L, M = re.match(r'L(\d+)_MNA(\d+)', cfg).groups()
        ck = f'{root}/flat_L{L}_M{M}_f{k}.json'
        if not os.path.exists(ck): continue
        s = json.load(open(ck))
        outer_sets[k] = set(s)
        allsc.update(s)
    chk(f'  {t}: the five held-out folds are disjoint',
        sum(len(v) for v in outer_sets.values()) == len(set().union(*outer_sets.values())))
    chk(f'  {t}: held-out folds cover the common variant set',
        len(set().union(*outer_sets.values())) == j['n_common'],
        f'{len(set().union(*outer_sets.values()))} vs {j["n_common"]}')


print('=' * 100)
print('7. REPORT_ADDITIONS.md TEXT AGAINST THE CSVs  (added 2026-10-08)')
print('=' * 100)
A = open(f'{O}/REPORT_ADDITIONS.md').read()
rec7 = pd.read_csv(f'{O}/recalibration.csv')
cal7 = pd.read_csv(f'{O}/calibration.csv')
op7 = pd.read_csv(f'{O}/operating_points.csv')
op7 = op7[op7.setting == 'full']

for lab, col, want in [('raw Pa', 'brier_raw_pa', '0.0938'), ('Platt', 'brier_platt', '0.0803'),
                       ('isotonic', 'brier_isotonic', '0.0787'),
                       ('base-rate forecast', 'brier_baserate', '0.0832')]:
    got = f'{rec7[col].mean():.4f}'
    chk(f'  mean Brier, {lab}', got == want, f'csv {got}, report {want}')
    chk(f'  report row "| {lab} |" carries it', f'| {lab} | {want} |' in A
        or f'| {lab} | **{want}** |' in A)
chk('  "beats base rate" column: 2 / 10 / 12 of 12',
    (int(rec7.raw_beats_baserate.sum()), int(rec7.iso_beats_baserate.sum()),
     int((rec7.brier_platt < rec7.brier_baserate).sum())) == (2, 10, 12),
    f'{int(rec7.raw_beats_baserate.sum())} / {int(rec7.iso_beats_baserate.sum())} / '
    f'{int((rec7.brier_platt < rec7.brier_baserate).sum())}')
dmax = float((rec7.auc_withinfold_raw - rec7.auc_withinfold_iso).max())
dmean = float((rec7.auc_withinfold_iso - rec7.auc_withinfold_raw).mean())
chk('  isotonic within-fold loss: worst cell -0.053, mean -0.013',
    (round(dmax, 3), round(dmean, 3)) == (0.053, -0.013), f'worst {dmax:.4f}, mean {dmean:.4f}')
_wi = (rec7.auc_withinfold_raw - rec7.auc_withinfold_iso).idxmax()
chk('  the worst isotonic cell is NT5C2 full',
    (rec7.loc[_wi, 'target'], rec7.loc[_wi, 'setting']) == ('NT5C2', 'full'),
    f'{rec7.loc[_wi, "target"]} {rec7.loc[_wi, "setting"]}')
chk('  report names NT5C2 full as that cell', 'NT5C2 `full`**' in A)
chk('  report states the -0.053 is the worst cell, not the mean',
    'worst single cell' in A and 'mean **\u22120.013**' in A)
pm = float((rec7.auc_withinfold_platt - rec7.auc_withinfold_raw).mean())
pa_ = float((rec7.auc_withinfold_platt - rec7.auc_withinfold_raw).abs().mean())
pf = float((rec7[rec7.setting == 'full'].auc_withinfold_platt
            - rec7[rec7.setting == 'full'].auc_withinfold_raw).abs().max())
chk('  Platt within-fold: mean +0.00008, |mean| 0.00012, 0 in all full cells',
    (f'{pm:.5f}', f'{pa_:.5f}', pf == 0.0) == ('0.00008', '0.00012', True),
    f'{pm:.5f} / {pa_:.5f} / full-max {pf}')
chk('  pooled isotonic AUC loss reaches 0.098',
    round(float((rec7.auc_raw_pa - rec7.auc_isotonic).max()), 3) == 0.098)
chk('  raw Pa overestimates by up to 5.9x, worst cell ALK full',
    round(float((cal7.pa_mean / cal7.prevalence).max()), 1) == 5.9
    and cal7.loc[(cal7.pa_mean / cal7.prevalence).idxmax(), 'target'] == 'ALK')

OPROW = {'MultiPASS': ('0.467', '0.400', '0.465', '0.253'),
         'full_gbm': ('0.317', '0.400', '0.467', '0.194'),
         'nearest_pos': ('0.317', '0.275', '0.480', '0.084'),
         'pos_only': ('0.283', '0.308', '0.413', '0.128'),
         'subst_lr': ('0.133', '0.108', '0.1205', '0.015')}
for mod, want in OPROW.items():
    g = op7[op7.model == mod]
    got = tuple(f'{g[c].mean():.{len(w.split(".")[1])}f}' for c, w in zip(
        ['prec_top10', 'prec_top20', 'sens_at_spec95', 'sens_at_spec99'], want))
    chk(f'  operating points, {mod} (mean over the 6 full-setting targets)', got == want,
        f'csv {got}, report {want}')
    chk(f'  report row for {mod}', f'| {mod} | ' + ' | '.join(want) + ' |' in A)
chk('  operating-points means are over 6 targets', int(op7.target.nunique()) == 6)


print('=' * 100)
print('8. CONFIGURATION-SELECTION OPTIMISM ACROSS ALL ELEVEN TARGETS  (added 2026-10-10)')
print('=' * 100)
S = pd.read_csv(f'{O}/nested_summary.csv')
chk('  all eleven targets have a result', len(S) == 11, f'{len(S)} rows')
SZ = {}
for f in os.listdir(f'{_SAVDR_ROOT}/inputs/variants'):
    if f.endswith('.csv') and not f.startswith('_'):
        dd = pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/{f}')
        SZ[f[:-4]] = (len(dd), int(dd.y.sum()))
for _, r in S.iterrows():
    t = r.target
    chk(f'  {t[:12]}: n and n_pos match the variant table',
        (int(r.n), int(r.n_pos)) == SZ[t], f'table says {SZ[t]}')
    f1, f2 = f'{O}/nested_{t}.json', f'{O}/nested_{t}_flat.json'
    j = json.load(open(f1)) if os.path.exists(f1) else json.load(open(f2))
    fa = j.get('flat_auc') or {}
    chk(f'  {t[:12]}: trainable + untrainable == 6 configs',
        len(fa) + len(j.get('flat_not_trainable') or []) == 6,
        f'{len(fa)} + {len(j.get("flat_not_trainable") or [])}')
    if fa:
        v = list(fa.values())
        chk(f'  {t[:12]}: grid min/median/max match the summary',
            (round(min(v), 4), round(float(np.median(v)), 4), round(max(v), 4))
            == (round(r.grid_min, 4), round(r.grid_median, 4), round(r.grid_max, 4)))
        chk(f'  {t[:12]}: best_config is the argmax of the grid',
            (j.get('flat_best_config') or max(fa, key=fa.get)) == max(fa, key=fa.get))
    if not np.isnan(r.optimism):
        # The stored optimism is round(unrounded_flat - unrounded_nested, 4), which is not
        # always the difference of the two rounded columns: subtracting two independently
        # 4-dp values can be off by 1e-4 in either direction. ALK is the case in point --
        # 0.9287500 - 0.9277462 = 0.0010038 -> 0.0010, while the printed columns give
        # 0.9288 - 0.9277 = 0.0011, because the two roundings move in opposite directions.
        # The tolerance is the size of that artefact and no larger, so a genuine
        # disagreement in the fourth decimal still fails. verify_optimism.py uses the same
        # bound against the same checkpoints.
        gap = abs(round(float(r.grid_max) - float(r.nested_auc), 4) - round(float(r.optimism), 4))
        chk(f'  {t[:12]}: optimism == flat_max - nested, within the double-rounding artefact',
            gap <= 2.0001e-4,
            f'{r.grid_max} - {r.nested_auc} vs stored {r.optimism}, gap {gap:.6f}')
        chk(f'  {t[:12]}: every picked config is in the grid',
            all(c in j['grid'] for c in j['nested_picked'].values()))
        chk(f'  {t[:12]}: five outer folds each picked a config', len(j['nested_picked']) == 5)
for t in ('BTK_ENST00000308731', 'FLT3'):
    j = json.load(open(f'{O}/nested_{t}_flat.json'))
    chk(f'  {t[:12]}: no configuration trainable at all', len(j['flat_auc']) == 0)
    chk(f'  {t[:12]}: all six listed as not trainable', len(j['flat_not_trainable']) == 6)
    log = open(f'{_SAVDR_ROOT}/results/data/nested_{t}.log').read()
    chk(f'  {t[:12]}: the log records six refusals', log.count('declined to train') == 6,
        f'{log.count("declined to train")} refusals')
sub = S[S.spread.notna() & (S.n_pos > 10)]
chk('  the target with most positives has the narrowest grid spread',
    sub.loc[sub.n_pos.idxmax(), 'target'] == sub.loc[sub.spread.idxmin(), 'target'],
    f'most positives: {sub.loc[sub.n_pos.idxmax(), "target"]}, '
    f'narrowest: {sub.loc[sub.spread.idxmin(), "target"]}')
chk('  measured optimism exists for exactly six targets', int(S.optimism.notna().sum()) == 6,
    f'{int(S.optimism.notna().sum())} targets')
chk('  and those six are ABL1, ALK, EGFR, KIT, NT5C2, SMO',
    sorted(S[S.optimism.notna()].target.str.split('_').str[0])
    == ['ABL1', 'ALK', 'EGFR', 'KIT', 'NT5C2', 'SMO'],
    ', '.join(sorted(S[S.optimism.notna()].target.str.split('_').str[0])))
chk('  optimism values are 0.0010 / 0.0019 / 0.0060 / 0.0194 / 0.0363 / 0.1337',
    sorted(round(float(x), 4) for x in S.optimism.dropna())
    == [0.0010, 0.0019, 0.0060, 0.0194, 0.0363, 0.1337],
    ' / '.join(f'{x:.4f}' for x in sorted(S.optimism.dropna())))
chk('  the three with no inner-fold fit are MET, ESR1, IDH1',
    sorted(S[S.stage.str.startswith('all 6')].target.str.split('_').str[0])
    == ['ESR1', 'IDH1', 'MET'])
chk('  the two where nothing trains at all are BTK, FLT3',
    sorted(S[S.stage.str.startswith('no configuration')].target.str.split('_').str[0])
    == ['BTK', 'FLT3'])
# The drift guard. nested_summary.csv was maintained by hand until 2026-10-17 and had gone
# stale: ALK, EGFR and KIT were recorded as "flat only" after their nested runs had finished,
# so the table disagreed with its own checkpoints and with Table S8. It is now derived, and
# this check fails if the committed file and the checkpoints ever part company again.
import io as _io, csv as _csv, importlib.util as _ilu
_spec = _ilu.spec_from_file_location('_mns', f'{O}/make_nested_summary.py')
_mns = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_mns)
_buf = _io.StringIO()
_w = _csv.DictWriter(_buf, _mns.COLS); _w.writeheader()
_w.writerows([_mns.row(t) for t in _mns.ORDER])
chk('  nested_summary.csv is exactly what the checkpoints derive',
    _buf.getvalue().replace('\r\n', '\n') ==
    open(f'{O}/nested_summary.csv').read().replace('\r\n', '\n'),
    'run make_nested_summary.py')
chk('  the probe file covers the five small-positive targets',
    len(json.load(open(f'{O}/nested_probe.json'))) == 5)

print('=' * 100)
print(f'{len(ok)} checks passed, {len(bad)} failed')
if bad:
    print('FAILURES:')
    for b in bad: print('  ' + b)
print("""
NOT verified by this script:
  * MultiPASS's own internals -- a closed Win32 binary.
  * That MNA 9 and MNA 4 are the right levels for a cross-protein base. Both were run and the
    difference is not Holm-significant anywhere, but the grid is two points, not a sweep.
  * The nested configuration search runs on a 6-point grid and one outer seed. It covers the six
    targets where it is computable; on MET, ESR1 and IDH1 every configuration declines to train on
    the inner folds, and on BTK and FLT3 none trains even on the flat grid, so for those five the
    optimism is unquantified rather than small.
  * No power curve is claimed anywhere; only the SE of the AUC difference was validated.
  * The CIViC therapy matching depends on drug-name normalisation; the class 'tyrosine kinase
    inhibitor - NS' has no CIViC member list and was excluded rather than guessed.
""")
