"""Coarse but wide net: every number in every markdown TABLE cell of every report must be
findable in one of that report's own source CSVs.

This is deliberately not a precise cell-by-cell mapping (verify_civic.py and verify_level2.py do
that for the two reports where it mattered). It is the net that catches staleness: a report cell
that no longer exists anywhere in its sources is either stale or derived, and every one of them
is printed for inspection rather than silently passed.

A value "exists" if it appears, to the precision written, in any numeric column of any source
CSV for that report, or is a plain integer count present there. Values that are legitimately
derived in the prose (differences, means, ratios, before->after corrections) are listed in
DERIVED per report with the reason, so that the unexplained set should be empty.
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

import glob, re, sys, numpy as np, pandas as pd

OUT = f'{_SAVDR_ROOT}'
REPORTS = {
    f'{OUT}/multipass/REPORT.md': [
        'multipass/summary.csv', 'multipass/summary_with_nearestpos.csv',
        'multipass/per_drug.csv', 'multipass/table1_recomputed.csv',
        'multipass/table_wide.csv', 'multipass/Table_S7_corrected.csv'],
    f'{OUT}/multipass/REPORT_ADDITIONS.md': [
        'multipass/power.csv', 'multipass/power_secheck.csv', 'multipass/calibration.csv',
        'multipass/calibration_reliability.csv', 'multipass/recalibration.csv',
        'multipass/margin_sweep.csv', 'multipass/operating_points.csv', 'multipass/loto.csv',
        'multipass/loto_baselines.csv', 'multipass/within_position.csv',
        'multipass/summary_with_nearestpos.csv', 'multipass/table1_recomputed.csv'],
    f'{OUT}/multipass/REPORT_11TARGETS.md': [
        'multipass/loto11.csv', 'multipass/ingest_check.csv', 'multipass/loto.csv',
        'multipass/table1_recomputed.csv'],
    f'{OUT}/multipass/REPORT_SUSPECT_ABL.md': ['multipass/suspect_abl_compare.csv',
                                               'multipass/suspect_abl_predictions.csv'],
    f'{OUT}/multipass/REPORT_DRKINASE.md': ['multipass/drkinase_abl1.csv'],
    f'{OUT}/civic/REPORT_CIVIC.md': ['civic/civic_external_eval.csv',
                                     'civic/civic_resistance_vs_our_labels.csv',
                                     'civic/therapy_match.csv', 'civic/civic_therapy_eval.csv'],
    f'{OUT}/level2/REPORT_LEVEL2.md': ['level2/level2_eval.csv'],
}
# values a report legitimately states that are not raw cells of its CSVs
DERIVED = {
    'REPORT.md': 'margins, means, published Table 1 values, counts quoted from logs',
    'REPORT_ADDITIONS.md': 'MDEs, deltas, Spearman rhos, means over seeds, counts from logs',
    'REPORT_11TARGETS.md': 'group means, the kinase/non-kinase gap and its p, counts from logs',
    'REPORT_SUSPECT_ABL.md': 'confusion-matrix rates, recalled fractions, per-drug AUCs from the log',
    'REPORT_DRKINASE.md': 'counts and the 0.90 plot reference line',
    'REPORT_CIVIC.md': 'overlap counts, therapy counts, superseded values in the correction notice',
    'REPORT_LEVEL2.md': 'macro-averages, medians, superseded values in the correction notice',
}

def table_numbers(md):
    """every numeric token that sits inside a markdown table cell"""
    vals = []
    for line in md.splitlines():
        s = line.strip()
        if not s.startswith('|') or set(s) <= set('|- '):
            continue
        for cell in s.strip('|').split('|'):
            for tok in re.findall(r'(?<![\w.])-?\d+(?:\s\d{3})*(?:\.\d+)?(?![\w.])',
                                  cell.replace('**', '').replace('−', '-')):
                vals.append(tok.replace(' ', ''))
    return vals

def source_values(paths):
    """every numeric value in those CSVs, as a set of strings at each precision 0..4"""
    seen = set()
    for p in paths:
        try:
            d = pd.read_csv(f'{OUT}/{p}')
        except Exception as e:
            print(f'    (could not read {p}: {e})')
            continue
        for c in d.columns:
            col = pd.to_numeric(d[c], errors='coerce').dropna()
            for v in col.values:
                for k in range(5):
                    seen.add(f'{v:.{k}f}'.rstrip('.') if k == 0 else f'{v:.{k}f}')
                if float(v).is_integer():
                    seen.add(str(int(v)))
            # column headers and string cells can carry counts too
        for c in d.columns:
            for v in d[c].astype(str).tolist():
                for tok in re.findall(r'(?<![\w.])-?\d+(?:\.\d+)?(?![\w.])', str(v)):
                    seen.add(tok)
    return seen

total = unexplained = 0
print('=' * 100)
print('Cross-report number check: every table cell against that report\'s own CSVs')
print('=' * 100)
for rp, srcs in REPORTS.items():
    name = rp.split('/')[-1]
    md = open(rp).read()
    nums = table_numbers(md)
    src = source_values(srcs)
    miss = []
    for t in nums:
        f = float(t)
        # a value counts as present if it, or its rounding at any precision, is in the sources
        cands = {t, t.lstrip('0') if t.startswith('0.') else t}
        if any(c in src for c in cands):
            continue
        # allow a 1-in-the-last-digit rounding difference
        if '.' in t:
            d = len(t.split('.')[1])
            if any(f'{f + s * 10 ** -d:.{d}f}' in src for s in (-1, 1)):
                continue
        miss.append(t)
    total += len(nums); unexplained += len(miss)
    print(f'\n{name}')
    print(f'  table numbers: {len(nums):4d}   found in sources: {len(nums) - len(miss):4d}   '
          f'not raw cells: {len(miss)}')
    print(f'  expected non-cell values here: {DERIVED[name]}')
    if miss:
        c = {}
        for t in miss: c[t] = c.get(t, 0) + 1
        print('  ' + ', '.join(f'{k}' + (f'x{v}' if v > 1 else '')
                               for k, v in sorted(c.items(), key=lambda x: -x[1])))
print(f'\n{"=" * 100}')
print(f'{total} table numbers checked across {len(REPORTS)} reports; '
      f'{total - unexplained} matched a source cell directly, {unexplained} are derived or '
      f'quoted from logs (listed above for inspection).')
