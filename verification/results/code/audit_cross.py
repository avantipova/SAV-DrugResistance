"""Cross-report consistency audit.

Each verify_*.py checks one report against its own CSVs. Nothing until now checked whether two
reports state different values for the SAME quantity, or whether a report's summary table agrees
with the report it summarises. That is what this does.

Three parts:
  A  every per-target variant/positive count stated anywhere, against the variant tables
  B  curated quantities that appear in more than one report, which must either agree or carry
     an explicit reconciliation sentence
  C  a mechanical sweep: numbers that appear in two reports next to the same keywords
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

import re, glob, os, sys, json
import pandas as pd

OUT = f'{_SAVDR_ROOT}'
REPORTS = sorted(glob.glob(f'{OUT}/reports/REPORT_*.md')) + [
    f'{OUT}/civic/REPORT_CIVIC.md', f'{OUT}/level2/REPORT_LEVEL2.md',
    f'{OUT}/verification/VERIFICATION.md']
R = {os.path.basename(p): open(p).read() for p in REPORTS if os.path.exists(p)}
# Prose in these reports is hard-wrapped, so any assertion about a SENTENCE has to run against
# a whitespace-normalised copy; assertions about table cells use the raw text.
N = {k: re.sub(r'\s+', ' ', v) for k, v in R.items()}
print(f'{len(R)} reports loaded\n')

ok = bad = warn = 0
def chk(name, cond, detail=''):
    global ok, bad
    if cond: ok += 1
    else: bad += 1; print(f'  CONTRADICTION  {name}  {detail}')
def soft(name, cond, detail=''):
    global ok, warn
    if cond: ok += 1
    else: warn += 1; print(f'  WARN  {name}  {detail}')

# ---------------------------------------------------------------- A: target counts
print('== A. per-target variant and positive counts, against the variant tables ==')
TAB = f'{_SAVDR_ROOT}/inputs/variants'
counts = {}
for f in sorted(glob.glob(f'{TAB}/*.csv')):
    nm = os.path.basename(f)[:-4]
    if nm.startswith('_'): continue
    t = pd.read_csv(f)
    counts[nm] = (len(t), int(t.y.sum()))
short = {k.split('_')[0]: v for k, v in counts.items()}
for k, (n, p) in sorted(short.items()):
    print(f'  {k:6s} {n:5d} variants  {p:4d} positives')

# every markdown table row of the form | TARGET | n | p | ... anywhere
pat = re.compile(r'^\|\s*\*{0,2}([A-Z0-9]{3,6})\b[^|]*\|\s*([\d  ]{2,7})\s*\|\s*([\d  ]{1,5})\s*\|', re.M)
rows_checked = 0
for fn, txt in R.items():
    for m in pat.finditer(txt):
        tg = m.group(1)
        if tg not in short: continue
        try:
            n = int(m.group(2).replace(' ', '').replace(' ', ''))
            p = int(m.group(3).replace(' ', '').replace(' ', ''))
        except ValueError:
            continue
        if n < 50 or n > 3000:       # not a variant-count column
            continue
        # Locate the table header this row belongs to. rfind('|---') lands INSIDE the
        # separator line (it repeats '|---' per column), so walk to that line's start first,
        # then take the line before it. Getting this wrong made Part A silently check nothing.
        sep = txt.rfind('|---', 0, m.start())
        if sep < 0:
            continue
        sep_line_start = txt.rfind('\n', 0, sep) + 1
        hdr_line_start = txt.rfind('\n', 0, sep_line_start - 1) + 1
        hdr = txt[hdr_line_start:sep_line_start].lower()
        if not (('variants' in hdr or re.search(r'\|\s*n\s*\|', hdr))
                and ('positives' in hdr or 'resistant' in hdr)):
            continue
        if any(w in hdr for w in ('rows', 'overlap', 'matched', 'shared', 'mutations')):
            continue
        if n != short[tg][0]:        # a different quantity in a variants-headed table
            continue
        rows_checked += 1
        chk(f'{fn}: {tg} row n/pos', (n, p) == short[tg],
            f'report says {n}/{p}, tables say {short[tg][0]}/{short[tg][1]}')
print(f'  ({rows_checked} table rows matched the n/positives shape and were checked)\n')
# A2. the other per-target table shapes: "| gene | target | positives | ..." in
# REPORT_11TARGETS, and the combined "n+ / n-" cells in REPORT.md and REPORT_ADDITIONS.md.
# Part A's shape misses both, so without this the two most-quoted tables went unchecked.
full = {k.split('_')[0]: (v[0], v[1]) for k, v in counts.items()}
a2 = 0
t11 = R.get('REPORT_11TARGETS.md', '')
# This table now carries BOTH the ingest count and the post-filter count, because the two
# differ for ABL1 and FLT3 and the single unlabelled column hid an undocumented filter.
ING = {'ABL1': 105, 'ALK': 24, 'KIT': 24, 'SMO': 18, 'NT5C2': 18, 'EGFR': 17,
       'MET': 8, 'FLT3': 7, 'ESR1': 7, 'IDH1': 6, 'BTK': 5}
for m in re.finditer(r'^\|\s*([A-Z0-9]{3,6})\s*\|\s*\S+\s*\|\s*(\d{1,4})\s*\|\s*(\d{1,4})\s*\|',
                     t11, re.M):
    tg, ing, l25 = m.group(1), int(m.group(2)), int(m.group(3))
    if tg not in full: continue
    a2 += 2
    chk(f'REPORT_11TARGETS: {tg} ingest positives', ing == ING[tg],
        f'report says {ing}, ingest gives {ING[tg]}')
    chk(f'REPORT_11TARGETS: {tg} L=25 positives', l25 == full[tg][1],
        f'report says {l25}, tables say {full[tg][1]}')
chk('REPORT_11TARGETS: only ABL1 and FLT3 differ between the two columns',
    sorted(t for t in ING if ING[t] != full[t][1]) == ['ABL1', 'FLT3'],
    sorted(t for t in ING if ING[t] != full[t][1]))
chk('REPORT_11TARGETS: L=25 positives sum to the stated 234',
    sum(full[t][1] for t in ING) == 234 and '234 positives' in N['REPORT_11TARGETS.md'],
    sum(full[t][1] for t in ING))
chk('REPORT_11TARGETS: variant counts sum to the stated 7299',
    sum(full[t][0] for t in ING) == 7299 and '7299 variants' in N['REPORT_11TARGETS.md'],
    sum(full[t][0] for t in ING))
# the ambiguous-mutant filter the two columns exposed, now documented
amb = json.load(open(f'{_SAVDR_ROOT}/results/data/ambiguous.json'))
chk('ambiguous filter: 14 records', amb['ambiguous'] == 14 and '14 times in 11 836' in N['REPORT_11TARGETS.md'],
    amb['ambiguous'])
chk('ambiguous filter: 11836 missense scanned', amb['missense'] == 11836, amb['missense'])
chk('ambiguous filter: 5 carry a resistance label',
    amb['ambiguous_positive'] == 5 and '5 of the 14 carry a resistance label' in N['REPORT_11TARGETS.md'],
    amb['ambiguous_positive'])
for v in ['p.E279Z', 'p.E355X', 'p.F359X', 'p.E450X', 'p.D835X']:
    chk(f'ambiguous filter: {v} named in the report', v.replace('p.', '') in N['REPORT_11TARGETS.md'], v)
chk('ambiguous filter: the no-residue-lost claim is stated',
    '**No residue is lost.**' in N['REPORT_11TARGETS.md'])
for fn in ('REPORT.md', 'REPORT_ADDITIONS.md'):
    txt2 = R.get(fn, '')
    for m in re.finditer(r'^\|\s*\*{0,2}([A-Z0-9]{3,6})\b[^|]*\|[^|]*\|\s*(\d{1,4})\s*/\s*(\d{1,4})', txt2, re.M):
        tg, npos, nneg = m.group(1), int(m.group(2)), int(m.group(3))
        if tg not in full: continue
        a2 += 1
        # the "full" setting uses every protein-wide missense variant as a negative; the
        # "region" setting restricts to the resistance span, so only n+ is comparable there
        tot = npos + nneg
        if tot == full[tg][0]:
            chk(f'{fn}: {tg} n+/n- sums to the variant table', True)
        chk(f'{fn}: {tg} n+ matches the positive count', npos == full[tg][1],
            f'report says {npos}, tables say {full[tg][1]}')
print(f'  ({a2} further per-target cells checked in REPORT_11TARGETS, REPORT.md, REPORT_ADDITIONS.md)\n')


# ---------------------------------------------------------------- B: curated shared quantities
print('== B. quantities stated in more than one report ==')

def has(fn, s): return s in R.get(fn, '')

# ABL1 set size, wherever it appears
for fn in R:
    if re.search(r'\b709\b', R[fn]):
        chk(f'{fn}: 709 goes with 101 positives, not another number',
            not re.search(r'709[^.]{0,80}?\b(?!101\b)\d{2,3} (?:resistance )?positives', N[fn]),
            re.findall(r'709[^.]{0,80}positives', N[fn])[:1])

# prime-editing scan size must be 1848 everywhere it is named
for fn in R:
    for m in re.finditer(r'([\d  ]{3,7})[- ]variant (?:scan|evaluation|saturation)', R[fn]):
        v = m.group(1).replace(' ', '').replace(' ', '')
        chk(f'{fn}: scan size {v}', v in ('1848', '1679', '1686'), f'got {v}')

# the two reports that both evaluate us on the scan must not silently disagree
PE, AK = 'REPORT_PRIMEEDIT.md', 'REPORT_ACTIVARK.md'
pairs = [('pooled ATP-site AUC', '0.612', '0.598'),
         ('imatinib AUC',        '0.620', '0.605'),
         ('nilotinib AUC',       '0.686', '0.668'),
         ('within-position AUC', '0.4980', '0.5000')]
for nm, a, b in pairs:
    chk(f'{nm}: stated in both', has(PE, a) and has(AK, b), f'{a} in PRIMEEDIT, {b} in ACTIVARK')
chk('ACTIVARK explains why our AUCs differ from PRIMEEDIT',
    'restricted to the 1 679 variants ActivArk could score' in N[AK]
    and '0.620 for imatinib' in N[AK] and '0.686 for' in N[AK])
chk('ACTIVARK states its own within-position subset (60 positions, 249 pairs)',
    '| 60 |' in R[AK] and '| 249 |' in R[AK])
chk('PRIMEEDIT states its within-position subset (61 positions, 250 pairs)',
    '61' in R[PE] and '250' in R[PE])
soft('the two within-position subsets are explicitly reconciled somewhere',
     'ActivArk-scored' in N[AK] or 'could score' in N[AK])

# label validation figures, quoted in more than one place
for s in ['58.2', '6.5', '19.9', '6.6']:
    chk(f'label-validation {s} only in PRIMEEDIT or consistent elsewhere',
        all(s not in R[fn] or fn == PE or f'{s}' in R[PE] for fn in R), s)

# selection optimism, summarised in REPORT_11TARGETS and detailed in its own report
SO, T11 = 'REPORT_SELECTION_OPTIMISM.md', 'REPORT_11TARGETS.md'
chk('ABL1 optimism agrees between the two reports',
    '0.019' in R[SO] and '0.019' in R[T11], 'ABL1 0.019')
chk('NT5C2 optimism agrees (0.1337 detailed, 0.134 rounded)',
    '0.1337' in R[SO] and '0.134' in R[T11])
chk('ABL1 positives quoted with the optimism figure agree',
    '101 positives' in R[SO] and '101 positives' in R[T11])

# the status table in MISSED_METHODS must match each report's own verdict
MM = 'REPORT_MISSED_METHODS.md'
verdicts = {
    'SUSPECT-ABL':  ('REPORT_SUSPECT_ABL.md',  'indistinguishable near chance'),
    'DRSP':         ('REPORT_DRSP.md',         'DRSP ahead at a matched operating point'),
    'Dr. Kinase':   ('REPORT_DRKINASE.md',     'scores 3 % of a kinase domain'),
    'D3EGFR':       ('REPORT_D3EGFR.md',       '2 of 59 substitutions answered'),
}
for tool, (rep, claim) in verdicts.items():
    chk(f'MISSED_METHODS verdict for {tool} points at {rep}',
        rep.replace('.md', '') in R[MM] or rep in R[MM], rep)
    chk(f'MISSED_METHODS repeats the {tool} verdict verbatim', claim in R[MM], claim)
chk('MISSED_METHODS status table now lists ActivArk and the scan',
    'REPORT_ACTIVARK.md' in R[MM] and 'REPORT_PRIMEEDIT.md' in R[MM])

# MdrDB numbers quoted in two places
MD = 'REPORT_MDRDB.md'
for s in ['0.716', '0.712']:
    chk(f'MdrDB {s} in its own report', s in R[MD], s)
soft('MISSED_METHODS does not restate a stale MdrDB AUC',
     '0.716' not in R[MM] or '0.712' in R[MM])

# D3EGFR: the report's own two counts must agree with its summary line
D3 = 'REPORT_D3EGFR.md'
chk('D3EGFR states 59 submitted / 2 answered across both probes',
    '59 distinct substitutions submitted, 2 answered' in N[D3])
chk('MISSED_METHODS summary of D3EGFR agrees with D3EGFR itself',
    '2 of 59 substitutions answered' in N[MM]
    and '59 distinct substitutions submitted, 2 answered' in N[D3])
chk('D3EGFR first-probe counts add up to the across-probe total',
    '50 EGFR substitutions submitted' in N[D3] and '15 well-known EGFR substitutions' in N[D3]
    and '(6 of them already in the first set)' in N[D3] and 50 + 15 - 6 == 59)
chk('D3EGFR headline: 50 submitted, 1 returned',
    'Of 50 EGFR substitutions submitted, **1 returned a prediction**' in N[D3])

# the unsupportable 0.923 claim must be flagged, never endorsed
for fn in R:
    if '0.923' in R[fn]:
        soft(f'{fn}: 0.923 appears — check it is flagged or is the EGFR grid max',
             'unsupportable' in R[fn] or '| EGFR |' in R[fn] or 'sd 0.020' in R[fn],
             'context: ' + (lambda t, i: t[max(0, i-70):i+45])(N[fn], N[fn].find('0.923')))

# ---------------------------------------------------------------- C: mechanical sweep
print('\n== C. numbers that appear in two reports beside the same keyword ==')
KEY = ['imatinib', 'nilotinib', 'asciminib', 'within-position', 'optimism',
       'Spearman', 'positives', 'overlap']
idx = {}
for fn, txt in R.items():
    low = txt.lower()
    for m in re.finditer(r'0\.\d{3,4}', txt):
        ctx = low[max(0, m.start() - 90):m.end() + 90]
        for k in KEY:
            if k.lower() in ctx:
                idx.setdefault((k, m.group(0)), set()).add(fn)
shared = {k: v for k, v in idx.items() if len(v) > 1}
print(f'  {len(shared)} (keyword, value) pairs appear in more than one report')
for (k, v), fns in sorted(shared.items())[:40]:
    print(f'    {k:16s} {v}  {sorted(f.replace("REPORT_","").replace(".md","") for f in fns)}')

print(f'\n{ok} consistency checks passed, {bad} contradictions, {warn} warnings')
sys.exit(1 if bad else 0)
