"""Recompute the published selection metric honestly, from the project's own spreadsheet.

No MultiPASS runs: total_table_avg_iap holds the per-drug Average IAP for every
(peptide length x MNA level x charge) configuration already, 12262 rows.

The sheet's `avg_iap` has two defects, both verified against its own cells:

  1. it averages ONLY the non-zero drug classes, so a configuration that made MultiPASS decline
     a hard drug class gets a better score for failing. On ABL1_ENST00000318560 the top-ranked
     L23/MNA7 scores 0.8706 with Bosutinib excluded and 0.7255 with it counted as 0 -- which
     moves it from 1st of 5 to 4th of 5.
  2. it averages `Not resistant` in as if it were a drug, adding the large easy negative class
     to a mean that is presented as per-drug performance.

Three metrics are computed for every configuration:
    sheet_rule   mean over non-zero drug columns + Not resistant      (what was published)
    zeros_count  mean over ALL drug columns + Not resistant, declines scored 0
    drugs_only   mean over ALL drug columns, declines scored 0, Not resistant EXCLUDED
and the per-target best configuration under each, so the cost of the two defects is explicit.
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

import zipfile, sys, numpy as np, pandas as pd
import xml.etree.ElementTree as ET

F = sys.argv[1] if len(sys.argv) > 1 else \
    '/root/.claude/uploads/2c3b11c7-e846-5338-83a9-42b57da727e2/923ce595-total_table_avg_iap_final_2.xlsx'
NS = '{http://purl.oclc.org/ooxml/spreadsheetml/main}'
NSR = '{http://purl.oclc.org/ooxml/officeDocument/relationships}'
z = zipfile.ZipFile(F)
sst = [''.join(t.text or '' for t in si.iter(NS + 't'))
       for si in ET.fromstring(z.read('xl/sharedStrings.xml'))]
rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
wbx = ET.fromstring(z.read('xl/workbook.xml'))
sheets = {s.get('name'): rels[s.get(NSR + 'id')] for s in wbx.iter(NS + 'sheet')}


def read(name):
    root = ET.fromstring(z.read('xl/' + sheets[name].lstrip('/').replace('xl/', '')))
    out = []
    for row in root.iter(NS + 'row'):
        v = []
        for c in row.iter(NS + 'c'):
            t, vv, isx = c.get('t'), c.find(NS + 'v'), c.find(NS + 'is')
            if t == 's' and vv is not None: v.append(sst[int(vv.text)])
            elif isx is not None: v.append(''.join(x.text or '' for x in isx.iter(NS + 't')))
            elif vv is not None: v.append(vv.text)
            else: v.append('')
        out.append(v)
    return out


rows = read('total_table_avg_iap')
hdr = [h.strip() for h in rows[0]]
df = pd.DataFrame(rows[1:], columns=hdr + [''] * (max(0, max(len(r) for r in rows[1:]) - len(hdr))))
NONDRUG = {'protein', 'peptide_length', 'mna_length', 'charge', 'avg_iap', 'for_web',
           'Выбор моделей', ''}
drugcols = [c for c in hdr if c not in NONDRUG and c != 'Not resistant']
print(f'rows: {len(df)}   drug columns: {len(drugcols)}')


def num(s):
    try: return float(s)
    except Exception: return np.nan


for c in drugcols + ['Not resistant', 'peptide_length', 'mna_length', 'charge', 'avg_iap']:
    if c in df.columns: df[c] = df[c].map(num)
df = df[df.protein.astype(str).str.len() > 3].copy()
df = df.dropna(subset=['peptide_length', 'mna_length'])

# The sheet stores an inapplicable drug as 0 -- the same value it uses for "MultiPASS declined
# to predict this class" -- so a 0 is ambiguous. The label file is not the right arbiter either:
# ABL1's labels name five drugs but that run only modelled two, so three are 0 everywhere.
#
# Data-driven rule: a class was ATTEMPTED for a protein if it has a non-zero IAP in at least one
# of that protein's configurations. Within a configuration, a 0 on an attempted class is then a
# genuine decline; a class never attempted is ignored entirely.
#
# One further guard. "Non-zero somewhere" alone is not enough: row 2835 of this sheet is a
# paste-misaligned ABL1_ENST00000318560 row (L15/MNA3) carrying three ALK-inhibitor values
# (Crizotinib 0.7746, Alectinib 0.9281, Lorlatinib 0.8868) for a protein that has no
# ALK-inhibitor labels, and the sheet averages all nine numbers into its avg_iap
# (0.858244... = 7.7242/9), which is why that one row's avg_iap does not reproduce from its own
# ABL1 columns. So the attempted set is intersected with the drugs the target's label file
# actually names, which makes it immune to single-row contamination.
import glob, os
labelled = {}
for p_ in glob.glob(f'{_SAVDR_ROOT}/inputs/mutations/*_mutations.csv'):
    stem = os.path.basename(p_)[:-len('_mutations.csv')]
    try:
        dd = pd.read_csv(p_)
    except Exception:
        continue
    labelled[stem.replace('_', '').upper()] = {str(x).strip().lower() for x in dd['Drug'].dropna()}

attempted = {}
for p, g in df.groupby('protein'):
    lab = labelled.get(str(p).replace('_', '').upper(), set())
    cols = []
    for c in drugcols:
        col = pd.to_numeric(g[c], errors='coerce')
        if (col.fillna(0) > 0).any() and c.strip().lower() in lab: cols.append(c)
    attempted[p] = cols
print('attempted drug classes per target: '
      + ', '.join(f'{p}={len(c)}' for p, c in sorted(attempted.items())))
missing_key = set()

res = []
for _, r in df.iterrows():
    cols = attempted.get(r.protein, [])
    if not cols:
        missing_key.add(r.protein); continue
    v = np.array([r[c] for c in cols], dtype=float)
    v = v[~np.isnan(v)]
    nr = r.get('Not resistant', np.nan)
    nonzero = v[v > 0]
    withnr = np.r_[v, [nr]] if not np.isnan(nr) else v
    nonzero_withnr = np.r_[nonzero, [nr]] if not np.isnan(nr) else nonzero
    res.append(dict(
        protein=r.protein, L=int(r.peptide_length), MNA=int(r.mna_length),
        charge=r.get('charge', np.nan), sheet_avg=r.get('avg_iap', np.nan),
        n_drug_classes=len(v), n_declined=int((v == 0).sum()),
        sheet_rule=float(nonzero_withnr.mean()) if len(nonzero_withnr) else np.nan,
        zeros_count=float(withnr.mean()) if len(withnr) else np.nan,
        drugs_only=float(v.mean()) if len(v) else np.nan))
R = pd.DataFrame(res)
if missing_key:
    print(f'no label file matched, skipped: {sorted(missing_key)}')



ok = R.dropna(subset=['sheet_avg', 'sheet_rule'])
d = (ok.sheet_rule - ok.sheet_avg).abs()
print(f'reproduces the sheet\'s own avg_iap: {(d < 5e-4).sum()}/{len(ok)} rows within 5e-4, '
      f'median |diff| {d.median():.6f}')

out = []
for p, g in R.groupby('protein'):
    row = dict(protein=p, n_configs=len(g), n_with_declines=int((g.n_declined > 0).sum()),
               max_drug_classes=int(g.n_drug_classes.max()))
    for metric in ('sheet_rule', 'zeros_count', 'drugs_only'):
        gg = g.dropna(subset=[metric])
        if gg.empty: continue
        b = gg.loc[gg[metric].idxmax()]
        row[f'{metric}_best'] = round(float(b[metric]), 4)
        row[f'{metric}_cfg'] = f'L{int(b.L)}/MNA{int(b.MNA)}'
    out.append(row)
O = pd.DataFrame(out).sort_values('protein')
O['drop_zeros_vs_sheet'] = (O.zeros_count_best - O.sheet_rule_best).round(4)
O['drop_drugsonly_vs_sheet'] = (O.drugs_only_best - O.sheet_rule_best).round(4)
O.to_csv(f'{_SAVDR_ROOT}/multipass/table1_recomputed.csv', index=False)
R.to_csv(f'{_SAVDR_ROOT}/multipass/table1_all_configs.csv', index=False)

pd.set_option('display.width', 250)
print()
print('=' * 120)
print('BEST CONFIGURATION PER TARGET UNDER EACH METRIC')
print('=' * 120)
print(O[['protein', 'n_configs', 'max_drug_classes', 'n_with_declines',
         'sheet_rule_best', 'sheet_rule_cfg', 'zeros_count_best', 'zeros_count_cfg',
         'drugs_only_best', 'drugs_only_cfg',
         'drop_zeros_vs_sheet', 'drop_drugsonly_vs_sheet']].to_string(index=False))
print()
print(f'configurations whose score is inflated by at least one declined class: '
      f'{int((R.n_declined > 0).sum())} of {len(R)} '
      f'({(R.n_declined > 0).mean():.1%})')
multi = O[O.max_drug_classes > 1]
print(f'\non multi-drug targets ({len(multi)}): counting declines as 0 costs '
      f'{-multi.drop_zeros_vs_sheet.mean():.4f} AUC on average '
      f'(worst {-multi.drop_zeros_vs_sheet.min():.4f}); also dropping "Not resistant" costs '
      f'{-multi.drop_drugsonly_vs_sheet.mean():.4f} on average')
print(f'configuration changes under the honest metric: '
      f'{int((O.sheet_rule_cfg != O.drugs_only_cfg).sum())} of {len(O)} targets')
