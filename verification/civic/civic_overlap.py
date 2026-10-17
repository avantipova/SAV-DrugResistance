"""CIViC as an independent label source for the six analysed targets.

CIViC release 01-Oct-2026, AcceptedClinicalEvidenceSummaries.tsv (CC0).
Only evidence_type=Predictive, significance=Resistance, evidence_direction=Supports is counted
as a CIViC resistance call. Fusion partners (BCR::ABL1, EML4::ALK, ...) are normalised to the
kinase, and compound molecular profiles ("Fusion AND ABL1 T315I") are mined for every
<gene> <X###Y> token, so fusion-context resistance mutations are not lost.
Every CIViC substitution's wild-type residue is checked against the reference FASTA the
pipeline itself reads, which is what establishes that CIViC uses the same numbering.
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

import csv, re, collections, pandas as pd
U = f'{_SAVDR_ROOT}/external/source-data'
F = U + '/civic/01-Oct-2026-AcceptedClinicalEvidenceSummaries.tsv'
TARGET = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'EGFR': 'EGFR',
          'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}
FASTA = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'EGFR': 'EGFR',
         'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}
FUS = re.compile(r'\b(?:BCR|EML4|NPM1|LMNA|TPM3|ATIC|CLTC|STRN|KIF5B)::(ABL1|ALK)\b')
PAT = re.compile(r'\b(ABL1|ALK|EGFR|KIT|NT5C2|SMO)\s+([A-Z])(\d+)([A-Z])\b')

def first_record(p):
    s, started = [], False
    for line in open(p):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)

SEQ = {g: first_record(f'{U}/fasta/{f}.faa') for g, f in FASTA.items()}

sup = collections.defaultdict(lambda: collections.defaultdict(set))    # gene -> variant -> therapies
lvl = collections.defaultdict(lambda: collections.defaultdict(set))
against = collections.defaultdict(set)
for r in csv.DictReader(open(F, newline='', encoding='utf-8'), delimiter='\t'):
    if r['significance'] != 'Resistance' or r['evidence_type'] != 'Predictive': continue
    for g, wt, pos, mu in PAT.findall(FUS.sub(r'\1', r['molecular_profile'])):
        v = f'p.{wt}{pos}{mu}'
        if r['evidence_direction'] == 'Supports':
            sup[g][v].add(r['therapies']); lvl[g][v].add(r['evidence_level'])
        else:
            against[g].add(v)

sc = pd.read_csv(f'{_SAVDR_ROOT}/multipass/scores_all.csv.gz')
out = []
for g, t in TARGET.items():
    d = sc[sc.target == t].drop_duplicates('aa')[['aa', 'y']]
    ours = dict(zip(d.aa, d.y))
    seq = SEQ[g]
    for v in sorted(sup[g], key=lambda s: int(re.search(r'\d+', s).group())):
        wt, pos, mu = v[2], int(v[3:-1]), v[-1]
        out.append(dict(
            gene=g, target=t, variant=v,
            wt_matches_reference=(pos <= len(seq) and seq[pos - 1] == wt),
            civic_levels=','.join(sorted(lvl[g][v])),
            civic_therapies='; '.join(sorted(x for x in sup[g][v] if x)),
            civic_also_contradicted=v in against[g],
            in_our_variant_set=v in ours,
            our_label={1: 'positive', 0: 'negative'}.get(ours.get(v), 'absent')))
df = pd.DataFrame(out)
df.to_csv(f'{_SAVDR_ROOT}/civic/civic_resistance_vs_our_labels.csv', index=False)

print('distinct CIViC resistance SAVs (Supports):', len(df))
print('wild-type residue mismatches against our reference FASTA:',
      int((~df.wt_matches_reference).sum()))
print()
piv = df.pivot_table(index='gene', columns='our_label', values='variant',
                     aggfunc='count', fill_value=0)
for c in ('positive', 'negative', 'absent'):
    if c not in piv: piv[c] = 0
piv['total'] = piv.sum(axis=1)
piv['contested_in_civic'] = df[df.civic_also_contradicted].groupby('gene').variant.count()
print(piv[['total', 'positive', 'negative', 'absent', 'contested_in_civic']].fillna(0).astype(int).to_string())
print()
print('CIViC resistance (Supports) that our labels call NOT resistant:')
for _, r in df[df.our_label == 'negative'].iterrows():
    print(f'  {r.gene:6s} {r.variant:10s} level={r.civic_levels:6s} {r.civic_therapies}')
print()
print('In CIViC but absent from our variant set:')
for _, r in df[df.our_label == 'absent'].iterrows():
    print(f'  {r.gene:6s} {r.variant:10s} level={r.civic_levels:6s} {r.civic_therapies}')
