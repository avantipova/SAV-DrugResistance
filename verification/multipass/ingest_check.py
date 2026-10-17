"""Validate every target BEFORE any expensive run.

A TARGET is a stem that has all three files: fasta/<stem>.faa, csv/<stem>.csv and
mutations/<stem>_mutations.csv. `ABL1` and `ABL1_ENST00000318560` are two DIFFERENT targets
with different numbering, so they are never mixed -- an earlier version of this script grouped
by gene name and silently paired one transcript's labels against another's missense set, which
showed up as ABL1 having 0 positives instead of 105.

Checks per target:
  * multi-record FASTA (KIT has 8 records, SMO 3, FLT3 11 -- only the first is canonical)
  * numbering: every p.X###Y wild-type residue against residue ### of the sequence. Where a gene
    has several candidate sequences, the one with the best agreement is reported, so the right
    isoform is chosen from the data rather than from external metadata.
  * label hygiene: labels that are not single substitutions, or absent from the missense set
  * feasibility: positives, and positions carrying both a positive and a negative (the limiting
    quantity for the within-position test)

Usage: [DATA=<dir>] python3 ingest_check.py [stem ...]
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

import glob, os, re, sys, collections, pandas as pd

U = os.environ.get('DATA', f'{_SAVDR_ROOT}/external/source-data')
SAV = re.compile(r'^p\.([A-Z])(\d+)([A-Z])$')
ONLY = set(sys.argv[1:])
MINPOS = 15
MINBOTH = 15


def records(path):
    out, name, s = [], None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None: out.append((name, ''.join(s)))
            name, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if name is not None: out.append((name, ''.join(s)))
    return out


# every stem that has labels + a missense table; the FASTA may be under a sibling stem
stems = sorted(os.path.basename(p)[:-len('_mutations.csv')]
               for p in glob.glob(f'{U}/mutations/*_mutations.csv'))
fastas = {os.path.basename(p)[:-4]: p for p in glob.glob(f'{U}/fasta/*.faa')}

rows = []
for stem in stems:
    if ONLY and stem not in ONLY: continue
    csvp = f'{U}/csv/{stem}.csv'
    labp = f'{U}/mutations/{stem}_mutations.csv'
    gene = stem.split('_')[0]
    # candidate sequences: the exact stem first, then any sibling of the same gene
    cand_stems = ([stem] if stem in fastas else []) + \
                 sorted(k for k in fastas if k != stem and k.split('_')[0] == gene)
    print('=' * 100)
    print(f'{stem}')
    print('=' * 100)
    if not os.path.exists(csvp):
        print('  NO missense csv -- cannot validate')
        rows.append(dict(target=stem, status='no_csv')); print(); continue
    if not cand_stems:
        print(f'  NO FASTA for gene {gene} -- cannot validate numbering or build peptides')
        g = pd.read_csv(csvp)
        g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
        lab = set(pd.read_csv(labp)['AA Mutation'].dropna().unique())
        usable = len([a for a in lab if SAV.match(str(a)) and a in set(g['AA Mutation'])])
        print(f'  (labels would give {usable} usable positives out of {g["AA Mutation"].nunique()} missense)')
        rows.append(dict(target=stem, status='no_fasta', n_missense=int(g['AA Mutation'].nunique()),
                         n_positives=usable))
        print(); continue

    g = pd.read_csv(csvp)
    g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
    missense = set(g['AA Mutation'].dropna().unique())
    labels = set(pd.read_csv(labp)['AA Mutation'].dropna().unique())
    d = pd.read_csv(labp)

    best, bestrate, multi = None, -1, 0
    for cs in cand_stems:
        r = records(fastas[cs])
        seq = r[0][1]
        if cs == cand_stems[0]: multi = int(len(r) > 1)
        okc = tot = 0
        for aa in missense:
            m = SAV.match(str(aa))
            if not m: continue
            tot += 1
            if int(m.group(2)) <= len(seq) and seq[int(m.group(2)) - 1] == m.group(1): okc += 1
        rate = okc / tot if tot else 0
        tagm = f'  records={len(r)}' + ('  <-- MULTI-RECORD, first record only' if len(r) > 1 else '')
        print(f'  {cs:26s} {len(seq):5d} aa  wild-type matches {okc:5d}/{tot:5d} = {rate:7.2%}{tagm}')
        if rate > bestrate: best, bestrate = cs, rate
    verdict = '** PERFECT **' if bestrate == 1 else '** IMPERFECT -- investigate before running **'
    print(f'  -> USE {best}.faa   {verdict}')

    not_sav = sorted(a for a in labels if not SAV.match(str(a)))
    not_in = sorted(a for a in labels if SAV.match(str(a)) and a not in missense)
    usable = [a for a in labels if SAV.match(str(a)) and a in missense]
    drugs = sorted({str(x).strip().lower() for x in d['Drug'].dropna()})
    posn = collections.defaultdict(lambda: [0, 0])
    for aa in missense:
        m = SAV.match(str(aa))
        if m: posn[int(m.group(2))][0 if aa in labels else 1] += 1
    both = sum(1 for v in posn.values() if v[0] and v[1])

    print(f'  labels {len(labels)}: {len(usable)} usable, {len(not_sav)} not single substitutions, '
          f'{len(not_in)} absent from the missense set')
    if not_sav: print(f'    not single: {", ".join(map(str, not_sav[:10]))}' + (' ...' if len(not_sav) > 10 else ''))
    if not_in: print(f'    absent    : {", ".join(map(str, not_in[:10]))}' + (' ...' if len(not_in) > 10 else ''))
    print(f'  drugs ({len(drugs)}): {", ".join(drugs)}')
    print(f'  {len(usable)} pos / {len(missense) - len(usable)} neg; {both} positions carry both')
    print(f'  per-target row: {"OK" if len(usable) >= MINPOS else f"UNDERPOWERED (<{MINPOS} positives)"}'
          f'   |   within-position test: {"OK" if both >= MINBOTH else "too few positions"}'
          f'   |   LOTO training data: always usable')
    print()
    rows.append(dict(target=stem, status='ok' if bestrate == 1 else 'numbering_mismatch',
                     use_fasta=best + '.faa', wt_agreement=round(bestrate, 5),
                     multi_record=multi, n_missense=len(missense), n_positives=len(usable),
                     n_dropped=len(not_sav) + len(not_in), n_drugs=len(drugs),
                     drugs=';'.join(drugs), positions_with_both=both))

df = pd.DataFrame(rows)
df.to_csv(f'{_SAVDR_ROOT}/multipass/ingest_check.csv', index=False)
cols = [c for c in ['target', 'status', 'use_fasta', 'wt_agreement', 'multi_record', 'n_missense',
                    'n_positives', 'n_dropped', 'n_drugs', 'positions_with_both'] if c in df.columns]
print('=' * 100)
print(df[cols].to_string(index=False))
okdf = df[df.status == 'ok'] if 'status' in df else df
print(f'\nusable targets: {len(okdf)}   '
      f'with >={MINPOS} positives: {int((okdf.n_positives >= MINPOS).sum()) if len(okdf) else 0}   '
      f'no FASTA: {int((df.status == "no_fasta").sum())}   '
      f'numbering problems: {int((df.status == "numbering_mismatch").sum())}')
if len(okdf):
    print(f'total positions carrying both a positive and a negative: '
          f'{int(okdf.positions_with_both.sum())}  (was 34 with six targets)')
