"""Table S7, rebuilt and verified by sequence checksum rather than by length.

For every FASTA the pipeline actually reads, the CRC64 of the first record is computed and
matched against the CRC64 UniProt reports for the cited accession. Where it does not match,
the sequence is identified through UniParc (query checksum:<CRC64>), which returns the
UniProtKB isoform that holds exactly that sequence. RefSeq protein<->mRNA pairings are taken
from the UniProt entry's own cross-references, so no NCBI lookup is involved.
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

import glob, os, csv
from Bio.SeqUtils.CheckSum import crc64
U = f'{_SAVDR_ROOT}/inputs/fasta'

def first_record(p):
    hdr, seq, started = None, [], False
    for line in open(p):
        if line.startswith('>'):
            if started: break
            hdr = line[1:].strip(); started = True; continue
        if started: seq.append(line.strip())
    return hdr, ''.join(seq)

# CRC64 as reported by UniProt for each cited accession (rest.uniprot.org, fields=sequence)
UNIPROT_CRC = {
    'P00519': ('85FE6C1C0E483EA2', 1130), 'Q06187': ('DF06B5D1FEC257CC', 659),
    'P00533': ('D8A2A50B4EFB6ED2', 1210), 'Q9UM73': ('0733D6C4FD212F41', 1620),
    'P08581': ('9CF896D1273905C3', 1390), 'P03372': ('5455C57AB0CCCAA7', 595),
    'P49902': ('4C27D762575E0EA2', 561),  'P10721': ('81B0CD76817F3454', 976),
    'Q99835': ('8B4C459B34D13F83', 787),  'P36888': ('6C1995718F352ECE', 993),
    'O75874': ('60428B0B6E5851DC', 414),
    'A0A087WZL3': ('7A43DA6A78687E85', 1243),
}
# sequences whose CRC64 matched no cited canonical: resolved through UniParc
UNIPARC = {
    'DF9D4512F78FFE52': ('UPI000013E4DE', 'P00519-2',  'NP_009297.2',    'NM_007313.3'),
    '8C582E0CD4D6280F': ('UPI000021216F', 'Q06187-2',  'NP_001274273.1', 'NM_001287344.2'),
    'C5F64DCBBC13CC88': ('UPI000014033E', 'P08581-2',  'NP_001120972.1', 'NM_001127500.3'),
    '9CEE26B795A041B2': ('UPI000387BFB2', 'E9PFD7 (INACTIVE: deleted from Ensembl)', '', ''),
}
# RefSeq protein <-> mRNA pairings quoted from the UniProt entries' own xrefs
REFSEQ = {
    'NP_005148.2': ('NM_005157.6', 'P00519-1'), 'NP_009297.2': ('NM_007313.3', 'P00519-2'),
    'NP_000052.1': ('NM_000061.3', 'Q06187-1'), 'NP_001274273.1': ('NM_001287344.2', 'Q06187-2'),
    'NP_005219.2': ('NM_005228.5', 'P00533-1'),
    'NP_000236.2': ('NM_000245.4', 'P08581-1'), 'NP_001120972.1': ('NM_001127500.3', 'P08581-2'),
}
# what Table S7 currently prints, per row
S7 = {
 'ABL1': 'P00519', 'ABL1_ENST00000318560': 'NM_005157.4', 'ALK': 'Q9UM73',
 'ALK_ENST00000618119': 'A0A087WZL3', 'BTK': 'Q06187', 'BTK_ENST00000308731': 'Q06187.3',
 'EGFR': 'P00533', 'EGFR_ENST00000454757': 'NM_005228.3', 'ESR1': 'P03372',
 'ESR1_ENST00000206249': 'P03372', 'FLT3': 'P36888', 'IDH1_ENST00000415913': 'O75874',
 'KIT': 'P10721', 'MET': 'P08581', 'MET_ENST00000397752': 'NM_000245.4',
 'NT5C2': 'P49902', 'NT5C2_ENST00000404739': 'P49902', 'SMO': 'Q99835',
}
ANALYSED = {'ABL1_ENST00000318560', 'ALK', 'EGFR', 'KIT', 'NT5C2', 'SMO'}

rows, by_crc = [], {}
for f in sorted(glob.glob(U + '/*.faa')):
    t = os.path.basename(f)[:-4]
    hdr, seq = first_record(f)
    c = crc64(seq).replace('CRC-', '')
    by_crc.setdefault(c, []).append(t)
    s7 = S7.get(t, '')
    base = s7.split('.')[0] if s7.startswith(('P', 'Q', 'O', 'A')) else ''
    if base in UNIPROT_CRC and UNIPROT_CRC[base][0] == c:
        correct, how = base, 'CRC64 == UniProt canonical'
    else:
        hit = next((a for a, (cc, _) in UNIPROT_CRC.items() if cc == c), None)
        if hit:
            correct, how = hit, 'CRC64 == UniProt canonical (different accession than S7 prints)'
        elif c in UNIPARC:
            correct, how = UNIPARC[c][1], f'CRC64 -> UniParc {UNIPARC[c][0]}'
        else:
            correct, how = '', 'unresolved'
    np_id = hdr.split()[0] if hdr and hdr.startswith('NP_') else ''
    nm_id, iso = REFSEQ.get(np_id, ('', ''))
    verdict = 'ok' if s7 == correct else 'FIX'
    rows.append(dict(row=t, analysed=t in ANALYSED, length=len(seq), crc64=c,
                     fasta_header=hdr[:60], s7_prints=s7, correct_uniprot=correct,
                     refseq_protein=np_id, refseq_mrna=nm_id, uniprot_isoform=iso,
                     verified_by=how, verdict=verdict))

with open(f'{_SAVDR_ROOT}/multipass/Table_S7_corrected.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

for r in rows:
    flag = '  ' if r['verdict'] == 'ok' else '->'
    star = '*' if r['analysed'] else ' '
    print(f"{flag}{star}{r['row']:24s} {r['length']:5d}  S7={r['s7_prints']:12s} "
          f"correct={r['correct_uniprot']:12s} {r['verified_by']}")
print('\nidentical sequences sharing a CRC64 (double-counted rows):')
for c, ts in by_crc.items():
    if len(ts) > 1: print('   ', c, '->', ', '.join(ts))
