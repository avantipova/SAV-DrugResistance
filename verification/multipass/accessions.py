"""Audit Table S7 against the FASTA files the pipeline actually used.

The FASTA headers carry the RefSeq protein accession and the sequence gives its length,
so the identity of each target's sequence can be settled offline, without UniProt.
Reference lengths for the UniProt entries listed in Table S7 are taken from the FASTA
files themselves where the same accession is claimed twice (the BTK case), and otherwise
flagged for a one-line manual check."""
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

import os, re, pandas as pd
U = f'{_SAVDR_ROOT}/external/source-data'
S7 = {  # Table S7 as published: target -> (Ensembl transcript, accession given)
 'ABL1': ('ENST00000372348.9', 'P00519'),
 'ABL1_ENST00000318560': ('ENST00000318560.5', 'NM_005157.4'),
 'ALK': ('ENST00000389048.8', 'Q9UM73'),
 'ALK_ENST00000618119': ('ENST00000618119.4', 'A0A087WZL3'),
 'BTK': ('ENST00000621635.4', 'Q06187'),
 'BTK_ENST00000308731': ('ENST00000308731.8', 'Q06187.3'),
 'EGFR': ('ENST00000275493.7', 'P00533'),
 'EGFR_ENST00000454757': ('ENST00000454757.6', 'NM_005228.3'),
 'ESR1': ('ENST00000440973.5', 'P03372'),
 'ESR1_ENST00000206249': ('ENST00000206249.8', 'P03372'),
 'FLT3': ('ENST00000241453.12', 'P36888'),
 'IDH1_ENST00000415913': ('ENST00000415913.5', 'O75874'),
 'KIT': ('ENST00000288135.6', 'P10721'),
 'MET': ('ENST00000318493.11', 'P08581'),
 'MET_ENST00000397752': ('ENST00000397752.8', 'NM_000245.4'),
 'NT5C2': ('ENST00000343289.9', 'P49902'),
 'NT5C2_ENST00000404739': ('ENST00000404739.8', 'P49902'),
 'SMO': ('ENST00000249373.8', 'Q99835'),
}
rows = []
for t, (ens, acc) in S7.items():
    p = f'{U}/fasta/{t}.faa'
    if not os.path.exists(p):
        rows.append(dict(target=t, note='fasta not staged')); continue
    hdr, seq, recs = None, [], 0
    for line in open(p):
        if line.startswith('>'):
            recs += 1
            if hdr is None: hdr = line[1:].strip()
        elif recs == 1: seq.append(line.strip())
    seq = ''.join(seq)
    refseq = hdr.split()[0] if hdr else ''
    rows.append(dict(target=t, transcript=ens, s7_accession=acc,
                     kind='RefSeq' if acc.startswith('N') else 'UniProt',
                     fasta_refseq=refseq, fasta_len=len(seq), fasta_records=recs,
                     isoform=(re.search(r'\[isoform=([^\]]+)\]', hdr).group(1) if hdr and 'isoform=' in hdr else '')))
d = pd.DataFrame(rows)
# accessions claimed for more than one target, with different sequence lengths -> impossible
d['acc_base'] = d.s7_accession.str.split('.').str[0]
dup = d.groupby('acc_base').filter(lambda g: len(g) > 1)
d['conflict'] = ''
for acc, g in dup.groupby('acc_base'):
    if g.fasta_len.nunique() > 1:
        d.loc[g.index, 'conflict'] = f'{acc} claimed for {g.fasta_len.nunique()} different lengths: ' + \
            ', '.join(f'{r.target}={r.fasta_len}aa' for r in g.itertuples())
d.drop(columns=['acc_base']).to_csv('09_accession_audit.csv', index=False)
pd.set_option('display.width', 220); pd.set_option('display.max_colwidth', 60)
print(d[['target','s7_accession','kind','fasta_refseq','fasta_len','isoform','conflict']].to_string(index=False))
