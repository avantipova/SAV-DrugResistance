"""Rebuild the pipeline's SDF for one target, exactly as aas.ipynb cells 3/12/23 do.

cell 3 : missense rows of csv/, drug columns set to 1 where the AA change is in mutations/,
         peptide = window of <length> centred on the mutation with the MUTANT residue substituted
         (extract_substring); rows with |peptide| != length are dropped.
cell 12: groupby('AA Mutation').head(18)  -- at most one row per length per variant.
cell 23: one SDF per peptide length; class field per molecule.
Class = the drug the variant is resistant to, else 'Not resistant' (multiclass, as in final_table_web).
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

import sys, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
U = f'{_SAVDR_ROOT}/external/source-data'
name, L = sys.argv[1], int(sys.argv[2])
out = sys.argv[3] if len(sys.argv) > 3 else f'{name}_{L}.sdf'

seq = ''.join(l.strip() for l in open(f'{U}/fasta/{name}.faa') if not l.startswith('>'))
g = pd.read_csv(f'{U}/csv/{name}.csv')
g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
m = pd.read_csv(f'{U}/mutations/{name}_mutations.csv')

drug = {}
for aa, d in zip(m['AA Mutation'], m['Drug']):
    drug.setdefault(aa, set()).add(str(d))

def extract(pep, mut, length):           # verbatim port of extract_substring
    parts = mut.split('.')
    pos = int(''.join(c for c in parts[1] if c.isdigit())) - 1
    if pos >= len(pep): return ''
    s = pos - (length // 2); e = s + length
    if s < 0: return ''
    p = list(pep); p[pos] = parts[1][-1]
    f = p[s:e]
    return ''.join(f) if len(f) == length else ''

w = Chem.SDWriter(out); n = skipped = 0
for aa in g['AA Mutation'].dropna().unique():
    pep = extract(seq, aa, L)
    if not pep: continue
    mol = Chem.MolFromSequence(pep)
    if mol is None: skipped += 1; continue
    ds = sorted(drug.get(aa, []))
    mol.SetProp('_Name', aa)
    mol.SetProp('Drug', '\n'.join(ds) if ds else 'Not resistant')
    mol.SetProp('AllDrugs', ';'.join(ds)); mol.SetProp('Peptide', pep)
    w.write(mol); n += 1
w.close()
print(f'{out}: {n} molecules, {skipped} unparsable by RDKit')
