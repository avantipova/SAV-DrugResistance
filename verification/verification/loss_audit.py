"""Why n_pos is smaller than the number of resistance mutations COSMIC lists, per target."""
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

import pandas as pd, re, csv
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
U = f'{_SAVDR_ROOT}/external/source-data'
CFG = {'ABL1_ENST00000318560': 17, 'ALK': 33, 'EGFR': 37, 'KIT': 23, 'NT5C2': 25, 'SMO': 31}

def first_record(path):
    s, started = [], False
    for line in open(path):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)

def extract(pep, mut, length):
    parts = mut.split('.'); pos = int(''.join(c for c in parts[1] if c.isdigit())) - 1
    if pos >= len(pep): return ''
    s = pos - length // 2; e = s + length
    if s < 0: return ''
    p = list(pep); p[pos] = parts[1][-1]; f = p[s:e]
    return ''.join(f) if len(f) == length else ''

rows, detail = [], []
for t, L in CFG.items():
    seq = first_record(f'{U}/fasta/{t}.faa')
    g = pd.read_csv(f'{U}/csv/{t}.csv')
    g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
    cs = set(g['AA Mutation'].dropna().unique())
    res = set(pd.read_csv(f'{U}/mutations/{t}_mutations.csv')['AA Mutation'].dropna())
    c = dict(kept=0, not_single_substitution=0, not_in_cosmic_missense=0,
             window_off_protein_end=0, unparseable_residue_code=0)
    for aa in sorted(res):
        if not re.match(r'^p\.[A-Z]\d+[A-Z]$', str(aa)):
            c['not_single_substitution'] += 1; detail.append((t, aa, 'not a single substitution')); continue
        if aa not in cs:
            c['not_in_cosmic_missense'] += 1; detail.append((t, aa, 'absent from the COSMIC missense export')); continue
        pep = extract(seq, aa, L)
        if not pep:
            c['window_off_protein_end'] += 1; detail.append((t, aa, f'{L}-residue window runs off the protein end')); continue
        if Chem.MolFromSequence(pep) is None:
            c['unparseable_residue_code'] += 1; detail.append((t, aa, 'ambiguous/stop residue code, no peptide')); continue
        c['kept'] += 1
    rows.append(dict(target=t, resistance_mutations_listed=len(res), **c))
df = pd.DataFrame(rows)
df.to_csv(f'{_SAVDR_ROOT}/verification/positive_loss_audit.csv', index=False)
with open(f'{_SAVDR_ROOT}/verification/positives_excluded.csv', 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['target', 'aa_mutation', 'reason']); w.writerows(detail)
print(df.to_string(index=False))
print(f'\ntotal listed {df.resistance_mutations_listed.sum()}, kept as positives {df.kept.sum()}, '
      f'excluded {df.resistance_mutations_listed.sum() - df.kept.sum()}')
print('\nexcluded, in full:')
for t, aa, r in detail: print(f'   {t:24s} {aa:20s} {r}')
