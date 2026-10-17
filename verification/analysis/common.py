"""Shared loaders for SAV-DR re-analysis. Data: ../data/{csv,mutations,fasta}."""
import os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
AA = set('ACDEFGHIKLMNPQRSTVWY')
MISSENSE_RE = re.compile(r'^p\.([A-Z])(\d+)([A-Z])$')

def targets():
    return sorted(f[:-4] for f in os.listdir(os.path.join(DATA, 'csv')) if f.endswith('.csv'))

def read_fasta(name):
    recs, cur = [], None
    with open(os.path.join(DATA, 'fasta', name + '.faa')) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cur = [line[1:], []]; recs.append(cur)
            elif line and cur is not None:
                cur[1].append(line)
    return [(h, ''.join(s)) for h, s in recs]

def parse_missense(hgvs):
    """(wt, pos, mut) for a clean single-residue missense, else None (X, *, =, del, ins...)."""
    m = MISSENSE_RE.match(str(hgvs))
    if not m or m.group(1) not in AA or m.group(3) not in AA or m.group(1) == m.group(3):
        return None
    return m.group(1), int(m.group(2)), m.group(3)

def load_all_variants(name):
    g = pd.read_csv(os.path.join(DATA, 'csv', name + '.csv'))
    return g

def load_resistance(name):
    p = os.path.join(DATA, 'mutations', name + '_mutations.csv')
    return pd.read_csv(p) if os.path.exists(p) else None
