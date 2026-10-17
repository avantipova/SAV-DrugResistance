"""Build variant tables for the 11-gene set, one transcript per gene.

Differences from the original 6-target pipeline:

  * one transcript per gene, chosen by ingest_check.py as the sequence whose residues agree with
    the export's own p.X###Y wild-type letters (ABL1 on the wrong isoform agrees 7.23%, on the
    right one 100.00%, so this is decisive).
  * a WILD-TYPE GUARD: any variant whose wild-type letter disagrees with the reference is
    dropped. On the original six this changes nothing (all 100%); on FLT3 it drops p.V571I,
    whose reference residue is R at 571 and which is not a resistance label.
  * a single common window length for every gene. The original per-target lengths came from the
    hand-curated `best_comb` shortlist, which the new genes have no entry in; and for a
    cross-protein base a homogeneous window is the defensible choice anyway.

Usage: DATA=<dir> python3 build11.py <L>
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

import os, re, sys, json, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

U = os.environ.get('DATA', f'{_SAVDR_ROOT}/inputs')
L = int(sys.argv[1]) if len(sys.argv) > 1 else 25
OUT = f'/home/claude/tab11_L{L}'
os.makedirs(OUT, exist_ok=True)
SAV = re.compile(r'^p\.([A-Z])(\d+)([A-Z])$')

# stem -> fasta stem, as validated by ingest_check.py (one transcript per gene)
TARGETS = {
    'ABL1_ENST00000318560': 'ABL1_ENST00000318560',
    'ALK': 'ALK',
    'BTK_ENST00000308731': 'BTK_ENST00000308731',
    'EGFR': 'EGFR',
    'ESR1_ENST00000206249': 'ESR1_ENST00000206249',
    'FLT3': 'FLT3',
    'IDH1_ENST00000415913': 'IDH1_ENST00000415913',
    'KIT': 'KIT',
    'MET': 'MET',
    'NT5C2': 'NT5C2',
    'SMO': 'SMO',
}


def first_record(path):
    """Only the first record: KIT has 8, SMO 3, FLT3 11, and the rest are other isoforms."""
    s, started = [], False
    for line in open(path):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)


def window(seq, mut, length):
    m = SAV.match(mut)
    pos = int(m.group(2)) - 1
    if pos >= len(seq): return ''
    s = pos - length // 2
    e = s + length
    if s < 0: return ''
    p = list(seq); p[pos] = m.group(3)
    f = p[s:e]
    return ''.join(f) if len(f) == length else ''


summary = []
for stem, fst in TARGETS.items():
    seq = first_record(f'{U}/fasta/{fst}.faa')
    g = pd.read_csv(f'{U}/csv/{stem}.csv')
    g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
    res = set(pd.read_csv(f'{U}/mutations/{stem}_mutations.csv')['AA Mutation'].dropna())
    drugmap = {}
    md = pd.read_csv(f'{U}/mutations/{stem}_mutations.csv')
    for _, r in md.iterrows():
        if pd.notna(r['AA Mutation']) and pd.notna(r['Drug']):
            drugmap.setdefault(r['AA Mutation'], set()).add(str(r['Drug']).strip())

    rows = []
    n_notsav = n_wtbad = n_nowin = n_nomol = 0
    for aa in g['AA Mutation'].dropna().unique():
        m = SAV.match(str(aa))
        if not m: n_notsav += 1; continue
        wt, pos = m.group(1), int(m.group(2))
        if pos > len(seq) or seq[pos - 1] != wt:        # the wild-type guard
            n_wtbad += 1; continue
        pep = window(seq, aa, L)
        if not pep: n_nowin += 1; continue
        if Chem.MolFromSequence(pep) is None: n_nomol += 1; continue
        rows.append(dict(aa=aa, pos=pos, pep=pep, y=int(aa in res),
                         drugs=';'.join(sorted(drugmap.get(aa, ()))) or ''))
    df = pd.DataFrame(rows)
    df.to_csv(f'{OUT}/{stem}.csv', index=False)
    summary.append(dict(target=stem, fasta=fst, seq_len=len(seq), L=L, n=len(df),
                        n_pos=int(df.y.sum()), n_neg=int((1 - df.y).sum()),
                        dropped_not_sav=n_notsav, dropped_wt_mismatch=n_wtbad,
                        dropped_window=n_nowin, dropped_rdkit=n_nomol,
                        n_positions=int(df.pos.nunique())))
    print(f'{stem:24s} seq={len(seq):5d}  {len(df):5d} variants  '
          f'{int(df.y.sum()):4d}+/{int((1-df.y).sum()):5d}-  '
          f'dropped: wt={n_wtbad} window={n_nowin} rdkit={n_nomol} notsav={n_notsav}', flush=True)

s = pd.DataFrame(summary)
s.to_csv(f'{OUT}/_summary.csv', index=False)
print()
print(f'TOTAL: {s.n.sum()} variants, {s.n_pos.sum()} positives across {len(s)} targets')
print(f'wild-type-guard drops: {s.dropped_wt_mismatch.sum()} '
      f'({", ".join(f"{r.target}={r.dropped_wt_mismatch}" for _, r in s.iterrows() if r.dropped_wt_mismatch)})')
json.dump(list(TARGETS), open(f'{OUT}/_targets.json', 'w'))
