"""Filter literature-harvested candidate resistance substitutions (2025-06 .. 2026-12,
Europe PMC abstracts) down to ones that can serve as a prospective test set.

Four mechanical gates, in order:
  1. the wild-type residue must match OUR reference FASTA at that position -- this is what
     catches variants co-mentioned from a different gene (BRAF V600E under ALK, KRAS G12C, ...)
  2. the variant must exist in our ingested missense set for that target (otherwise the
     platform has no peptide window for it and cannot score it)
  3. the variant must NOT already be a training positive -- that is what makes it prospective
  4. it must be a single substitution with a specified mutant residue
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

import pandas as pd, re, json, os, glob

CAND = {
 'ABL1':   ['V379I','Y167H','T117P','T315I','Y253H','E255K'],
 'ALK':    ['I1171N','G1202R','V600G','G12C','L858R','L1196M','R110L','F1174L'],
 'BTK':    ['C481S','T474M','T474I','V600E','A428D','T316A','L528W','V416L','C481R','C481F',
            'C481Y','D539A','Y545N','L265P','D539G','T474A'],
 'EGFR':   ['C797S','T790M','L858R','V804F','V843I','S768I','S257L','Q701L','L747S','V802F',
            'C238F','H168P','L718Q','L844V','G12D','I195T','V600E','L816Q','D761N'],
 'ESR1':   ['Y537S','D538G','L370F','E380Q','L536H','L536R','Y537N'],
 'FLT3':   ['F691L','D835Y','E76K','V579A','D835V','Y842C','G13R'],
 'IDH1':   ['R132H','V600E'],
 'KIT':    ['G680R','D842V'],
 'MET':    ['C797S','V600E','I38T','S627N','S257L','W548L','T790M','L718Q','L844V','E340D',
            'D1228N','G12C','D842V'],
 'NT5C2':  ['D431V','D224G'],
 'SMO':    [],
}
TAB = {'ABL1':'ABL1_ENST00000318560','ALK':'ALK','BTK':'BTK_ENST00000308731','EGFR':'EGFR',
       'ESR1':'ESR1_ENST00000206249','FLT3':'FLT3','IDH1':'IDH1_ENST00000415913','KIT':'KIT',
       'MET':'MET','NT5C2':'NT5C2','SMO':'SMO'}
U=f'{_SAVDR_ROOT}/inputs'
def first_record(path):
    s,st=[],False
    for l in open(path):
        if l.startswith('>'):
            if st: break
            st=True; continue
        if st: s.append(l.strip())
    return ''.join(s)

rows=[]
for gene, variants in CAND.items():
    name=TAB[gene]
    seq=first_record(f'{U}/fasta/{name}.faa')
    tab=pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/{name}.csv')
    have=set(tab.aa.str.replace('p.','',regex=False))
    pos_set=set(tab[tab.y==1].aa.str.replace('p.','',regex=False))
    for v in variants:
        m=re.match(r'^([A-Z])(\d+)([A-Z])$',v)
        wt,p,mt=m.group(1),int(m.group(2)),m.group(3)
        ref = seq[p-1] if 1<=p<=len(seq) else None
        gate1 = (ref==wt)
        gate2 = v in have
        gate3 = v not in pos_set
        rows.append(dict(gene=gene, variant=v, pos=p, wt=wt, mt=mt,
            ref_residue=ref, wt_matches_reference=gate1,
            in_our_missense_set=gate2, already_a_training_positive=(not gate3),
            prospective=bool(gate1 and gate2 and gate3)))
D=pd.DataFrame(rows)
D.to_csv(f'{_SAVDR_ROOT}/results/data/prospective_candidates.csv',index=False)
print(f'{len(D)} candidate (gene, variant) pairs harvested\n')
print('gate 1 - wild-type residue matches our reference FASTA:')
bad=D[~D.wt_matches_reference]
print(f'  rejected {len(bad)} of {len(D)}: '
      + ', '.join(f'{r.gene}/{r.variant}(ref {r.ref_residue})' for _,r in bad.iterrows()))
ok1=D[D.wt_matches_reference]
print(f'\ngate 2 - present in our ingested missense set: {int(ok1.in_our_missense_set.sum())} of {len(ok1)}')
miss=ok1[~ok1.in_our_missense_set]
if len(miss): print('  not ingested: '+', '.join(f'{r.gene}/{r.variant}' for _,r in miss.iterrows()))
ok2=ok1[ok1.in_our_missense_set]
print(f'\ngate 3 - already a training positive: {int(ok2.already_a_training_positive.sum())} of {len(ok2)}')
tp=ok2[ok2.already_a_training_positive]
if len(tp): print('  already known to the model: '+', '.join(f'{r.gene}/{r.variant}' for _,r in tp.iterrows()))
P=D[D.prospective]
print(f'\n=== {len(P)} candidates survive all three gates ===')
for g,gg in P.groupby('gene'):
    print(f'  {g:6s} {", ".join(gg.variant)}')
