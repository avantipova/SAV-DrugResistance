"""Score the prospective resistance substitutions with the deployed model configuration.

For each target the manuscript's Table 1 names the configuration the platform serves. Those
configurations are reproduced here (peptide length, MNA level, uncharged), the base is trained on
ALL of that target's variants for the matching drug -- as the deployed model is -- and the
prospective variants are scored. No cross-validation is needed: these variants were reported
after the COSMIC v100 freeze and are absent from the training labels, so they are held out by time.

Both Pa (resistant) and Pi (not resistant) are read, because the platform's decision rule is
Pa > Pi.
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

import os, sys, glob, json, shutil, subprocess, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
# target -> (table name, Table 1 peptide length, Table 1 MNA level, drug string to match)
CFG = {
 'BTK':  ('BTK_ENST00000308731',  9,  3, 'Ibrutinib'),
 'EGFR': ('EGFR',                37, 30, 'Osimertinib'),
 'ESR1': ('ESR1_ENST00000206249',11, 12, 'Endocrine therapy'),
 'FLT3': ('FLT3',                 7,  5, 'Quizartinib'),
}
TEST = {'BTK': ['T474I','L528W'], 'EGFR': ['V802F','L844V'],
        'ESR1': ['L536R','L370F'], 'FLT3': ['Y842C']}
U=f'{_SAVDR_ROOT}/inputs'
def first_record(path):
    s,st=[],False
    for l in open(path):
        if l.startswith('>'):
            if st: break
            st=True; continue
        if st: s.append(l.strip())
    return ''.join(s)
def window(seq,pos,mt,L):
    i=pos-1; s=i-L//2; e=s+L
    if s<0 or e>len(seq): return ''
    p=list(seq); p[i]=mt
    return ''.join(p[s:e])
def read_pa_pi(op):
    """Returns (Pa for Resistant, Pa for Not resistant)."""
    if not os.path.exists(op): return (None,None)
    pa=pi=None
    for line in open(op):
        p=line.rstrip('\n').split('\t')
        if len(p)>=3:
            cls=p[2].strip()
            if cls=='Resistant': pa=float(p[0])
            elif cls=='Not resistant': pi=float(p[0])
    return (pa,pi)

rows=[]
for gene,(name,L,MNA,drug) in CFG.items():
    work=f'/home/claude/prosp/{gene}'; os.makedirs(work,exist_ok=True)
    for e in ('OLMPASSdoMNABase.exe','OLMPASSofSDF.exe'):
        if not os.path.exists(f'{work}/{e}'): shutil.copy(f'{_SAVDR_ROOT}/results/data/{e}',work)
    seq=first_record(f'{U}/fasta/{name}.faa')
    tab=pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/{name}.csv')
    ex=tab.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$')
    tab['pos']=ex[1].astype(int); tab['mt']=ex[2]
    tab['pep']=[window(seq,p,m,L) for p,m in zip(tab.pos,tab.mt)]
    tab['is_drug']=[drug.lower() in str(d).lower() for d in tab.drugs]
    tab=tab[(tab.pep!='') ]
    tab=tab[[Chem.MolFromSequence(p) is not None for p in tab.pep]]
    npos=int((tab.y.astype(bool)&tab.is_drug).sum())
    print(f'{gene}: L={L} MNA={MNA} drug={drug} -> {len(tab)} trainable variants, '
          f'{npos} positives for this drug', flush=True)
    if not os.path.exists(f'{work}/B.MSAS'):
        w=Chem.SDWriter(f'{work}/train.sdf')
        for _,r in tab.iterrows():
            m=Chem.MolFromSequence(r.pep); m.SetProp('_Name',r.aa)
            m.SetProp('Drug','Resistant' if (r.y and r.is_drug) else 'Not resistant')
            w.write(m)
        w.close()
        open(f'{work}/t.txt','w',newline='').write(
            f'ClassLimit     3\r\nBaseCreate     {MNA}  B\r\nBaseAddNewData train.sdf  Drug\r\n'
            f'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
        for _ in range(3):
            for st in glob.glob(f'{work}/B.*')+glob.glob(f'{work}/OLMPASS2*'):
                if not st.endswith('.exe'): os.remove(st)
            subprocess.run(['wine','OLMPASSdoMNABase.exe','t.txt'],cwd=work,env=WINE,
                           stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=1800)
            if os.path.exists(f'{work}/B.MSAS'): break
        else:
            print(f'  {gene}: MultiPASS declined to train at L={L} MNA={MNA}', flush=True)
            for v in TEST[gene]:
                rows.append(dict(gene=gene,variant=v,L=L,MNA=MNA,drug=drug,n_pos=npos,
                                 Pa=None,Pi=None,call='model would not train'))
            continue
    # score the prospective variants plus, as a positive control, the target's own known positives
    import re
    todo=[]
    for v in TEST[gene]:
        m=re.match(r'^([A-Z])(\d+)([A-Z])$',v)
        pep=window(seq,int(m.group(2)),m.group(3),L)
        todo.append((v,pep,'prospective'))
    for _,r in tab[tab.y.astype(bool)&tab.is_drug].iterrows():
        todo.append((r.aa.replace('p.',''),r.pep,'training positive'))
    names={}
    w=Chem.SDWriter(f'{work}/batch.sdf')
    for j,(v,pep,kind) in enumerate(todo):
        if not pep or Chem.MolFromSequence(pep) is None: continue
        mm=Chem.MolFromSequence(pep); mm.SetProp('_Name',v)
        fn=f'p{j:04d}.txt'; mm.SetProp('FileName',fn); names[(v,kind)]=fn; w.write(mm)
    w.close()
    subprocess.run(['wine','OLMPASSofSDF.exe','batch.sdf','/BB.MSAS','/C'],cwd=work,env=WINE,
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=1800)
    for (v,kind),fn in names.items():
        pa,pi=read_pa_pi(f'{work}/{fn}')
        rows.append(dict(gene=gene,variant=v,L=L,MNA=MNA,drug=drug,n_pos=npos,kind=kind,
            Pa=pa,Pi=pi,call=('resistant' if (pa is not None and pi is not None and pa>pi)
                              else ('not resistant' if pa is not None else 'no result'))))
D=pd.DataFrame(rows); D.to_csv(f'{_SAVDR_ROOT}/results/data/prospective_scored.csv',index=False)
print('\n=== prospective variants ===')
P=D[D.kind=='prospective'] if 'kind' in D else D
print(P[['gene','variant','L','MNA','n_pos','Pa','Pi','call']].to_string(index=False))
if 'kind' in D:
    print('\n=== the same models on their own training positives (sanity control) ===')
    for g,gg in D[D.kind=='training positive'].groupby('gene'):
        hit=(gg.call=='resistant').sum()
        print(f'  {g:5s} {hit}/{len(gg)} of its training positives called resistant; '
              f'mean Pa {gg.Pa.mean():.3f}')
