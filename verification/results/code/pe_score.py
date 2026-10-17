"""Score every missense variant of the ABL1 saturation prime-editing scan with the
leave-one-target-out MultiPASS base trained WITHOUT ABL1 (loto11_L25_mna9_ABL1_*/B.MSAS).

The model has never seen an ABL1 variant, so the full 1848-variant scan is a genuine
external test set -- and unlike COSMIC-derived sets it is not enriched for observed
variants, because prime editing made all of them by design.
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

import os, shutil, subprocess, json, pandas as pd, numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

L = 25
LOTO = '/home/claude/loto11_L25_mna9_ABL1_ENST00000318560'
WORK = '/home/claude/pe_score'
WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
os.makedirs(WORK, exist_ok=True)
for f in ('B.MSAS', 'OLMPASSofSDF.exe'):
    src = f'{LOTO}/{f}' if os.path.exists(f'{LOTO}/{f}') else f'{_SAVDR_ROOT}/results/data/{f}'
    if not os.path.exists(f'{WORK}/{f}'): shutil.copy(src, WORK)

seq = ''.join(l.strip() for l in
              open(f'{_SAVDR_ROOT}/inputs/fasta/ABL1_ENST00000318560.faa') if not l.startswith('>'))

def window(pos, mt, length=L):
    i = pos - 1
    s = i - length // 2; e = s + length
    if s < 0 or e > len(seq): return ''
    p = list(seq); p[i] = mt
    return ''.join(p[s:e])

sc = pd.read_csv(f'{_SAVDR_ROOT}/external/primeedit/RS_Imatinib.csv')
sc = sc[sc.mut_type == 'Missense'].copy()
ex = sc.AA_var.str.extract(r'^([A-Z])(\d+)([A-Z])$')
sc['wt'], sc['pos'], sc['mt'] = ex[0], ex[1].astype(int), ex[2]
sc['pep'] = [window(p, m) for p, m in zip(sc.pos, sc.mt)]
bad = sc[sc.pep == '']
print(f'scan missense {len(sc)}; dropped for window running off the protein: {len(bad)}')
sc = sc[sc.pep != ''].copy()
sc = sc[[Chem.MolFromSequence(p) is not None for p in sc.pep]]
print(f'scorable {len(sc)}')

# reproduce our own pep for the overlap, as a check that the window is built identically
ours = pd.read_csv(f'{LOTO}/scores.csv')
oe = ours.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$')
ours['wt'], ours['pos'], ours['mt'] = oe[0], oe[1].astype(int), oe[2]
ours['rep'] = [window(p, m) for p, m in zip(ours.pos, ours.mt)]
n_ok = int((ours.rep == ours.pep).sum())
print(f'our pep reproduced from the FASTA: {n_ok}/{len(ours)}')
assert n_ok == len(ours), ours[ours.rep != ours.pep][['aa','pep','rep']].head()

part = f'{WORK}/partial.tsv'
done = {}
if os.path.exists(part):
    for line in open(part):
        p = line.rstrip('\n').split('\t')
        if len(p) == 2: done[p[0]] = float(p[1])

def read_pa(op):
    if not os.path.exists(op): return None
    for line in open(op):
        p = line.rstrip('\n').split('\t')
        if len(p) == 3 and p[2].strip() == 'Resistant':
            return float(p[0])
    return 0.0

todo = [(k, p) for k, p in zip(sc.AA_var, sc.pep) if k not in done]
print(f'to score: {len(todo)}', flush=True)
CH = 250
with open(part, 'a', buffering=1) as pf:
    for i in range(0, len(todo), CH):
        chunk = todo[i:i + CH]
        names = {}
        w = Chem.SDWriter(f'{WORK}/batch.sdf')
        for j, (k, pep) in enumerate(chunk):
            m = Chem.MolFromSequence(pep); m.SetProp('_Name', k)
            fn = f'm{i + j:06d}.txt'; m.SetProp('FileName', fn); names[k] = fn
            w.write(m)
        w.close()
        subprocess.run(['wine', 'OLMPASSofSDF.exe', 'batch.sdf', '/BB.MSAS', '/C'], cwd=WORK,
                       env=WINE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=1800)
        miss = 0
        for k, fn in names.items():
            v = read_pa(f'{WORK}/{fn}')
            if v is None: v, miss = 0.0, miss + 1
            done[k] = v; pf.write(f'{k}\t{v}\n')
            os.path.exists(f'{WORK}/{fn}') and os.remove(f'{WORK}/{fn}')
        print(f'  {i + len(chunk)}/{len(todo)} (no result file: {miss})', flush=True)

out = pd.DataFrame({'AA_var': list(done), 'multipass_pa': list(done.values())})
out.to_csv(f'{_SAVDR_ROOT}/results/data/pe_scan_pa.csv', index=False)
print('wrote', len(out), 'scores')
