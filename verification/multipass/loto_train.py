"""Train one leave-one-target-out MNA base. Separated from scoring because a cross-target base
takes longer than a single foreground call allows, so this is run detached and polled.

Writes B.done only after B.MSAS exists, so a killed run never leaves a half-written base that
looks usable. Refuses to touch an existing finished base.
"""
import os, sys, glob, shutil, subprocess, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
TARGETS = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
           'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
MNA, held = int(sys.argv[1]), sys.argv[2]
work = f'loto_{held}_mna{MNA}'
os.makedirs(work, exist_ok=True)
for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofMOL.exe'):
    if not os.path.exists(f'{work}/{e}'): shutil.copy(e, work)

if os.path.exists(f'{work}/B.done') and os.path.exists(f'{work}/B.MSAS'):
    print(f'{held}: base already trained'); raise SystemExit

frames = []
for t, (L, M) in TARGETS.items():
    d = pd.read_csv(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_{L}_{M}_s0/scores.csv')
    d['target'] = t
    frames.append(d[['target', 'aa', 'pos', 'pep', 'y']])
D = pd.concat(frames, ignore_index=True)
tr = D[D.target != held].reset_index(drop=True)

sdf = 'train.sdf'
if not os.path.exists(f'{work}/{sdf}'):
    w = Chem.SDWriter(f'{work}/{sdf}')
    for _, r in tr.iterrows():
        m = Chem.MolFromSequence(r.pep)
        m.SetProp('_Name', f'{r.target}|{r.aa}')
        m.SetProp('Drug', 'Resistant' if r.y else 'Not resistant')
        w.write(m)
    w.close()

task = 't.txt'
open(f'{work}/{task}', 'w', newline='').write(
    f'ClassLimit     3\r\nBaseCreate     {MNA}  B\r\nBaseAddNewData {sdf}  Drug\r\n'
    f'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')

for stale in glob.glob(f'{work}/B.*') + glob.glob(f'{work}/OLMPASS2*'):
    if not stale.endswith('.exe'): os.remove(stale)
print(f'{held}: training on {len(tr)} variants ({int(tr.y.sum())} resistant) '
      f'from {tr.target.nunique()} targets, MNA {MNA}', flush=True)
subprocess.run(['wine', 'OLMPASSdoMNABase.exe', task], cwd=work, env=WINE,
               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=7200)
if not os.path.exists(f'{work}/B.MSAS'):
    raise SystemExit(f'{held}: B.MSAS never written')
open(f'{work}/B.done', 'w').close()
os.remove(f'{work}/{sdf}')          # 86 MB, reconstructible
print(f'{held}: base done, {os.path.getsize(f"{work}/B.MSAS")/1e6:.0f} MB', flush=True)
