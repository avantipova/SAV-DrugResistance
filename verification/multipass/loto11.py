"""Leave-one-target-out across 11 proteins, common window length.

Train an MNA base on every variant of the other ten proteins, score every variant of the
held-out one. No folds and no seeds: the split is by protein.

Two modes so the long base training can run detached from the bounded scoring step:
    python3 loto11.py train <MNA> <held>
    python3 loto11.py score <MNA> <held>

Base and per-molecule checkpoints throughout; B.done is written only after B.MSAS exists, so a
killed run never leaves a half-trained base that looks usable. train.sdf is deleted on success
-- a truncated SDF left by a killed writer silently trained one earlier base on 4093 of 4152
molecules, which the verification caught by reading the substance count out of B.HST.
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

import os, sys, glob, json, shutil, subprocess, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
L = int(os.environ.get('L', 25))
TAB = f'/home/claude/tab11_L{L}'
TARGETS = json.load(open(f'{TAB}/_targets.json'))
MODE, MNA, HELD = sys.argv[1], int(sys.argv[2]), sys.argv[3]
WORK = f'/home/claude/loto11_L{L}_mna{MNA}_{HELD}'
EXE = f'{_SAVDR_ROOT}/results/data'


def table(t):
    return pd.read_csv(f'{TAB}/{t}.csv')


def mol(t, r, cls=None):
    m = Chem.MolFromSequence(r.pep)
    m.SetProp('_Name', f'{t}|{r.aa}')
    if cls: m.SetProp('Drug', cls)
    return m


os.makedirs(WORK, exist_ok=True)
for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofMOL.exe'):
    if not os.path.exists(f'{WORK}/{e}'): shutil.copy(f'{EXE}/{e}', WORK)

if MODE == 'train':
    if os.path.exists(f'{WORK}/B.done') and os.path.exists(f'{WORK}/B.MSAS'):
        print(f'{HELD}: base already trained'); raise SystemExit
    sdf = 'train.sdf'
    ntr = 0
    if not os.path.exists(f'{WORK}/{sdf}'):
        w = Chem.SDWriter(f'{WORK}/{sdf}')
        for t in TARGETS:
            if t == HELD: continue
            d = table(t)
            for _, r in d.iterrows():
                w.write(mol(t, r, 'Resistant' if r.y else 'Not resistant'))
                ntr += 1
        w.close()
    open(f'{WORK}/t.txt', 'w', newline='').write(
        f'ClassLimit     3\r\nBaseCreate     {MNA}  B\r\nBaseAddNewData {sdf}  Drug\r\n'
        f'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
    for stale in glob.glob(f'{WORK}/B.*') + glob.glob(f'{WORK}/OLMPASS2*'):
        if not stale.endswith('.exe'): os.remove(stale)
    print(f'{HELD}: training on {ntr} variants from {len(TARGETS)-1} proteins, MNA {MNA}', flush=True)
    subprocess.run(['wine', 'OLMPASSdoMNABase.exe', 't.txt'], cwd=WORK, env=WINE,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=10800)
    if not os.path.exists(f'{WORK}/B.MSAS'):
        raise SystemExit(f'{HELD}: B.MSAS never written')
    open(f'{WORK}/B.done', 'w').close()
    os.remove(f'{WORK}/{sdf}')
    print(f'{HELD}: base done, {os.path.getsize(f"{WORK}/B.MSAS")/1e6:.0f} MB', flush=True)
    raise SystemExit

# ---- score ----------------------------------------------------------------------------------
if not (os.path.exists(f'{WORK}/B.done') and os.path.exists(f'{WORK}/B.MSAS')):
    raise SystemExit(f'{HELD}: no trained base -- run train first')
te = table(HELD)
part = f'{WORK}/partial.tsv'
done = {}
if os.path.exists(part):
    for line in open(part):
        p = line.rstrip('\n').split('\t')
        if len(p) == 2: done[int(p[0])] = float(p[1])
todo = [i for i in range(len(te)) if i not in done]
print(f'{HELD}: {len(done)} cached, {len(todo)} to score', flush=True)


def predict(i):
    f, o = f'q{i}.sdf', f'q{i}.txt'
    w = Chem.SDWriter(f'{WORK}/{f}'); w.write(mol(HELD, te.iloc[i])); w.close()
    subprocess.run(['wine', 'OLMPASSofMOL.exe', f, o, '/BB.MSAS', '/C'], cwd=WORK, env=WINE,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    pa = np.nan
    op = f'{WORK}/{o}'
    if os.path.exists(op):
        for line in open(op):
            p = line.rstrip('\n').split('\t')
            if len(p) == 3 and p[2].strip() == 'Resistant':
                pa = float(p[0]); break
    os.remove(f'{WORK}/{f}')
    os.path.exists(op) and os.remove(op)
    return i, pa


scores = np.full(len(te), np.nan)
for i, v in done.items(): scores[i] = v
if todo:
    with open(part, 'a', buffering=1) as pf, ThreadPoolExecutor(int(os.environ.get('WORKERS', 4))) as ex:
        for i, pa in ex.map(predict, todo):
            scores[i] = pa
            pf.write(f'{i}\t{pa}\n')
if np.isnan(scores).any():
    print(f'{HELD}: {int(np.isnan(scores).sum())} still unscored', flush=True)
    raise SystemExit
te['multipass_pa'] = scores
te.to_csv(f'{WORK}/scores.csv', index=False)
from sklearn.metrics import roc_auc_score
print(f'{HELD}: LOTO AUC = {roc_auc_score(te.y, np.nan_to_num(scores)):.3f} '
      f'on {int(te.y.sum())}+/{int((1-te.y).sum())}-', flush=True)
