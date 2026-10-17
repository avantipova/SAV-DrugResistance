"""Within-target position-grouped CV with per-variant predictions saved.

The ActivArk comparison in REPORT_ACTIVARK_3TARGETS.md used leave-one-target-out scores, i.e. a
model that never saw the target, against ActivArk, which was trained on 441 kinases including
these three. This recomputes our side with the model the web platform actually serves: trained
on the target itself, scored out-of-fold.

Everything except what the model was trained on is held identical to the LOTO run -- L = 25,
MNA = 9, the same executable, the same peptide windows -- so the only variable is target
awareness. Folds are StratifiedGroupKFold(5) grouped by residue position, as elsewhere.
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
from rdkit import Chem, RDLogger
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog('rdApp.*')
WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
NAME = sys.argv[1]
ROOT = f'/home/claude/wt_{NAME}'
os.makedirs(ROOT, exist_ok=True)
for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofSDF.exe'):
    if not os.path.exists(f'{ROOT}/{e}'): shutil.copy(f'{_SAVDR_ROOT}/results/data/{e}', ROOT)

tab = pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/{NAME}.csv').set_index('aa')
Y, G = tab.y.values, tab.pos.values
AA = np.array(tab.index)
print(f'{NAME}: {len(tab)} variants, {int(Y.sum())} positives, '
      f'{len(set(G))} positions', flush=True)

def mol(aa, pep, cls=None):
    m = Chem.MolFromSequence(pep); m.SetProp('_Name', aa)
    if cls: m.SetProp('Drug', cls)
    return m

def read_pa(op):
    if not os.path.exists(op): return None
    for line in open(op):
        p = line.rstrip('\n').split('\t')
        if len(p) == 3 and p[2].strip() == 'Resistant': return float(p[0])
    return 0.0

def block(k, train_aa, test_aa):
    ck = f'{ROOT}/fold{k}.json'
    if os.path.exists(ck): return json.load(open(ck))
    work = f'{ROOT}/f{k}'; os.makedirs(work, exist_ok=True)
    for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofSDF.exe'):
        if not os.path.exists(f'{work}/{e}'): shutil.copy(f'{ROOT}/{e}', work)
    if not (os.path.exists(f'{work}/B.done') and os.path.exists(f'{work}/B.MSAS')):
        w = Chem.SDWriter(f'{work}/train.sdf')
        for aa in train_aa:
            r = tab.loc[aa]
            w.write(mol(aa, r.pep, 'Resistant' if r.y else 'Not resistant'))
        w.close()
        open(f'{work}/t.txt', 'w', newline='').write(
            'ClassLimit     3\r\nBaseCreate     9  B\r\nBaseAddNewData train.sdf  Drug\r\n'
            'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
        for _ in range(3):
            for stale in glob.glob(f'{work}/B.*') + glob.glob(f'{work}/OLMPASS2*'):
                if not stale.endswith('.exe'): os.remove(stale)
            subprocess.run(['wine', 'OLMPASSdoMNABase.exe', 't.txt'], cwd=work, env=WINE,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
            if os.path.exists(f'{work}/B.MSAS'): break
        else:
            json.dump({'__failed__': True}, open(ck, 'w')); return {'__failed__': True}
        open(f'{work}/B.done', 'w').close()
    out = {}
    CH = 400
    for i in range(0, len(test_aa), CH):
        chunk = test_aa[i:i + CH]
        names = {}
        w = Chem.SDWriter(f'{work}/batch.sdf')
        for j, aa in enumerate(chunk):
            m = mol(aa, tab.loc[aa].pep); fn = f'm{i + j:06d}.txt'
            m.SetProp('FileName', fn); names[aa] = fn; w.write(m)
        w.close()
        subprocess.run(['wine', 'OLMPASSofSDF.exe', 'batch.sdf', '/BB.MSAS', '/C'], cwd=work,
                       env=WINE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
        for aa, fn in names.items():
            v = read_pa(f'{work}/{fn}')
            out[aa] = 0.0 if v is None else v
            os.path.exists(f'{work}/{fn}') and os.remove(f'{work}/{fn}')
        os.path.exists(f'{work}/batch.sdf') and os.remove(f'{work}/batch.sdf')
    json.dump(out, open(ck, 'w'))
    for f in glob.glob(f'{work}/B.MS*') + glob.glob(f'{work}/B.MBFA') + glob.glob(f'{work}/train.sdf'):
        os.remove(f)
    return out

pred = {}
for k, (tr, te) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(tab, Y, G)):
    r = block(k, AA[tr], AA[te])
    if r.get('__failed__'):
        print(f'  fold {k}: MultiPASS declined to train', flush=True); continue
    pred.update(r)
    print(f'  fold {k}: {len(tr)} train / {len(te)} test done', flush=True)
S = pd.DataFrame({'aa': list(pred), 'wt_pa': list(pred.values())}).set_index('aa')
S = S.join(tab[['pos', 'y', 'drugs']])
S.to_csv(f'{_SAVDR_ROOT}/results/data/within_target_{NAME}.csv')
print(f'{NAME}: saved {len(S)} out-of-fold predictions; '
      f'whole-protein AUC {roc_auc_score(S.y, S.wt_pa):.4f}', flush=True)
