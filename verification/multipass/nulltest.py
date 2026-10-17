"""Shuffled-label null for MultiPASS: the identical fold pipeline on permuted labels.
Usage: nulltest.py <target> <L> <MNA> [region] <n_perm>"""
import os, sys, subprocess, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem, RDLogger
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog('rdApp.*')
import importlib.util
spec = importlib.util.spec_from_file_location('f', 'folds.py')
NAME, L, MNA = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
REGION = 'region' in sys.argv[4:]
NPERM = int([a for a in sys.argv[4:] if a.isdigit()][0])
src = f'fold_{NAME}_{L}_{MNA}' + ('_region' if REGION else '') + '_s0/scores.csv'
d = pd.read_csv(src)
WORK = f'null_{NAME}_{L}_{MNA}' + ('_region' if REGION else ''); os.makedirs(WORK, exist_ok=True)
import shutil
for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofMOL.exe'): shutil.copy(e, WORK)
WINE = dict(os.environ, WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
def mol(r, cls=None):
    m = Chem.MolFromSequence(r.pep); m.SetProp('_Name', str(r.aa))
    if cls: m.SetProp('Drug', cls)
    return m
def predict(a):
    i, r, base = a
    f, o = f'q{i}.sdf', f'q{i}.txt'
    w = Chem.SDWriter(f'{WORK}/{f}'); w.write(mol(r)); w.close()
    subprocess.run(['wine', 'OLMPASSofMOL.exe', f, o, f'/B{base}', '/C'], cwd=WORK, env=WINE,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    pa, o = np.nan, f'{WORK}/{o}'
    if os.path.exists(o):
        for line in open(o):
            p = line.rstrip('\n').split('\t')
            if len(p) == 3 and p[2].strip() == 'Resistant': pa = float(p[0]); break
    os.remove(f'{WORK}/{f}'); os.path.exists(o) and os.remove(o)
    return i, pa
rng = np.random.default_rng(11); aucs = []
for perm in range(NPERM):
    y = rng.permutation(d.y.values); s = np.full(len(d), np.nan)
    for k, (tr, te) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(d, y, d.pos.values)):
        base = f'N{k}'
        sdf = f'tr{k}.sdf'; w = Chem.SDWriter(f'{WORK}/{sdf}')
        for j in tr: w.write(mol(d.iloc[j], 'Resistant' if y[j] else 'Not resistant'))
        w.close()
        open(f'{WORK}/t{k}.txt', 'w', newline='').write(
            f'ClassLimit     3\r\nBaseCreate     {MNA}  {base}\r\nBaseAddNewData {sdf}  Drug\r\n'
            f'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
        subprocess.run(['wine', 'OLMPASSdoMNABase.exe', f't{k}.txt'], cwd=WORK, env=WINE,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
        with ThreadPoolExecutor(8) as ex:
            for i, pa in ex.map(predict, [(i, d.iloc[i], base + '.MSAS') for i in te]): s[i] = pa
    ok = ~np.isnan(s); aucs.append(roc_auc_score(y[ok], s[ok]))
    print(f'perm {perm}: AUC {aucs[-1]:.3f}', flush=True)
obs = roc_auc_score(d.y.values, d.multipass_pa.values)
print(f'observed {obs:.3f} | null mean {np.mean(aucs):.3f} max {np.max(aucs):.3f} '
      f'p={(1 + sum(a >= obs for a in aucs)) / (len(aucs) + 1):.3f}')
pd.DataFrame(dict(perm=range(len(aucs)), auc=aucs)).to_csv(f'{WORK}/null.csv', index=False)
