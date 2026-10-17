"""MultiPASS evaluated on the SAME position-grouped folds as the baseline models.

Binary task per target: 'Resistant' (to any drug) vs 'Not resistant'.
For each fold: train an MSAR/MSAS base on the training variants only (OLMPASSdoMNABase),
then score every held-out variant with OLMPASSofMOL /C and take Pa('Resistant').
Nothing from the test fold ever enters the base, so the AUC is out-of-fold and directly
comparable with 02_baselines.py (same split, same seed, same variants).
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

import os, sys, glob, subprocess, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem, RDLogger
from sklearn.model_selection import StratifiedGroupKFold
RDLogger.DisableLog('rdApp.*')

U = f'{_SAVDR_ROOT}/external/source-data'
WINE = dict(os.environ, WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
NAME, L, MNA = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
REGION = 'region' in sys.argv[4:]
SEED = int([a for a in sys.argv[4:] if a.isdigit()][0]) if any(a.isdigit() for a in sys.argv[4:]) else 0
WORK = f'fold_{NAME}_{L}_{MNA}' + ('_region' if REGION else '') + f'_s{SEED}'; os.makedirs(WORK, exist_ok=True)
import shutil
for e in ('OLMPASSdoMNABase.exe','OLMPASSofMOL.exe'): shutil.copy(e, WORK)

def first_record(path):
    """FIX: KIT (8 records), SMO (3) and FLT3 (11) are multi-record FASTA files. Reading them
    with a plain "skip the > lines" concatenates every isoform into one string, which lets a
    peptide window run past the end of the canonical sequence and into the next record's
    N-terminus instead of being rejected. Only the first record is the canonical sequence."""
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

seq = first_record(f'{U}/fasta/{NAME}.faa')
g = pd.read_csv(f'{U}/csv/{NAME}.csv'); g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
res = set(pd.read_csv(f'{U}/mutations/{NAME}_mutations.csv')['AA Mutation'])

rows = []
for aa in g['AA Mutation'].dropna().unique():
    pep = extract(seq, aa, L)
    if not pep or Chem.MolFromSequence(pep) is None: continue
    rows.append(dict(aa=aa, pos=int(''.join(c for c in aa.split('.')[1] if c.isdigit())),
                     pep=pep, y=int(aa in res)))
df = pd.DataFrame(rows)
if REGION:
    lo, hi = df.loc[df.y == 1, 'pos'].min() - 10, df.loc[df.y == 1, 'pos'].max() + 10
    df = df[(df.pos >= lo) & (df.pos <= hi)].reset_index(drop=True)
y, grp = df.y.values, df.pos.values
print(f'{NAME} L={L} MNA={MNA} region={REGION}: {y.sum()} pos / {(1-y).sum()} neg, {df.pos.nunique()} positions', flush=True)

def mol(r, cls=None):
    m = Chem.MolFromSequence(r.pep); m.SetProp('_Name', r.aa)
    if cls: m.SetProp('Drug', cls)
    return m

SHM = os.environ.get('SHM') == '1'
def BASEREF(base):
    if not SHM: return base
    dst = f'/dev/shm/{WORK}_{base}'
    if not os.path.exists(dst):
        import shutil as _sh; _sh.copy(f'{WORK}/{base}', dst + '.tmp'); os.replace(dst + '.tmp', dst)
    return 'Z:' + dst.replace('/', chr(92))

def predict(args):
    i, r, base = args
    f, o = f'q{i}.sdf', f'q{i}.txt'
    w = Chem.SDWriter(f'{WORK}/{f}'); w.write(mol(r)); w.close()
    subprocess.run(['wine', 'OLMPASSofMOL.exe', f, o, f'/B{BASEREF(base)}', '/C'], cwd=WORK,
                   env=WINE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    pa = np.nan
    o = f'{WORK}/{o}'
    if os.path.exists(o):
        for line in open(o):
            p = line.rstrip('\n').split('\t')
            if len(p) == 3 and p[2].strip() == 'Resistant':
                pa = float(p[0]); break
    os.remove(f'{WORK}/{f}'); os.path.exists(o) and os.remove(o)
    return i, pa

scores = np.full(len(df), np.nan)
for k, (tr, te) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(df, y, grp)):
    ck = f'{WORK}/fold{k}_scores.npy'
    if os.path.exists(ck):                      # already done before a restart
        scores[te] = np.load(ck); print(f'  fold {k}: resumed from checkpoint', flush=True); continue
    base = f'B{k}'
    sdf = f'train{k}.sdf'; w = Chem.SDWriter(f'{WORK}/{sdf}')
    for _, r in df.iloc[tr].iterrows(): w.write(mol(r, 'Resistant' if r.y else 'Not resistant'))
    w.close()
    task = f't{k}.txt'
    open(f'{WORK}/{task}', 'w', newline='').write(f'ClassLimit     3\r\nBaseCreate     {MNA}  {base}\r\n'
        f'BaseAddNewData {sdf}  Drug\r\nBaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
    # EGFR's base takes ~5 min to train and 289 MB on disk; a mid-fold restart must not redo it.
    # The marker is written only after .MSAS appeared, so a half-written base is never reused.
    trained = os.path.exists(f'{WORK}/{base}.done') and os.path.exists(f'{WORK}/{base}.MSAS')
    if trained:
        print(f'  fold {k}: reusing cached base {base}.MSAS', flush=True)
    else:
        for attempt in range(3):
            # MultiPASS auto-renames the base when files of that name are already there (leftovers from a
            # killed run), so BaseSaveAsSAS lands on an auto-named .MSAS and the fold scores nothing.
            for stale in glob.glob(f'{WORK}/{base}.*') + glob.glob(f'{WORK}/OLMPASS2*'):
                if not stale.endswith('.exe'): os.remove(stale)   # never touch the MultiPASS binaries
            subprocess.run(['wine', 'OLMPASSdoMNABase.exe', task], cwd=WORK, env=WINE,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
            if os.path.exists(f'{WORK}/{base}.MSAS'): break
            print(f'  fold {k}: {base}.MSAS missing after training, retry {attempt+1}', flush=True)
        else:
            raise RuntimeError(f'{base}.MSAS never written -- refusing to score against a missing base')
        open(f'{WORK}/{base}.done', 'w').close()
    # per-molecule checkpoint: the container is reclaimed a couple of minutes after each turn,
    # so a fold that takes longer than that has to survive being killed mid-scoring.
    part = f'{WORK}/fold{k}_partial.tsv'
    done = {}
    if os.path.exists(part):
        for line in open(part):
            p = line.split('\t')
            if len(p) == 2: done[int(p[0])] = float(p[1])
    todo = [i for i in te if i not in done]
    for i, v in done.items(): scores[i] = v
    if todo: print(f'  fold {k}: {len(done)} cached, {len(todo)} to score', flush=True)
    with open(part, 'a', buffering=1) as pf, ThreadPoolExecutor(int(os.environ.get('WORKERS', 8))) as ex:
        for i, pa in ex.map(predict, [(i, df.iloc[i], base + '.MSAS') for i in todo]):
            scores[i] = pa
            pf.write(f'{i}\t{pa}\n')
    np.save(ck, scores[te])
    print(f'  fold {k}: trained on {len(tr)}, scored {len(te)}, missing {np.isnan(scores[te]).sum()}', flush=True)

df['multipass_pa'] = scores
df.to_csv(f'{WORK}/scores.csv', index=False)
ok = ~np.isnan(scores)
from sklearn.metrics import roc_auc_score, average_precision_score
print('MultiPASS out-of-fold AUC', round(roc_auc_score(y[ok], scores[ok]), 3),
      'AUPRC', round(average_precision_score(y[ok], scores[ok]), 3),
      f'on {ok.sum()}/{len(y)} variants', flush=True)
