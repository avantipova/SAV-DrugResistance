"""MultiPASS with the project's own PER-DRUG multiclass labels, on the same position-grouped folds.

The binary run (folds.py) collapses every drug into "Resistant". Table 1 of the paper does not:
it reports one IAP per drug. This script writes the drug list into the SDF class field exactly
as build_sdf.py does (newline-separated, "Not resistant" otherwise), so one base per fold
predicts Pa for every drug at once, and the per-drug AUC is read from that.

For a given drug: positives = variants annotated resistant TO THAT DRUG; negatives = variants
with no resistance annotation at all. Variants resistant only to a different drug are dropped
rather than counted as negatives, since their status for this drug is unknown.

Original header follows.
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


"""MultiPASS evaluated on the SAME position-grouped folds as the baseline models.

Binary task per target: 'Resistant' (to any drug) vs 'Not resistant'.
For each fold: train an MSAR/MSAS base on the training variants only (OLMPASSdoMNABase),
then score every held-out variant with OLMPASSofMOL /C and take Pa('Resistant').
Nothing from the test fold ever enters the base, so the AUC is out-of-fold and directly
comparable with 02_baselines.py (same split, same seed, same variants).
"""
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
WORK = f'drug_{NAME}_{L}_{MNA}' + ('_region' if REGION else '') + f'_s{SEED}'; os.makedirs(WORK, exist_ok=True)
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
_m = pd.read_csv(f'{U}/mutations/{NAME}_mutations.csv')
res = set(_m['AA Mutation'])
DRUGS = {}
for _aa, _d in zip(_m['AA Mutation'], _m['Drug']):
    DRUGS.setdefault(_aa, set()).add(str(_d))

rows = []
for aa in g['AA Mutation'].dropna().unique():
    pep = extract(seq, aa, L)
    if not pep or Chem.MolFromSequence(pep) is None: continue
    rows.append(dict(aa=aa, pos=int(''.join(c for c in aa.split('.')[1] if c.isdigit())),
                     pep=pep, y=int(aa in res),
                     drugs=';'.join(sorted(DRUGS.get(aa, [])))))
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

def drugclass(r):
    """multiclass label exactly as the project's own SDF builder writes it:
    newline-separated drug names for a resistant variant, else 'Not resistant'"""
    return '\n'.join(r.drugs.split(';')) if r.drugs else 'Not resistant'

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
    pa = {}
    o = f'{WORK}/{o}'
    if os.path.exists(o):
        for line in open(o):
            p = line.rstrip('\n').split('\t')
            if len(p) == 3:
                try: pa[p[2].strip()] = float(p[0])
                except ValueError: pass
    os.remove(f'{WORK}/{f}'); os.path.exists(o) and os.remove(o)
    return i, pa

scores = [None] * len(df)
for k, (tr, te) in enumerate(StratifiedGroupKFold(5, shuffle=True, random_state=SEED).split(df, y, grp)):
    ck = f'{WORK}/fold{k}_scores.json'
    if os.path.exists(ck):                      # already done before a restart
        import json as _j
        for _i, _v in _j.load(open(ck)).items(): scores[int(_i)] = _v
        print(f'  fold {k}: resumed from checkpoint', flush=True); continue
    base = f'B{k}'
    sdf = f'train{k}.sdf'; w = Chem.SDWriter(f'{WORK}/{sdf}')
    for _, r in df.iloc[tr].iterrows(): w.write(mol(r, drugclass(r)))
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
    import json as _j
    part = f'{WORK}/fold{k}_partial.tsv'
    done = {}
    if os.path.exists(part):
        for line in open(part):
            p = line.rstrip('\n').split('\t', 1)
            if len(p) == 2:
                try: done[int(p[0])] = _j.loads(p[1])
                except Exception: pass
    todo = [i for i in te if i not in done]
    for i, v in done.items(): scores[i] = v
    if todo: print(f'  fold {k}: {len(done)} cached, {len(todo)} to score', flush=True)
    with open(part, 'a', buffering=1) as pf, ThreadPoolExecutor(int(os.environ.get('WORKERS', 8))) as ex:
        for i, pa in ex.map(predict, [(i, df.iloc[i], base + '.MSAS') for i in todo]):
            scores[i] = pa
            pf.write(f'{i}\t' + _j.dumps(pa) + '\n')
    _j.dump({int(i): scores[i] for i in te}, open(ck, 'w'))
    print(f'  fold {k}: trained on {len(tr)}, scored {len(te)}, '
          f'missing {sum(scores[i] is None for i in te)}', flush=True)

acts = sorted({a for d in scores if d for a in d})
for a in acts:
    df['pa::' + a] = [np.nan if d is None else d.get(a, np.nan) for d in scores]
df.to_csv(f'{WORK}/scores_drug.csv', index=False)
from sklearn.metrics import roc_auc_score
print('activities predicted:', ', '.join(acts), flush=True)
for a in acts:
    if a == 'Not resistant': continue
    y_a = np.array([int(a in (r.drugs.split(';') if r.drugs else [])) for r in df.itertuples()])
    neg = (df.y.values == 0)
    keep = (y_a == 1) | neg                 # drop variants resistant only to OTHER drugs
    s_a = df['pa::' + a].values
    ok = keep & np.isfinite(s_a)
    if y_a[ok].sum() >= 3 and (1 - y_a[ok]).sum() >= 3:
        print(f'   {a:38s} n+={int(y_a[ok].sum()):4d} n-={int((1-y_a[ok]).sum()):5d} '
              f'AUC={roc_auc_score(y_a[ok], s_a[ok]):.3f}', flush=True)
