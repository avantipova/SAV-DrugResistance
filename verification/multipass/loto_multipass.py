"""Leave-one-target-out for MultiPASS: train on five proteins, score the held-out sixth.

This is the analysis `08_loto.csv` was supposed to be. That file was stale (written 3 Oct,
before the multi-record FASTA fix of 4 Oct, with KIT at 943 negatives instead of 936 and SMO at
464 instead of 453), its EGFR row silently fell back to the `region` variant table, and it had
no MultiPASS row at all.

Design. One MNA base per held-out target, trained on every variant of the other five targets,
labelled Resistant / Not resistant exactly as build_sdf.py does it; then every variant of the
held-out target is scored with OLMPASSofMOL /C and Pa('Resistant') is read off.

No cross-validation and no seeds: the split is by protein, so there is nothing random to
average over. The peptides keep each target's own window length (17/33/37/23/25/31) because
that is what the main analysis used, and a length cue cannot inflate the held-out AUC -- every
test molecule of a given held-out target has the same length, so length carries no information
*within* the target, which is where the AUC is computed. The MNA level has to be a single
choice for a cross-target base; 9 is used (the level ABL1 and KIT use), and MNA 4 is run as a
sensitivity check.

`nearest_pos` and `pos_only` are reported but are not meaningful here by construction: residue
numbering is not comparable across proteins. They are the floor, not a competitor.
"""
import os, sys, glob, shutil, subprocess, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
TARGETS = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
           'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
MNA = int(sys.argv[1]) if len(sys.argv) > 1 else 9
HELD = sys.argv[2] if len(sys.argv) > 2 else None

frames = []
for t, (L, M) in TARGETS.items():
    f = f'fold_{t}_{L}_{M}_s0/scores.csv'
    if not glob.glob(f):
        raise SystemExit(f'missing {f} -- the full-setting run must exist for every target')
    d = pd.read_csv(f)
    d['target'] = t
    frames.append(d[['target', 'aa', 'pos', 'pep', 'y']])
D = pd.concat(frames, ignore_index=True)
print('variants per target:', D.groupby('target').size().to_dict(), flush=True)
print('positives per target:', D[D.y == 1].groupby('target').size().to_dict(), flush=True)


def mol(r, cls=None):
    m = Chem.MolFromSequence(r.pep)
    m.SetProp('_Name', f'{r.target}|{r.aa}')
    if cls: m.SetProp('Drug', cls)
    return m


def predict(args):
    i, r, work, base = args
    f, o = f'q{i}.sdf', f'q{i}.txt'
    w = Chem.SDWriter(f'{work}/{f}'); w.write(mol(r)); w.close()
    subprocess.run(['wine', 'OLMPASSofMOL.exe', f, o, f'/B{base}', '/C'], cwd=work, env=WINE,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
    pa = np.nan
    op = f'{work}/{o}'
    if os.path.exists(op):
        for line in open(op):
            p = line.rstrip('\n').split('\t')
            if len(p) == 3 and p[2].strip() == 'Resistant':
                pa = float(p[0]); break
    os.remove(f'{work}/{f}')
    os.path.exists(op) and os.remove(op)
    return i, pa


for held in ([HELD] if HELD else list(TARGETS)):
    work = f'loto_{held}_mna{MNA}'
    os.makedirs(work, exist_ok=True)
    for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofMOL.exe'):
        if not os.path.exists(f'{work}/{e}'): shutil.copy(e, work)
    out = f'{work}/scores.csv'
    if os.path.exists(out):
        print(f'{held}: already done'); continue

    tr = D[D.target != held].reset_index(drop=True)
    te = D[D.target == held].reset_index(drop=True)
    base = 'B'
    if not (os.path.exists(f'{work}/{base}.done') and os.path.exists(f'{work}/{base}.MSAS')):
        print(f'{held}: no trained base yet -- run loto_train.py {MNA} {held} first', flush=True)
        continue

    part = f'{work}/partial.tsv'
    done = {}
    if os.path.exists(part):
        for line in open(part):
            p = line.split('\t')
            if len(p) == 2: done[int(p[0])] = float(p[1])
    todo = [i for i in range(len(te)) if i not in done]
    print(f'{held}: {len(done)} cached, {len(todo)} to score', flush=True)
    scores = np.full(len(te), np.nan)
    for i, v in done.items(): scores[i] = v
    with open(part, 'a', buffering=1) as pf, ThreadPoolExecutor(int(os.environ.get('WORKERS', 8))) as ex:
        for i, pa in ex.map(predict, [(i, te.iloc[i], work, base + '.MSAS') for i in todo]):
            scores[i] = pa
            pf.write(f'{i}\t{pa}\n')
    if np.isnan(scores).any():
        print(f'{held}: still {int(np.isnan(scores).sum())} unscored, not writing scores.csv',
              flush=True)
        continue
    te['multipass_pa'] = scores
    te.to_csv(out, index=False)
    from sklearn.metrics import roc_auc_score
    print(f'{held}: LOTO MultiPASS AUC = {roc_auc_score(te.y, scores):.3f} '
          f'on {int(te.y.sum())}+/{int((1-te.y).sum())}-', flush=True)
