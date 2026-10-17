"""How much of the published AUC is selection optimism?

REPORT.md's caveat: "Table 1 reports the CV value of the configuration chosen by maximising
that same CV." The honest number re-does the (window length L, MNA level) search *inside* each
training fold, so the configuration is never chosen using the variants it is scored on.

  flat_max      best config by the outer CV itself -- the quantity Table 1 reports
  nested        per outer fold, the config is chosen by an inner CV on that fold's training
                variants only, then retrained on the whole training fold and applied to the
                held-out fold. The configs differ between outer folds; that is the point.
  optimism      flat_max - nested

Grid: L in {17, 25, 33} x MNA in {4, 9}. Outer 5 folds, inner 3 folds, both
StratifiedGroupKFold grouped by residue position. One outer seed (0): the question here is
selection bias, not fold-assignment variance, which the main table already reports as sd.

Every (config, fold) base and every molecule score is checkpointed, because a single
foreground call cannot finish this.
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
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog('rdApp.*')

# The data root is configurable because the five targets added later (BTK, ESR1, FLT3, IDH1,
# MET) live in the merged 11-target ingest, not in the original six-target upload.
U = os.environ.get('DATA_ROOT', f'{_SAVDR_ROOT}/external/source-data')
WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
GRID = [(L, M) for L in (17, 25, 33) for M in (4, 9)]
NAME = sys.argv[1]
ROOT = os.environ.get('NESTED_ROOT', f'nested_{NAME}')
os.makedirs(ROOT, exist_ok=True)
for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofMOL.exe', 'OLMPASSofSDF.exe'):
    if not os.path.exists(f'{ROOT}/{e}'): shutil.copy(e, ROOT)


def first_record(path):
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


SEQ = first_record(f'{U}/fasta/{NAME}.faa')
g = pd.read_csv(f'{U}/csv/{NAME}.csv')
g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
RES = set(pd.read_csv(f'{U}/mutations/{NAME}_mutations.csv')['AA Mutation'])
AAL = [a for a in g['AA Mutation'].dropna().unique()]


def table(L):
    """Variant table at window length L. The variant SET depends on L (windows that run off
    either end of the protein are dropped), so every comparison below is made on the
    intersection across the grid -- otherwise configs would be scored on different variants."""
    rows = []
    for aa in AAL:
        pep = extract(SEQ, aa, L)
        if not pep or Chem.MolFromSequence(pep) is None: continue
        rows.append(dict(aa=aa, pos=int(''.join(c for c in aa.split('.')[1] if c.isdigit())),
                         pep=pep, y=int(aa in RES)))
    return pd.DataFrame(rows).set_index('aa')


TABS = {L: table(L) for L in sorted({L for L, _ in GRID})}
COMMON = sorted(set.intersection(*[set(t.index) for t in TABS.values()]))
print(f'{NAME}: variants per L ' + ', '.join(f'L{L}={len(t)}' for L, t in TABS.items())
      + f'  -> common set {len(COMMON)}', flush=True)
BASE = TABS[min(TABS)].loc[COMMON]
Y, G = BASE.y.values, BASE.pos.values
print(f'{NAME}: {Y.sum()} pos / {(1-Y).sum()} neg on the common set', flush=True)


def mol(aa, pep, cls=None):
    m = Chem.MolFromSequence(pep); m.SetProp('_Name', aa)
    if cls: m.SetProp('Drug', cls)
    return m


def score_block(tag, L, M, train_aa, test_aa):
    """Train a base on train_aa, score test_aa. Returns {aa: Pa}. Fully checkpointed."""
    ck = f'{ROOT}/{tag}.json'
    if os.path.exists(ck):
        return json.load(open(ck))
    work = f'{ROOT}/{tag}'
    os.makedirs(work, exist_ok=True)
    for e in ('OLMPASSdoMNABase.exe', 'OLMPASSofMOL.exe', 'OLMPASSofSDF.exe'):
        if not os.path.exists(f'{work}/{e}'): shutil.copy(e, work)
    tab = TABS[L]
    if not (os.path.exists(f'{work}/B.done') and os.path.exists(f'{work}/B.MSAS')):
        w = Chem.SDWriter(f'{work}/train.sdf')
        for aa in train_aa:
            r = tab.loc[aa]
            w.write(mol(aa, r.pep, 'Resistant' if r.y else 'Not resistant'))
        w.close()
        open(f'{work}/t.txt', 'w', newline='').write(
            f'ClassLimit     3\r\nBaseCreate     {M}  B\r\nBaseAddNewData train.sdf  Drug\r\n'
            f'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
        for _ in range(3):
            for stale in glob.glob(f'{work}/B.*') + glob.glob(f'{work}/OLMPASS2*'):
                if not stale.endswith('.exe'): os.remove(stale)
            subprocess.run(['wine', 'OLMPASSdoMNABase.exe', 't.txt'], cwd=work, env=WINE,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
            if os.path.exists(f'{work}/B.MSAS'): break
        else:
            # Not an infrastructure failure: on small training folds MultiPASS declines to train
            # at all -- B.HST reports every activity type "Failed", 0 selected, average IAP
            # 0.0000. Record it so the caller can report "not computable" instead of crashing.
            json.dump({'__failed__': True}, open(ck, 'w'))
            print(f'  {tag}: MultiPASS declined to train on {len(train_aa)} variants '
                  f'({int(TABS[L].loc[list(train_aa)].y.sum())} resistant)', flush=True)
            return {'__failed__': True}
        open(f'{work}/B.done', 'w').close()

    part = f'{work}/partial.tsv'
    done = {}
    if os.path.exists(part):
        for line in open(part):
            p = line.rstrip('\n').split('\t')
            if len(p) == 2: done[p[0]] = float(p[1])
    todo = [a for a in test_aa if a not in done]

    def read_pa(op):
        """Pa for the 'Resistant' activity, 0.0 when the file is missing or has no such row."""
        if not os.path.exists(op): return 0.0
        for line in open(op):
            p = line.rstrip('\n').split('\t')
            if len(p) == 3 and p[2].strip() == 'Resistant':
                return float(p[0])
        return 0.0

    def batch(chunk, start):
        """Score a chunk of molecules with ONE OLMPASSofSDF call.

        OLMPASSofSDF.EXE <in.sdf> /B<base> [/C] writes one result file per substance, named by
        the SDF field given with /F (default `FileName`), so an index-based name is written into
        each record and mapped back afterwards. This is ~8x faster than launching
        OLMPASSofMOL.exe once per molecule -- the cost was wine process startup, not scoring --
        and it was verified to give byte-identical Pa values on 40 SMO molecules (max |diff| 0).
        """
        names = {}
        w = Chem.SDWriter(f'{work}/batch.sdf')
        for j, aa in enumerate(chunk):
            m = mol(aa, tab.loc[aa].pep)
            fn = f'm{start + j:06d}.txt'
            m.SetProp('FileName', fn)
            names[aa] = fn
            w.write(m)
        w.close()
        subprocess.run(['wine', 'OLMPASSofSDF.exe', 'batch.sdf', '/BB.MSAS', '/C'], cwd=work,
                       env=WINE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=1800)
        out = {}
        for aa, fn in names.items():
            op = f'{work}/{fn}'
            out[aa] = read_pa(op)
            os.path.exists(op) and os.remove(op)
        os.path.exists(f'{work}/batch.sdf') and os.remove(f'{work}/batch.sdf')
        return out

    if todo:
        CH = int(os.environ.get('CHUNK', 250))
        with open(part, 'a', buffering=1) as pf:
            for i in range(0, len(todo), CH):
                for aa, pa in batch(todo[i:i + CH], i).items():
                    done[aa] = pa
                    pf.write(f'{aa}\t{pa}\n')
    out = {a: done[a] for a in test_aa}
    json.dump(out, open(ck, 'w'))
    for f in glob.glob(f'{work}/B.MS*') + glob.glob(f'{work}/B.MBFA') + glob.glob(f'{work}/train.sdf'):
        os.remove(f)                      # keep the .json, drop the bulk
    return out


outer = list(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(BASE, Y, G))
AA = np.array(COMMON)

# ---------- flat: one full outer CV per configuration ----------------------------------------
flat, flat_failed = {}, []
for L, M in GRID:
    s, bad = {}, False
    for k, (tr, te) in enumerate(outer):
        r = score_block(f'flat_L{L}_M{M}_f{k}', L, M, AA[tr], AA[te])
        if r.get('__failed__'): bad = True; break
        s.update(r)
    if bad:
        flat_failed.append(f'L{L}_MNA{M}'); print(f'  flat L={L} MNA={M}: not trainable', flush=True)
        continue
    v = np.array([s[a] for a in COMMON])
    flat[(L, M)] = roc_auc_score(Y, v)
    print(f'  flat L={L} MNA={M}: AUC {flat[(L, M)]:.3f}', flush=True)

# With FLAT_ONLY the run stops after the flat grid. The grid spread alone -- how far apart the
# configurations are on the same variants -- is the cheap half of the question Table 1 raises,
# and it costs 30 blocks instead of 120. Used for the large targets, where the nested stage does
# not fit in the time available.
if os.environ.get('FLAT_ONLY'):
    fm = max(flat, key=flat.get) if flat else None
    json.dump(dict(target=NAME, flat_only=True, n_pos=int(Y.sum()), n_neg=int((1 - Y).sum()),
                   n_common=len(COMMON), grid=[f'L{L}_MNA{M}' for L, M in GRID],
                   flat_auc={f'L{L}_MNA{M}': round(v, 4) for (L, M), v in flat.items()},
                   flat_not_trainable=flat_failed,
                   flat_best_config=(f'L{fm[0]}_MNA{fm[1]}' if fm else None),
                   flat_max_auc=(round(flat[fm], 4) if fm else None),
                   flat_min_auc=(round(min(flat.values()), 4) if flat else None),
                   flat_median_auc=(round(float(np.median(list(flat.values()))), 4) if flat else None),
                   note='flat grid only; nested selection not run for this target'),
              open(f'{_SAVDR_ROOT}/multipass/nested_{NAME}_flat.json', 'w'), indent=1)
    print(f'{NAME}: flat grid only -- spread '
          f'{min(flat.values()):.3f}..{max(flat.values()):.3f}', flush=True)
    raise SystemExit(0)

# ---------- nested: choose the configuration inside each outer training fold -------------------
picked, nested_scores, inner_failed = {}, {}, 0
for k, (tr, te) in enumerate(outer):
    ytr, gtr = Y[tr], G[tr]
    inner = list(StratifiedGroupKFold(3, shuffle=True, random_state=0).split(BASE.iloc[tr], ytr, gtr))
    best, bestauc = None, -1
    for L, M in GRID:
        s, bad = {}, False
        for j, (itr, ite) in enumerate(inner):
            r = score_block(f'inner_L{L}_M{M}_o{k}_i{j}', L, M, AA[tr][itr], AA[tr][ite])
            if r.get('__failed__'): bad = True; break
            s.update(r)
        if bad:
            inner_failed += 1; continue
        v = np.array([s[a] for a in AA[tr]])
        a = roc_auc_score(ytr, v)
        print(f'  outer {k} inner L={L} MNA={M}: {a:.3f}', flush=True)
        if a > bestauc: best, bestauc = (L, M), a
    if best is None:
        print(f'  outer {k}: no configuration trainable on the inner folds -- '
              f'nested selection not computable for {NAME}', flush=True)
        json.dump(dict(target=NAME, nested_computable=False, inner_failed=inner_failed,
                       n_pos=int(Y.sum()), n_neg=int((1 - Y).sum()),
                       flat_auc={f'L{L}_MNA{M}': round(v, 4) for (L, M), v in flat.items()},
                       flat_not_trainable=flat_failed,
                       note='MultiPASS declines to train on the inner training folds: B.HST '
                            'reports every activity type Failed, 0 selected, average IAP 0.0000'),
                  open(f'{_SAVDR_ROOT}/multipass/nested_{NAME}.json', 'w'), indent=1)
        raise SystemExit(0)
    picked[k] = best
    print(f'  outer {k}: picked L={best[0]} MNA={best[1]} (inner AUC {bestauc:.3f})', flush=True)
    nested_scores.update(score_block(f'flat_L{best[0]}_M{best[1]}_f{k}', best[0], best[1],
                                     AA[tr], AA[te]))

nv = np.array([nested_scores[a] for a in COMMON])
nested_auc = roc_auc_score(Y, nv)
fm = max(flat, key=flat.get)
out = dict(target=NAME, n_pos=int(Y.sum()), n_neg=int((1 - Y).sum()), n_common=len(COMMON),
           grid=[f'L{L}_MNA{M}' for L, M in GRID],
           flat_auc={f'L{L}_MNA{M}': round(v, 4) for (L, M), v in flat.items()},
           flat_best_config=f'L{fm[0]}_MNA{fm[1]}', flat_max_auc=round(flat[fm], 4),
           flat_median_auc=round(float(np.median(list(flat.values()))), 4),
           nested_auc=round(nested_auc, 4),
           nested_picked={str(k): f'L{v[0]}_MNA{v[1]}' for k, v in picked.items()},
           optimism=round(flat[fm] - nested_auc, 4))
json.dump(out, open(f'{_SAVDR_ROOT}/multipass/nested_{NAME}.json', 'w'), indent=1)
print(json.dumps(out, indent=1), flush=True)
