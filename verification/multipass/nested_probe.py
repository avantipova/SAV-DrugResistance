"""Cheap probe: can the nested selection even be run for a target?

The expensive part of `nested_config.py` is scoring. Whether nested selection is COMPUTABLE at
all depends only on whether MultiPASS agrees to train a base on the inner training folds, which
needs no scoring. For targets with very few resistant variants it refuses -- B.HST reports every
activity type "Failed", 0 selected, average IAP 0.0000 -- and the optimism simply cannot be
estimated. This probes that in seconds instead of hours.
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
RDLogger.DisableLog('rdApp.*')

WINE = dict(os.environ,
            WINEPREFIX='/tmp/claude-0/-home-claude/2c3b11c7-e846-5338-83a9-42b57da727e2/scratchpad/wt/.wine',
            WINEARCH='win32', WINEDEBUG='-all')
TAB = f'{_SAVDR_ROOT}/inputs/variants'
out = []
for NAME in sys.argv[1:]:
    d = pd.read_csv(f'{TAB}/{NAME}.csv')
    Y, G = d.y.values, d.pos.values
    outer = list(StratifiedGroupKFold(5, shuffle=True, random_state=0).split(d, Y, G))
    tr, te = outer[0]
    inner = list(StratifiedGroupKFold(3, shuffle=True, random_state=0).split(
        d.iloc[tr], Y[tr], G[tr]))
    rec = dict(target=NAME, n=len(d), n_pos=int(Y.sum()), trials=[])
    for M in (4, 9):
        for j, (itr, _) in enumerate(inner):
            idx = tr[itr]
            work = f'probe_{NAME}_M{M}_i{j}'
            os.makedirs(work, exist_ok=True)
            shutil.copy('OLMPASSdoMNABase.exe', work)
            w = Chem.SDWriter(f'{work}/train.sdf')
            for _, r in d.iloc[idx].iterrows():
                m = Chem.MolFromSequence(r.pep)
                if m is None: continue
                m.SetProp('_Name', r.aa)
                m.SetProp('Drug', 'Resistant' if r.y else 'Not resistant')
                w.write(m)
            w.close()
            open(f'{work}/t.txt', 'w', newline='').write(
                f'ClassLimit     3\r\nBaseCreate     {M}  B\r\nBaseAddNewData train.sdf  Drug\r\n'
                f'BaseTraining\r\nBaseSave\r\nBaseSaveAsSAS\r\nBaseClose\r\n')
            subprocess.run(['wine', 'OLMPASSdoMNABase.exe', 't.txt'], cwd=work, env=WINE,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=1800)
            trained = os.path.exists(f'{work}/B.MSAS')
            npos = int(Y[idx].sum())
            rec['trials'].append(dict(mna=M, inner=j, n_train=len(idx), n_pos_train=npos,
                                      trained=trained))
            print(f'  {NAME} MNA{M} inner{j}: n={len(idx)} pos={npos} -> '
                  f'{"trained" if trained else "DECLINED"}', flush=True)
            shutil.rmtree(work, ignore_errors=True)
    # A configuration is usable only if it trains on EVERY inner fold -- one declined fold and
    # that configuration cannot be scored, so it drops out of the inner comparison. Nested
    # selection is therefore computable only if at least one MNA level trains on all three.
    rec['n_trained'] = sum(t['trained'] for t in rec['trials'])
    per_mna = {M: all(t['trained'] for t in rec['trials'] if t['mna'] == M)
               for M in sorted({t['mna'] for t in rec['trials']})}
    rec['all_inner_trained_per_mna'] = per_mna
    rec['computable'] = any(per_mna.values())
    out.append(rec)
    print(f'{NAME}: {rec["n_trained"]}/6 inner bases trained, per-MNA all-folds {per_mna} '
          f'-> nested selection {"computable" if rec["computable"] else "NOT computable"}',
          flush=True)
json.dump(out, open(f'{_SAVDR_ROOT}/multipass/nested_probe.json', 'w'), indent=1)
