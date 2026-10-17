"""Is the `region` collapse an artifact of the +-10 aa margin?

REPORT.md states the margin is a design choice with no sensitivity analysis. This sweeps it
over 5, 10, 20, 30, 50 aa and the unrestricted protein, for the three baselines plus the
untrained `nearest_pos` rule, each refit per margin on position-grouped folds (5 seeds).

MultiPASS cannot be refit per margin without hours of new compute per margin, so two things it
*can* contribute for free are reported:

  multipass_full_trained   the existing `full`-setting out-of-fold Pa, restricted to the
                           variants inside the margin. Trained protein-wide, scored inside the
                           window -- so it answers "can a protein-wide model rank within the
                           region", which is the question a user of the tool actually faces.
  multipass_region_trained the existing `region`-setting Pa, available at margin 10 only,
                           where it is the number in the main table.

Both are genuine out-of-fold scores; neither is refit here. The comparison between them at
margin 10 shows what retraining inside the region buys.
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

import glob, re, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats
RDLogger.DisableLog('rdApp.*')

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # helpers: first_record, baselines, delong, ci, AAS/KD/...

CFG = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
       'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
MARGINS = [5, 10, 20, 30, 50, None]      # None = unrestricted protein (the `full` setting)


def extract(pep, mut, length):
    parts = mut.split('.'); pos = int(''.join(c for c in parts[1] if c.isdigit())) - 1
    if pos >= len(pep): return ''
    s = pos - length // 2; e = s + length
    if s < 0: return ''
    p = list(pep); p[pos] = parts[1][-1]; f = p[s:e]
    return ''.join(f) if len(f) == length else ''


def variant_table(name, L):
    """Rebuilt exactly as folds.py builds it, including the multi-record FASTA fix."""
    seq = first_record(f'{U}/fasta/{name}.faa')
    g = pd.read_csv(f'{U}/csv/{name}.csv')
    g = g[g['Type'] == 'Substitution - Missense'].dropna(subset=['Position'])
    res = set(pd.read_csv(f'{U}/mutations/{name}_mutations.csv')['AA Mutation'])
    rows = []
    for aa in g['AA Mutation'].dropna().unique():
        pep = extract(seq, aa, L)
        if not pep or Chem.MolFromSequence(pep) is None: continue
        rows.append(dict(aa=aa, pos=int(''.join(c for c in aa.split('.')[1] if c.isdigit())),
                         pep=pep, y=int(aa in res)))
    return pd.DataFrame(rows)


def nearest_pos_scores(d, seed):
    y = d.y.values; s = np.zeros(len(y))
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=seed).split(y.reshape(-1, 1), y, d.pos.values):
        hot = np.unique(d.pos.values[tr][y[tr] == 1])
        s[te] = 0.0 if len(hot) == 0 else -np.abs(d.pos.values[te][:, None] - hot[None, :]).min(1)
    return s


def existing_mp(name, L, MNA, region):
    """Seed-averaged rank-normalised Pa from the runs already on disk, keyed by variant name."""
    tag = f'fold_{name}_{L}_{MNA}' + ('_region' if region else '') + '_s*'
    fs = sorted(glob.glob(f'{tag}/scores.csv'))
    if not fs: return None
    per = []
    for f in fs:
        d = pd.read_csv(f)
        v = stats.rankdata(np.nan_to_num(d.multipass_pa.values, nan=0.0)) / len(d)
        per.append(pd.Series(v, index=d.aa.values))
    m = pd.concat(per, axis=1)
    if m.isna().any().any(): m = m.dropna()
    return m.mean(1)


rows = []
for name, (L, MNA) in CFG.items():
    tab = variant_table(name, L)
    hot = tab.loc[tab.y == 1, 'pos']
    lo0, hi0 = hot.min(), hot.max()
    mp_full = existing_mp(name, L, MNA, False)
    mp_reg = existing_mp(name, L, MNA, True)
    for mg in MARGINS:
        d = tab if mg is None else tab[(tab.pos >= lo0 - mg) & (tab.pos <= hi0 + mg)]
        d = d.reset_index(drop=True)
        y = d.y.values
        if y.sum() < 5 or (1 - y).sum() < 20: continue
        acc = {}
        for seed in range(5):
            sc = baselines(d, name, seed)
            sc['nearest_pos'] = nearest_pos_scores(d, seed)
            for k, v in sc.items():
                acc.setdefault(k, []).append(stats.rankdata(v) / len(v))
        avg = {k: np.mean(v, 0) for k, v in acc.items()}
        if mp_full is not None and d.aa.isin(mp_full.index).all():
            avg['multipass_full_trained'] = stats.rankdata(mp_full.loc[d.aa].values) / len(d)
        if mg == 10 and mp_reg is not None and d.aa.isin(mp_reg.index).all():
            avg['multipass_region_trained'] = stats.rankdata(mp_reg.loc[d.aa].values) / len(d)
        for k, v in avg.items():
            rows.append(dict(target=name, margin=('full' if mg is None else mg),
                             span=f'{lo0}-{hi0}', n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                             model=k, auc=round(roc_auc_score(y, v), 3)))
        print(f'{name:22s} margin={str(mg):5s} {int(y.sum()):3d}+/{int((1-y).sum()):4d}-  '
              + '  '.join(f'{k}={roc_auc_score(y, v):.3f}' for k, v in avg.items()), flush=True)

df = pd.DataFrame(rows)
df.to_csv(f'{_SAVDR_ROOT}/multipass/margin_sweep.csv', index=False)

print('\n' + '=' * 110)
print('AUC vs margin (aa) -- "full" is the unrestricted protein')
print('=' * 110)
for mdl in ['multipass_full_trained', 'multipass_region_trained', 'nearest_pos', 'pos_only',
            'full_gbm', 'subst_lr']:
    sub = df[df.model == mdl]
    if sub.empty: continue
    print(f'\n{mdl}:')
    print(sub.pivot_table(index='target', columns='margin', values='auc').to_string())

print('\n' + '=' * 110)
print('n inside each margin')
print('=' * 110)
n = df[df.model == 'pos_only'].pivot_table(index='target', columns='margin',
                                           values=['n_pos', 'n_neg'], aggfunc='first')
print(n.to_string())
