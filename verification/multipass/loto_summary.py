"""Pool the leave-one-target-out results: MultiPASS, the three baselines, CIs and paired tests.

Every model is evaluated on exactly the same held-out variant table (the `full`-setting table of
each target), so all comparisons are paired on identical variants and DeLong applies.

Three questions:
  1. Does MultiPASS transfer to an unseen protein?           LOTO AUC vs 0.5
  2. Does it transfer better than the baselines?             LOTO MultiPASS vs LOTO baselines
  3. How much is lost by not having seen the protein?        LOTO MultiPASS vs within-protein
                                                             MultiPASS, same variants, paired
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

import glob, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])      # first_record, midrank, delong, ci, AAS/KD/VOL/CHG, U

TARGETS = {'ABL1_ENST00000318560': (17, 9), 'ALK': (33, 4), 'EGFR': (37, 30),
           'KIT': (23, 9), 'NT5C2': (25, 4), 'SMO': (31, 7)}
MNA_LEVELS = [9, 4]      # the cross-target base needs one MNA level; both are reported because
                         # the choice turned out to matter (KIT 0.504 -> 0.619, SMO 0.665 -> 0.807)

frames = []
for t, (L, M) in TARGETS.items():
    d = pd.read_csv(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{t}_{L}_{M}_s0/scores.csv')
    seq = first_record(f'{U}/fasta/{t}.faa')
    d['target'] = t; d['seqlen'] = len(seq); d['seq'] = seq
    frames.append(d[['target', 'aa', 'pos', 'pep', 'y', 'seqlen', 'seq']])
D = pd.concat(frames, ignore_index=True)
D['wt'] = [a.split('.')[1][0] for a in D.aa]
D['mut'] = [a.split('.')[1][-1] for a in D.aa]


def feats(d, with_pos, win=7):
    X = {}
    for a in AAS:
        X['wt_' + a] = [float(w == a) for w in d.wt]
        X['mut_' + a] = [float(m == a) for m in d.mut]
    X['d_hyd'] = [KD[m] - KD[w] for w, m in zip(d.wt, d.mut)]
    X['d_vol'] = [VOL[m] - VOL[w] for w, m in zip(d.wt, d.mut)]
    X['d_chg'] = [CHG.get(m, 0) - CHG.get(w, 0) for w, m in zip(d.wt, d.mut)]
    if with_pos:
        X['pos_norm'] = list(d.pos / d.seqlen)
        for a in AAS:
            X['win_' + a] = [s[max(0, p - 1 - win):p + win].count(a) for p, s in zip(d.pos, d.seq)]
    return pd.DataFrame(X).values


MK = {'pos_only': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.05, class_weight='balanced', random_state=0), False, True),
      'subst_lr': (lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight='balanced', max_iter=2000)), False, False),
      'full_gbm': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0), True, False)}

rows = []
for held, (L, M) in TARGETS.items():
    tr, te = D[D.target != held], D[D.target == held].reset_index(drop=True)
    y = te.y.values
    sc = {}
    for lvl in MNA_LEVELS:
        f = f'loto_{held}_mna{lvl}/scores.csv'
        if not glob.glob(f):
            print(f'{held}: no LOTO MultiPASS scores at MNA {lvl}'); continue
        mp = pd.read_csv(f)
        assert list(mp.aa) == list(te.aa), f'{held}: variant order differs at MNA {lvl}'
        sc[f'MultiPASS_loto_mna{lvl}'] = np.nan_to_num(mp.multipass_pa.values, nan=0.0)
    if not sc:
        continue
    ref = f'MultiPASS_loto_mna{MNA_LEVELS[0]}'
    for nm, (mk, wp, only_pos) in MK.items():
        Xtr = (tr.pos.values / tr.seqlen.values).reshape(-1, 1) if only_pos else feats(tr, wp)
        Xte = (te.pos.values / te.seqlen.values).reshape(-1, 1) if only_pos else feats(te, wp)
        m = mk(); m.fit(Xtr, tr.y.values)
        sc[nm + '_loto'] = m.predict_proba(Xte)[:, 1]
    # within-protein MultiPASS on the same variants: seed-averaged, from the main analysis
    per = []
    for g in sorted(glob.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_{held}_{L}_{M}_s*/scores.csv')):
        dd = pd.read_csv(g)
        per.append(pd.Series(stats.rankdata(np.nan_to_num(dd.multipass_pa.values, nan=0.0)) / len(dd),
                             index=dd.aa.values))
    sc['MultiPASS_within_protein'] = pd.concat(per, axis=1).mean(axis=1).loc[te.aa].values

    for k, v in sc.items():
        lo, hi = ci(y, v)
        line = dict(held_out=held, model=k, n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                    auc=round(roc_auc_score(y, v), 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                    auprc=round(average_precision_score(y, v), 3))
        if k != ref:
            _, p = delong(y, sc[ref], v)
            line['delong_p_vs_multipass_loto'] = float(f'{p:.3g}')
        rows.append(line)

df = pd.DataFrame(rows)
mask = df.delong_p_vs_multipass_loto.notna()
p = df.loc[mask, 'delong_p_vs_multipass_loto'].values
order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)
df.to_csv(f'{_SAVDR_ROOT}/multipass/loto.csv', index=False)

w = df.pivot_table(index='held_out', columns='model', values='auc')
lot = [f'MultiPASS_loto_mna{l}' for l in MNA_LEVELS if f'MultiPASS_loto_mna{l}' in w.columns]
cols = ['MultiPASS_within_protein'] + lot + ['pos_only_loto', 'full_gbm_loto', 'subst_lr_loto']
print('=' * 104)
print('LEAVE-ONE-TARGET-OUT: MultiPASS base trained on the other five proteins')
print('=' * 104)
print(w[cols].to_string())
print('\nmean over the six targets:')
print(w[cols].mean().round(3).to_string())

ref = lot[0]
wp = df[df.model == 'MultiPASS_within_protein'].set_index('held_out')
for k in lot:
    mpl = df[df.model == k].set_index('held_out')
    print(f'\n--- {k} ---')
    for t in w.index:
        r = mpl.loc[t]
        print(f'  {t:22s} AUC {r.auc:.3f}  CI [{r.ci_lo:.3f}, {r.ci_hi:.3f}]  '
              f'{"includes 0.5" if r.ci_lo <= 0.5 <= r.ci_hi else "excludes 0.5"}'
              f'   within-protein {wp.loc[t, "auc"]:.3f}  drop {r.auc - wp.loc[t, "auc"]:+.3f}')
    print(f'  mean {mpl.auc.mean():.3f}; CI includes 0.5 in '
          f'{int(((mpl.ci_lo <= 0.5) & (mpl.ci_hi >= 0.5)).sum())}/{len(mpl)}; '
          f'mean drop from within-protein {(mpl.auc - wp.auc).mean():+.3f}')

if len(lot) == 2:
    print(f'\n--- MNA sensitivity: {lot[0]} vs {lot[1]} (paired DeLong, same variants) ---')
    a = df[df.model == lot[1]].set_index('held_out')
    for t in w.index:
        print(f'  {t:22s} {w.loc[t, lot[0]]:.3f} -> {w.loc[t, lot[1]]:.3f}  '
              f'delta {w.loc[t, lot[1]] - w.loc[t, lot[0]]:+.3f}  '
              f'p={a.loc[t, "delong_p_vs_multipass_loto"]:.3g}  '
              f'Holm={a.loc[t, "delong_p_holm"]:.3g}')
    print(f'  mean |delta| = {(w[lot[1]] - w[lot[0]]).abs().mean():.3f}; '
          f'Holm-significant in {int((a.delong_p_holm < 0.05).sum())}/{len(a)}')

print(f'\n--- within-protein vs LOTO ({ref}), paired DeLong, Holm across {m} tests ---')
for t in w.index:
    print(f'  {t:22s} p={wp.loc[t, "delong_p_vs_multipass_loto"]:.3g}  '
          f'Holm={wp.loc[t, "delong_p_holm"]:.3g}')
print(f'  Holm-significant in {int((wp.delong_p_holm < 0.05).sum())} of {len(wp)}')

print(f'\n--- LOTO MultiPASS ({ref}) vs each baseline LOTO ---')
mpl = df[df.model == ref].set_index('held_out')
for nm in ['pos_only_loto', 'full_gbm_loto', 'subst_lr_loto']:
    g = df[df.model == nm].set_index('held_out')
    print(f'  {nm:16s} mean AUC {g.auc.mean():.3f}  '
          f'MultiPASS ahead in {int((mpl.auc.loc[g.index].values > g.auc.values).sum())}/6  '
          f'Holm-significant {int((g.delong_p_holm < 0.05).sum())}/6')
