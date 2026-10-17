"""Summarise the 11-protein leave-one-target-out run, and test the pattern it shows.

Baselines, cross-protein. `pos_only` is kept only as a floor: residue numbering is not
comparable between proteins, so it cannot carry information by construction. `subst_lr`
(substitution chemistry) and `full_gbm` (chemistry + normalised position + window composition)
are the meaningful cross-protein comparators.

The pattern to test: the seven protein kinases transfer, the four non-kinases do not. That is a
post-hoc grouping, so it is tested as such and reported with that caveat, not as a hypothesis
the run was designed to answer.
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

import glob, json, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy import stats

src = open('summarize.py').read()
exec(src.split('rows = []')[0])          # midrank, delong, ci, AAS/KD/VOL/CHG

L = 25
MNA = 9
TAB = f'/home/claude/tab11_L{L}'
TARGETS = json.load(open(f'{TAB}/_targets.json'))

# EC-class grouping, from the protein's family -- fixed before looking at the AUCs it explains
KINASE = {'ABL1_ENST00000318560', 'ALK', 'BTK_ENST00000308731', 'EGFR', 'FLT3', 'KIT', 'MET'}
OTHER = {'ESR1_ENST00000206249': 'nuclear receptor', 'IDH1_ENST00000415913': 'dehydrogenase',
         'NT5C2': "5'-nucleotidase", 'SMO': 'GPCR / Frizzled'}

tabs = {t: pd.read_csv(f'{TAB}/{t}.csv') for t in TARGETS}
SEQLEN = {r.target: r.seq_len for _, r in pd.read_csv(f'{TAB}/_summary.csv').iterrows()}
for t, d in tabs.items():
    d['target'] = t
    d['wt'] = [a.split('.')[1][0] for a in d.aa]
    d['mut'] = [a.split('.')[1][-1] for a in d.aa]
    d['seqlen'] = SEQLEN[t]
D = pd.concat(tabs.values(), ignore_index=True)


def feats(d, with_pos):
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
            X['win_' + a] = [p.count(a) for p in d.pep]
    return pd.DataFrame(X).values


MK = {'pos_only': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.05, class_weight='balanced', random_state=0), False, True),
      'subst_lr': (lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.3, class_weight='balanced', max_iter=2000)), False, False),
      'full_gbm': (lambda: HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05, class_weight='balanced', random_state=0), True, False)}

rows = []
for held in TARGETS:
    f = f'/home/claude/loto11_L{L}_mna{MNA}_{held}/scores.csv'
    if not glob.glob(f):
        print(f'{held}: no scores, skipped'); continue
    mp = pd.read_csv(f)
    te = tabs[held].reset_index(drop=True)
    assert list(mp.aa) == list(te.aa), f'{held}: variant order differs'
    tr = D[D.target != held]
    y = te.y.values
    sc = {'MultiPASS': np.nan_to_num(mp.multipass_pa.values, nan=0.0)}
    for nm, (mk, wp, only_pos) in MK.items():
        Xtr = (tr.pos.values / tr.seqlen.values).reshape(-1, 1) if only_pos else feats(tr, wp)
        Xte = (te.pos.values / te.seqlen.values).reshape(-1, 1) if only_pos else feats(te, wp)
        m = mk(); m.fit(Xtr, tr.y.values)
        sc[nm] = m.predict_proba(Xte)[:, 1]
    for k, v in sc.items():
        lo, hi = ci(y, v)
        line = dict(held_out=held, family='kinase' if held in KINASE else OTHER[held],
                    model=k, n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                    auc=round(roc_auc_score(y, v), 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                    auprc=round(average_precision_score(y, v), 3))
        if k != 'MultiPASS':
            _, p = delong(y, sc['MultiPASS'], v)
            line['delong_p_vs_multipass'] = float(f'{p:.3g}')
        rows.append(line)
    print(f'{held:24s} MultiPASS {sc["MultiPASS"] is not None and roc_auc_score(y, sc["MultiPASS"]):.3f}',
          flush=True)

df = pd.DataFrame(rows)
mask = df.delong_p_vs_multipass.notna()
p = df.loc[mask, 'delong_p_vs_multipass'].values
order = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
for rank, idx in enumerate(order):
    run = max(run, (m - rank) * p[idx]); adj[idx] = min(run, 1.0)
df.loc[mask, 'delong_p_holm'] = np.round(adj, 4)
df.to_csv(f'{_SAVDR_ROOT}/multipass/loto11.csv', index=False)

mp = df[df.model == 'MultiPASS'].sort_values('auc')
pd.set_option('display.width', 200)
print('\n' + '=' * 96)
print(f'LEAVE-ONE-TARGET-OUT, 11 proteins, common window L={L}, MNA {MNA}')
print('=' * 96)
print(mp[['held_out', 'family', 'n_pos', 'n_neg', 'auc', 'ci_lo', 'ci_hi', 'auprc']].to_string(index=False))
print(f'\nmean over 11 targets: {mp.auc.mean():.3f}')
print(f'CI includes 0.5 in {int(((mp.ci_lo <= 0.5) & (mp.ci_hi >= 0.5)).sum())} of {len(mp)}')

print('\nbaselines (cross-protein), mean AUC:')
for nm in ['pos_only', 'subst_lr', 'full_gbm']:
    g = df[df.model == nm]
    ahead = int((mp.set_index('held_out').auc.loc[g.held_out].values < g.auc.values).sum())
    print(f'  {nm:10s} {g.auc.mean():.3f}   beats MultiPASS in {ahead}/{len(g)}   '
          f'Holm-significant {int((g.delong_p_holm < 0.05).sum())}/{len(g)}')

print('\n' + '=' * 96)
print('POST-HOC GROUPING: protein kinases vs the rest')
print('=' * 96)
k = mp[mp.family == 'kinase']
o = mp[mp.family != 'kinase']
print('  kinases    : ' + ', '.join(f'{r.held_out.split("_")[0]} {r.auc:.3f}' for _, r in k.iterrows()))
print('  non-kinases: ' + ', '.join(f'{r.held_out.split("_")[0]} {r.auc:.3f} ({r.family})' for _, r in o.iterrows()))
print(f'  mean {k.auc.mean():.3f} vs {o.auc.mean():.3f}   difference {k.auc.mean()-o.auc.mean():+.3f}')
u = stats.mannwhitneyu(k.auc, o.auc, alternative='greater')
print(f'  Mann-Whitney U (one-sided, kinases > rest): U={u.statistic:.1f}, p={u.pvalue:.4f}')
print(f'  n: {len(k)} kinases, {len(o)} non-kinases -- this grouping was chosen AFTER seeing the')
print('  AUCs, so the p-value is descriptive, not a test of a pre-registered hypothesis.')
print('\n  confound to state: the two highest cells, FLT3 0.970 and MET 0.964, have 6 and 8')
print('  positives. Spearman rho between n_pos and AUC across the 11 targets: '
      f'{stats.spearmanr(mp.n_pos, mp.auc).statistic:+.2f}')
