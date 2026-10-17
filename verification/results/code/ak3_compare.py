"""ActivArk vs SAV-DrugResistance on ALK, KIT and EGFR.

Unlike the ABL1 comparison in REPORT_ACTIVARK.md, there is no experimental ground truth here:
the only labels are our own COSMIC/CIViC-derived ones. That favours us by construction -- our
LOTO model was trained on labels of this kind, albeit never on this target -- so this is a
replication check on whether the ABL1 ordering holds, not a neutral benchmark. Stated as such.
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

import pandas as pd, numpy as np, json
from scipy import stats
from sklearn.metrics import roc_auc_score
A=pd.read_csv(f'{_SAVDR_ROOT}/results/data/activark_3targets.csv')
rows=[]; detail={}
def boot(y,s,n=4000,seed=0):
    rng=np.random.default_rng(seed); y=np.asarray(y); s=np.asarray(s)
    ip,ineg=np.where(y==1)[0],np.where(y==0)[0]
    if len(ip)<3 or len(ineg)<3: return (np.nan,np.nan)
    v=[roc_auc_score(y[i],s[i]) for i in
       (np.r_[rng.choice(ip,len(ip),True),rng.choice(ineg,len(ineg),True)] for _ in range(n))]
    return float(np.percentile(v,2.5)),float(np.percentile(v,97.5))
def paired(y,s1,s2,n=4000,seed=0):
    rng=np.random.default_rng(seed); y=np.asarray(y); d=[]
    ip,ineg=np.where(y==1)[0],np.where(y==0)[0]
    for _ in range(n):
        i=np.r_[rng.choice(ip,len(ip),True),rng.choice(ineg,len(ineg),True)]
        d.append(roc_auc_score(y[i],s1[i])-roc_auc_score(y[i],s2[i]))
    d=np.array(d); return d.mean(),np.percentile(d,2.5),np.percentile(d,97.5),2*min((d<=0).mean(),(d>=0).mean())

for gene in ['ALK','KIT','EGFR']:
    loto=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/{gene}.csv')
    ex=loto.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$')
    loto['AA_var']=ex[0]+ex[1]+ex[2]; loto['pos']=ex[1].astype(int)
    a=A[A.gene==gene]
    M=loto.merge(a[['AA_var','RvN','D','mutType']],on='AA_var')
    assert len(M)==len(loto), (gene,len(M),len(loto))
    sc=M[M.RvN.notna()]
    print(f'\n=== {gene}: {len(M)} variants, {int(M.y.sum())} positives')
    print(f'    ActivArk scored {len(sc)} ({100*len(sc)/len(M):.0f}%), '
          f'holding {int(sc.y.sum())} of our {int(M.y.sum())} positives')
    print(f'    positions scored: {sc.pos.min()}-{sc.pos.max()}; '
          f"ActivArk's own resistance training variants here: {int((sc.mutType=='resistance').sum())}")
    base=sc[sc.mutType!='resistance']
    detail[gene]=dict(n=len(M), n_pos=int(M.y.sum()), scored=len(sc),
                      scored_pos=int(sc.y.sum()), pos_min=int(sc.pos.min()),
                      pos_max=int(sc.pos.max()),
                      ak_train_res=int((sc.mutType=='resistance').sum()),
                      base_n=len(base), base_pos=int(base.y.sum()))
    if base.y.sum()<3:
        print('    too few positives left after removing ActivArk training variants; no AUC')
        continue
    a1=roc_auc_score(base.y,base.multipass_pa); a2=roc_auc_score(base.y,base.RvN)
    l1,h1=boot(base.y.values,base.multipass_pa.values); l2,h2=boot(base.y.values,base.RvN.values)
    m,lo,hi,p=paired(base.y.values,base.RvN.values,base.multipass_pa.values)
    print(f'    ours {a1:.3f} ({l1:.3f}-{h1:.3f})   ActivArk {a2:.3f} ({l2:.3f}-{h2:.3f})   '
          f'AK-ours {m:+.3f} ({lo:+.3f}..{hi:+.3f}) p={p:.4f}')
    rows.append(dict(gene=gene,n=len(base),n_pos=int(base.y.sum()),
        ours=round(a1,4),ours_lo=round(l1,4),ours_hi=round(h1,4),
        ak=round(a2,4),ak_lo=round(l2,4),ak_hi=round(h2,4),
        delta=round(m,4),d_lo=round(lo,4),d_hi=round(hi,4),p=round(p,4)))
    # ActivArk's own caveat: drop strong deactivating calls (high D)
    for thr,lab in [(0.5,'D<0.5'),(0.3,'D<0.3')]:
        f=base[base.D<thr]
        if f.y.sum()>=3 and (1-f.y).sum()>=3:
            print(f'      {lab}: n={len(f):4d} pos={int(f.y.sum()):3d}  '
                  f'ours {roc_auc_score(f.y,f.multipass_pa):.3f}  '
                  f'ActivArk {roc_auc_score(f.y,f.RvN):.3f}')
    s=base.dropna(subset=['RvN','D'])
    print(f'      Spearman(RvN, D) = {stats.spearmanr(s.RvN,s.D).statistic:+.3f}; '
          f'Spearman(RvN, our Pa) = {stats.spearmanr(s.RvN,s.multipass_pa).statistic:+.3f}')
F=pd.DataFrame(rows); F.to_csv(f'{_SAVDR_ROOT}/results/data/ak3_compare.csv',index=False)
json.dump(dict(per_target=rows, coverage=detail), open('ak3_compare.json','w'), indent=1)
