"""ActivArk vs SAV-DrugResistance on ALK, KIT, EGFR, with the model the platform serves.

Replaces the leave-one-target-out side of REPORT_ACTIVARK_3TARGETS.md with within-target
position-grouped CV predictions (within_target_cv.py), so a target-aware model is compared with
a target-aware model. Everything else is unchanged: same variants, same ActivArk scores, same
exclusion of ActivArk's own curated resistance variants.
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
rows=[]; extra={}
for gene in ['ALK','KIT','EGFR']:
    W=pd.read_csv(f'{_SAVDR_ROOT}/results/data/within_target_{gene}.csv')
    ex=W.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); W['AA_var']=ex[0]+ex[1]+ex[2]
    L=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/{gene}.csv')
    exl=L.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); L['AA_var']=exl[0]+exl[1]+exl[2]
    M=W.merge(L[['AA_var','multipass_pa']],on='AA_var').merge(
        A[A.gene==gene][['AA_var','RvN','D','mutType']],on='AA_var')
    assert len(M)==len(W), (gene,len(M),len(W))
    whole_wt=roc_auc_score(M.y,M.wt_pa); whole_lo=roc_auc_score(M.y,M.multipass_pa)
    sc=M[M.RvN.notna()]
    base=sc[sc.mutType!='resistance']
    base=base.sort_values('AA_var').reset_index(drop=True)
    extra[gene]=dict(whole_within=round(whole_wt,4), whole_loto=round(whole_lo,4),
                     scored=len(sc), base_n=len(base), base_pos=int(base.y.sum()))
    print(f'\n=== {gene}: whole protein within-target {whole_wt:.4f}, LOTO {whole_lo:.4f}')
    if base.y.sum()<3:
        print(f'    only {int(base.y.sum())} positives outside ActivArk training; no AUC')
        continue
    a_w=roc_auc_score(base.y,base.wt_pa); a_l=roc_auc_score(base.y,base.multipass_pa)
    a_k=roc_auc_score(base.y,base.RvN)
    lw,hw=boot(base.y.values,base.wt_pa.values); lk,hk=boot(base.y.values,base.RvN.values)
    m,lo,hi,p=paired(base.y.values,base.RvN.values,base.wt_pa.values)
    print(f'    n={len(base)} pos={int(base.y.sum())}')
    print(f'    ours, within-target {a_w:.3f} ({lw:.3f}-{hw:.3f})')
    print(f'    ours, LOTO          {a_l:.3f}')
    print(f'    ActivArk            {a_k:.3f} ({lk:.3f}-{hk:.3f})')
    print(f'    ActivArk - ours(within-target) = {m:+.3f} ({lo:+.3f}..{hi:+.3f})  p={p:.4f}')
    rows.append(dict(gene=gene,n=len(base),n_pos=int(base.y.sum()),
        ours_wt=round(a_w,4),ours_wt_lo=round(lw,4),ours_wt_hi=round(hw,4),
        ours_loto=round(a_l,4), ak=round(a_k,4),ak_lo=round(lk,4),ak_hi=round(hk,4),
        delta=round(m,4),d_lo=round(lo,4),d_hi=round(hi,4),p=round(p,4)))
F=pd.DataFrame(rows); F.to_csv(f'{_SAVDR_ROOT}/results/data/ak3_fair.csv',index=False)
json.dump(dict(per_target=rows, whole=extra), open('ak3_fair.json','w'), indent=1)
