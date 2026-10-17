"""ActivArk vs SAV-DrugResistance on the prime-editing scan, applying ActivArk's own
documented caveat: its resistance predictor "lacks information about activating and
deactivating information", and its help page says users "may ignore its predictions when
other predictors predict the variant to be strongly deactivating". So the comparison is
repeated after removing the variants its own A/D predictor calls Deactivating with high
confidence, and again at a series of thresholds on its deactivating probability D."""
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
from sklearn.metrics import roc_auc_score
RAW=f'{_SAVDR_ROOT}/external/primeedit'
DRUGS=['Imatinib','Nilotinib','Dasatinib','Bosutinib','Ponatinib','Asciminib']
ATP=[d for d in DRUGS if d!='Asciminib']
ak=pd.read_csv(f'{_SAVDR_ROOT}/results/data/activark_abl1.csv')
pa=pd.read_csv(f'{_SAVDR_ROOT}/results/data/pe_scan_pa.csv')
C=pd.DataFrame({d: pd.read_csv(f'{RAW}/RS_{d}.csv').query("mut_type=='Missense'")
                .set_index('AA_var').Classification for d in DRUGS})
D=pa.merge(ak,on='AA_var').set_index('AA_var').join(C)
D['pos']=[int(''.join(c for c in i if c.isdigit())) for i in D.index]
D['atp_r']=(D[ATP]=='Resistant').any(axis=1); D['atp_s']=(D[ATP]=='Sensitive').all(axis=1)
D.to_csv(f'{_SAVDR_ROOT}/results/data/ak_merged.csv')

def boot(y,s,n=4000,seed=0):
    rng=np.random.default_rng(seed); y=np.asarray(y); s=np.asarray(s)
    ip,ineg=np.where(y==1)[0],np.where(y==0)[0]
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

base = D[D.RvN.notna() & (D.mutType!='resistance')]
scopes = [
 ('ActivArk-scored, its 7 training resistance variants removed', base),
 ('also removing Deactivating (High)', base[~((base.nda_class=='Deactivating')&(base.nda_conf=='High'))]),
 ('also removing any Deactivating call',  base[base.nda_class!='Deactivating']),
 ('Activating calls only',                base[base.nda_class=='Activating']),
]
rows=[]
for nm,sub in scopes:
    print(f'\n=== {nm}: {len(sub)} variants ===')
    for drug in ['Imatinib','Nilotinib','any ATP-site TKI']:
        if drug=='any ATP-site TKI':
            b=sub[sub.atp_r|sub.atp_s].copy(); b['y']=b.atp_r.astype(int)
        else:
            b=sub[sub[drug].isin(['Resistant','Sensitive'])].copy(); b['y']=(b[drug]=='Resistant').astype(int)
        if b.y.sum()<5 or (1-b.y).sum()<5:
            print(f'  {drug:20s} too few ({int(b.y.sum())} res)'); continue
        a1=roc_auc_score(b.y,b.multipass_pa); a2=roc_auc_score(b.y,b.RvN)
        l1,h1=boot(b.y.values,b.multipass_pa.values); l2,h2=boot(b.y.values,b.RvN.values)
        m,lo,hi,p=paired(b.y.values,b.RvN.values,b.multipass_pa.values)
        rows.append(dict(scope=nm,drug=drug,n=len(b),n_res=int(b.y.sum()),
            ours=round(a1,4),ours_lo=round(l1,4),ours_hi=round(h1,4),
            ak=round(a2,4),ak_lo=round(l2,4),ak_hi=round(h2,4),
            delta=round(m,4),d_lo=round(lo,4),d_hi=round(hi,4),p=round(p,4)))
        print(f'  {drug:20s} n={len(b):4d} res={int(b.y.sum()):3d}  ours {a1:.3f} ({l1:.3f}-{h1:.3f})  '
              f'ActivArk {a2:.3f} ({l2:.3f}-{h2:.3f})  AK-ours {m:+.3f} ({lo:+.3f}..{hi:+.3f}) p={p:.3f}')
pd.DataFrame(rows).to_csv(f'{_SAVDR_ROOT}/results/data/ak_caveat.csv',index=False)

print('\n=== continuous version: AUC as the deactivating probability is capped ===')
cur=[]
for thr in [1.01,0.9,0.8,0.7,0.6,0.5,0.4,0.3]:
    sub=base[base.D<thr]; b=sub[sub.atp_r|sub.atp_s].copy(); b['y']=b.atp_r.astype(int)
    if b.y.sum()<5 or (1-b.y).sum()<5: continue
    a1=roc_auc_score(b.y,b.multipass_pa); a2=roc_auc_score(b.y,b.RvN)
    cur.append(dict(D_below=thr,n=len(b),n_res=int(b.y.sum()),ours=round(a1,4),ak=round(a2,4)))
    print(f'  D < {thr:<4} n={len(b):4d} res={int(b.y.sum()):3d}  ours {a1:.3f}  ActivArk {a2:.3f}')
pd.DataFrame(cur).to_csv(f'{_SAVDR_ROOT}/results/data/ak_threshold.csv',index=False)

# is ActivArk's RvN simply tracking its own deactivating score?
from scipy import stats
s=base.dropna(subset=['RvN','D'])
print(f'\nSpearman(RvN, D) = {stats.spearmanr(s.RvN,s.D).statistic:+.3f}')
print(f'Spearman(RvN, our Pa) = {stats.spearmanr(s.RvN,s.multipass_pa).statistic:+.3f}')
b=base[base.atp_r|base.atp_s].copy(); b["y"]=b.atp_r.astype(int)
print(f'AUC of ActivArk deactivating prob D for experimental resistance: '
      f'{roc_auc_score(b.y,b.D):.3f} (so D runs {"against" if roc_auc_score(b.y,b.D)<0.5 else "with"} resistance)')
json.dump(dict(caveat=rows,threshold=cur,
    rho_RvN_D=round(float(stats.spearmanr(s.RvN,s.D).statistic),4),
    rho_RvN_ours=round(float(stats.spearmanr(s.RvN,s.multipass_pa).statistic),4),
    auc_D=round(float(roc_auc_score(b.y,b.D)),4)),
    open(f'{_SAVDR_ROOT}/results/data/ak_caveat.json','w'),indent=1)
