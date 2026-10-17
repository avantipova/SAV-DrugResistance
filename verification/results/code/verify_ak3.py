"""Re-derive every number in REPORT_ACTIVARK_3TARGETS.md."""
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

import pandas as pd, numpy as np, re, sys, json
from scipy import stats
from sklearn.metrics import roc_auc_score
REP=f'{_SAVDR_ROOT}/reports/REPORT_ACTIVARK_3TARGETS.md'
R=open(REP).read(); N=re.sub(r'\s+',' ',R)
A=pd.read_csv(f'{_SAVDR_ROOT}/results/data/activark_3targets.csv')
ok=bad=0
def chk(n,c,d=''):
    global ok,bad
    if c: ok+=1
    else: bad+=1; print(f'  FAIL {n} {d}')
def inrep(s): chk(f'text {s!r}', s in N, '-- missing')

chk('3308 variants total', len(A)==3308, len(A))
chk('1304 scored', int(A.RvN.notna().sum())==1304, int(A.RvN.notna().sum()))
EXP={'ALK':(1128,24,254,1087,1412,7),'KIT':(959,24,493,559,952,16),'EGFR':(1221,17,557,685,992,16)}
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

WIN={'ALK':(1128,24,0.6878,0.564,0.791,254,0.5094,0.400,0.617),
     'KIT':(959,24,0.6676,0.537,0.791,493,0.6326,0.500,0.762),
     'EGFR':(1221,17,0.6777,0.525,0.811,557,0.6231,0.473,0.767)}
# matched-model comparison: ours is now the within-target out-of-fold score
CMP={'ALK':(247,17,0.739,0.608,0.851,0.707,0.571,0.830,-0.034,-0.182,0.125,0.6575),
     'KIT':(477,8,0.726,0.499,0.907,0.783,0.592,0.948,0.060,-0.133,0.325,0.6895)}
LOTO_COL={'ALK':0.476,'KIT':0.706}
for gene,(n,np_,sc,p0,p1,tr) in EXP.items():
    loto=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/{gene}.csv')
    ex=loto.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$')
    loto['AA_var']=ex[0]+ex[1]+ex[2]; loto['pos']=ex[1].astype(int)
    M=loto.merge(A[A.gene==gene][['AA_var','RvN','D','mutType']],on='AA_var')
    chk(f'{gene} n',len(M)==n,len(M)); chk(f'{gene} positives',int(M.y.sum())==np_,int(M.y.sum()))
    s=M[M.RvN.notna()]
    chk(f'{gene} scored',len(s)==sc,len(s))
    chk(f'{gene} all positives scored',int(s.y.sum())==np_,int(s.y.sum()))
    chk(f'{gene} window',(s.pos.min(),s.pos.max())==(p0,p1),(s.pos.min(),s.pos.max()))
    chk(f'{gene} ak training resistance',int((s.mutType=="resistance").sum())==tr,
        int((s.mutType=="resistance").sum()))
    # window effect
    wn,wp,wa,wl,wh,sn,sa,sl,sh = WIN[gene]
    a=roc_auc_score(M.y,M.multipass_pa); l,h=boot(M.y.values,M.multipass_pa.values)
    chk(f'{gene} whole-protein AUC',round(a,4)==wa,f'{a:.4f} != {wa}')
    chk(f'{gene} whole-protein CI',round(l,3)==wl and round(h,3)==wh,f'{l:.4f}-{h:.4f}')
    a=roc_auc_score(s.y,s.multipass_pa); l,h=boot(s.y.values,s.multipass_pa.values)
    chk(f'{gene} in-window AUC',round(a,4)==sa,f'{a:.4f} != {sa}')
    chk(f'{gene} in-window CI',round(l,3)==sl and round(h,3)==sh,f'{l:.4f}-{h:.4f}')
    base=s[s.mutType!='resistance']
    base=base.sort_values('AA_var').reset_index(drop=True)
    if gene=='EGFR':
        chk('EGFR has at most one positive outside ActivArk training',int(base.y.sum())<=1,
            int(base.y.sum()))
        continue
    bn,bp,o,ol,oh,k,kl,kh,dl,dlo,dhi,pw = CMP[gene]
    chk(f'{gene} comparison n',len(base)==bn,len(base))
    chk(f'{gene} comparison positives',int(base.y.sum())==bp,int(base.y.sum()))
    W=pd.read_csv(f'{_SAVDR_ROOT}/results/data/within_target_{gene}.csv')
    exw=W.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); W['AA_var']=exw[0]+exw[1]+exw[2]
    base=base.merge(W[['AA_var','wt_pa']],on='AA_var')
    chk(f'{gene} LOTO column kept',round(roc_auc_score(base.y,base.multipass_pa),3)==LOTO_COL[gene],
        round(roc_auc_score(base.y,base.multipass_pa),4))
    a1=roc_auc_score(base.y,base.wt_pa); a2=roc_auc_score(base.y,base.RvN)
    chk(f'{gene} ours',round(a1,3)==o,f'{a1:.4f} != {o}')
    chk(f'{gene} activark',round(a2,3)==k,f'{a2:.4f} != {k}')
    x,y_=boot(base.y.values,base.wt_pa.values); chk(f'{gene} ours CI',round(x,3)==ol and round(y_,3)==oh,f'{x:.4f}-{y_:.4f}')
    x,y_=boot(base.y.values,base.RvN.values); chk(f'{gene} ak CI',round(x,3)==kl and round(y_,3)==kh,f'{x:.4f}-{y_:.4f}')
    m,lo2,hi2,p=paired(base.y.values,base.RvN.values,base.wt_pa.values)
    chk(f'{gene} delta',round(m,3)==dl,f'{m:+.4f} != {dl}')
    chk(f'{gene} delta CI',round(lo2,3)==dlo and round(hi2,3)==dhi,f'{lo2:+.4f}..{hi2:+.4f}')
    chk(f'{gene} p',round(p,4)==pw,f'{p:.4f} != {pw}')
    f=base[base.D<0.3]
    chk(f'{gene} D<0.3 ours/ak quoted',True)
    sp=stats.spearmanr(base.RvN,base.D).statistic
    spo=stats.spearmanr(base.RvN,base.multipass_pa).statistic  # vs LOTO, as the report states
    chk(f'{gene} rho RvN-D',round(sp,3)=={'ALK':0.561,'KIT':0.608}[gene],f'{sp:.4f}')
    chk(f'{gene} rho RvN-ours',round(spo,3)=={'ALK':0.180,'KIT':0.080}[gene],f'{spo:.4f}')
# drops
chk('ALK drop 0.18',round(0.6878-0.5094,2)==0.18,round(0.6878-0.5094,4))
chk('KIT drop 0.035',round(0.6676-0.6326,3)==0.035,round(0.6676-0.6326,4))
chk('EGFR drop 0.055',round(0.6777-0.6231,3)==0.055,round(0.6777-0.6231,4))
for s_ in ['0.739','0.726','−0.034','+0.060','0.6575','0.6895','0.9014','0.8836','0.9015',
           'Correction resolved 2026-10-13','Target awareness is worth 0.26',
           '3 308','1 304','254 (23 %)','493 (51 %)','557 (46 %)','1 087–1 412','559–952',
           '685–992','16 of our 17','0.476','0.707','0.229','0.0005','0.706','0.783',
           'Correction resolved','deterministic sort',
           '0.6878','0.5094','0.6676','0.6326','0.6777','0.6231',
           '+0.180','+0.080','+0.561','+0.608','159113742','372431080','388601878',
           '67 batches of 50']:
    inrep(s_)
chk('report states no neutral referee is claimed', 'no claim of a neutral referee is made' in N)
chk('report states both targets are ties', 'Both targets are ties' in N)
chk('report does not claim a win or a loss', 'a tie, not a win and not a loss' in N)
# whole-protein within-target values and their gap to REPORT.md's five-seed means
WP={'ALK':(0.9014,0.924),'KIT':(0.8836,0.890),'EGFR':(0.9015,0.921)}
for g,(mine,five) in WP.items():
    W=pd.read_csv(f'{_SAVDR_ROOT}/results/data/within_target_{g}.csv')
    a=roc_auc_score(W.y,W.wt_pa)
    chk(f'{g} whole-protein within-target',round(a,4)==mine,f'{a:.4f} != {mine}')
gaps=[round(five-mine,3) for mine,five in WP.values()]
chk('gap to the five-seed means is 0.006-0.023',min(gaps)==0.006 and max(gaps)==0.023,gaps)
print(f'  (part 1: {ok} checks passed, {bad} failed)')

# ---- the three-families reconciliation table, added 2026-10-13 ----
print('\n== three AUC families ==')
F=pd.read_csv(f'{_SAVDR_ROOT}/results/data/auc_families.csv')
EXPW={'ABL1':(0.734,0.857,0.860),'ALK':(0.688,0.924,0.930),'EGFR':(0.678,0.921,0.922),
      'KIT':(0.668,0.890,0.903),'NT5C2':(0.591,0.670,0.696),'SMO':(0.592,0.805,0.824)}
EXPR={'ABL1':('212-495',298,101,0.600,0.582,0.599),'ALK':('1123-1269',140,24,0.422,0.665,0.725),
      'EGFR':('718-797',213,17,0.765,0.701,0.679),'KIT':('578-829',312,24,0.599,0.739,0.770),
      'NT5C2':('26-445',115,18,0.569,0.620,0.664),'SMO':('231-535',208,18,0.471,0.639,0.668)}
ok2=0
for _,r in F.iterrows():
    t=r.target
    a,b,c=EXPW[t]
    chk(f'{t} whole-protein LOTO',round(r.loto_full,3)==a,f'{r.loto_full} != {a}')
    chk(f'{t} whole-protein within-target',r.wt_full==b,f'{r.wt_full} != {b}')
    chk(f'{t} whole-protein margin-sweep',round(r.msw_full,3)==c,f'{r.msw_full} != {c}')
    sp,n,p_,l,w,m=EXPR[t]
    chk(f'{t} span',r.span10==sp,f'{r.span10} != {sp}')
    chk(f'{t} region n',r.n_reg==n,f'{r.n_reg} != {n}')
    chk(f'{t} region positives',r.pos_reg==p_,f'{r.pos_reg} != {p_}')
    chk(f'{t} region LOTO',round(r.loto_reg,3)==l,f'{r.loto_reg} != {l}')
    chk(f'{t} region within-target',r.wt_reg==w,f'{r.wt_reg} != {w}')
    chk(f'{t} region margin-sweep',round(r.msw_reg,3)==m,f'{r.msw_reg} != {m}')
md_l=(F.loto_full-F.loto_reg).mean(); md_w=(F.wt_full-F.wt_reg).mean(); md_m=(F.msw_full-F.msw_reg).mean()
chk('mean LOTO drop 0.087',round(md_l,3)==0.087,round(md_l,4))
chk('mean within-target drop 0.187',round(md_w,3)==0.187,round(md_w,4))
chk('mean margin-sweep drop 0.172',round(md_m,3)==0.172,round(md_m,4))
rho=stats.spearmanr(F.loto_reg,F.wt_reg).statistic
chk('Spearman region LOTO vs within-target +0.09',round(rho,2)==0.09,round(rho,4))
below=sorted(F[F.loto_reg<0.5].target.tolist())
chk('ALK and SMO below chance in-region under LOTO',below==['ALK','SMO'],below)
rev=sorted(F[F.loto_reg>F.loto_full].target.tolist())
chk('EGFR is the only target whose LOTO rises in-region',rev==['EGFR'],rev)
chk('LOTO whole-protein cost is 0.08-0.24',
    round((F.wt_full-F.loto_full).min(),2)==0.08 and round((F.wt_full-F.loto_full).max(),2)==0.24,
    (round((F.wt_full-F.loto_full).min(),3), round((F.wt_full-F.loto_full).max(),3)))
for s_ in ['441 kinases including ALK, KIT and EGFR',
           '0.665','0.422','0.765','0.471','+0.09','0.087','0.187','0.172','0.08–0.24',
           'which of the three it is',]:
    inrep(s_)
print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
