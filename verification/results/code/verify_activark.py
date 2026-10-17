"""Independent re-derivation of every number in REPORT_ACTIVARK.md."""
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
RAW=f'{_SAVDR_ROOT}/external/primeedit'; REP=f'{_SAVDR_ROOT}/reports/REPORT_ACTIVARK.md'
R=open(REP).read()
DRUGS=['Imatinib','Nilotinib','Dasatinib','Bosutinib','Ponatinib','Asciminib']
ATP=[d for d in DRUGS if d!='Asciminib']
ok=bad=0
def chk(n,c,d=''):
    global ok,bad
    if c: ok+=1
    else: bad+=1; print(f'  FAIL {n} {d}')
def inrep(s): chk(f'text {s!r}', s in R, '-- missing from report')

ak=pd.read_csv(f'{_SAVDR_ROOT}/results/data/activark_abl1.csv')
pa=pd.read_csv(f'{_SAVDR_ROOT}/results/data/pe_scan_pa.csv')
C=pd.DataFrame({d: pd.read_csv(f'{RAW}/RS_{d}.csv').query("mut_type=='Missense'")
                .set_index('AA_var').Classification for d in DRUGS})
D=pa.merge(ak,on='AA_var').set_index('AA_var').join(C)
D['pos']=[int(''.join(c for c in i if c.isdigit())) for i in D.index]
D['atp_r']=(D[ATP]=='Resistant').any(axis=1); D['atp_s']=(D[ATP]=='Sensitive').all(axis=1)

# ---- coverage ----
chk('1848 submitted', len(D)==1848, len(D))
chk('1686 scored', int(D.RvN.notna().sum())==1686, int(D.RvN.notna().sum()))
chk('91.2 percent', abs(100*1686/1848-91.2)<0.05, f'{100*1686/1848:.2f}')
chk('162 unscored', int(D.RvN.isna().sum())==162, int(D.RvN.isna().sum()))
un=D[D.RvN.isna()]
chk('unscored span 184-211', (un.pos.min(),un.pos.max())==(184,211), (un.pos.min(),un.pos.max()))
sc=D[D.RvN.notna()]
chk('scored span 212-505', (sc.pos.min(),sc.pos.max())==(212,505), (sc.pos.min(),sc.pos.max()))
chk('D present iff RvN present', int((D.D.notna()!=D.RvN.notna()).sum())==0)
tr=set(D[D.mutType=='resistance'].index)
chk('7 training resistance variants', len(tr)==7, len(tr))
chk('training set identity',
    tr=={'D276G','E255K','F311I','F317L','G250E','T315I','Y253H'}, sorted(tr))
chk('2 increase + 1 neutral',
    int((D.mutType=='other').sum())==2 and int((D.mutType=='neutral').sum())==1,
    (int((D.mutType=='other').sum()), int((D.mutType=='neutral').sum())))
for v,p in {'T315I':0.976,'F317L':0.961,'E255K':0.917,'Y253H':0.906,'G250E':0.825,
            'F311I':0.783,'D276G':0.738}.items():
    chk(f'{v} RvN {p}', abs(D.loc[v,'RvN']-p)<5e-4, D.loc[v,'RvN'])
ours=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/ABL1_ENST00000318560.csv')
oe=ours.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); ours['AA_var']=oe[0]+oe[1]+oe[2]
o=ours.set_index('AA_var').join(D[['RvN']],how='inner')
chk('305 overlap', len(o)==305, len(o))
chk('91 of 91 positives scored',
    int(o[o.y==1].RvN.notna().sum())==91 and int(o.y.sum())==91,
    (int(o.y.sum()), int(o[o.y==1].RvN.notna().sum())))
chk('1 unscored is ATP-site resistant', int(un.atp_r.sum())==1, int(un.atp_r.sum()))

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

base=D[D.RvN.notna()&(D.mutType!='resistance')]
chk('1679 in primary scope', len(base)==1679, len(base))

WANT={'Nilotinib':(75,0.668,0.608,0.723,0.435,0.365,0.505,-0.233,-0.323,-0.144,None),
      'Imatinib' :(130,0.605,0.557,0.652,0.400,0.351,0.451,-0.205,-0.268,-0.141,None),
      'Ponatinib':(12,0.614,0.507,0.718,0.607,0.433,0.781,-0.007,-0.170,0.160,0.9255),
      'Dasatinib':(30,0.560,0.438,0.682,0.604,0.480,0.720,0.043,-0.100,0.178,0.5465),
      'Bosutinib':(12,0.660,0.495,0.810,0.752,0.606,0.877,0.091,-0.034,0.223,0.1515),
      'Asciminib':(215,0.487,0.444,0.530,0.456,0.414,0.497,-0.032,-0.089,0.025,0.2735)}
for drug,(nr,a1,l1,h1,a2,l2,h2,dl,dlo,dhi,pw) in WANT.items():
    b=base[base[drug].isin(['Resistant','Sensitive'])].copy(); b['y']=(b[drug]=='Resistant').astype(int)
    chk(f'{drug} n_res',int(b.y.sum())==nr,f'{int(b.y.sum())} != {nr}')
    g1=roc_auc_score(b.y,b.multipass_pa); g2=roc_auc_score(b.y,b.RvN)
    chk(f'{drug} ours',round(g1,3)==a1,f'{g1:.4f} != {a1}')
    chk(f'{drug} ak',round(g2,3)==a2,f'{g2:.4f} != {a2}')
    x1,y1=boot(b.y.values,b.multipass_pa.values); x2,y2=boot(b.y.values,b.RvN.values)
    chk(f'{drug} ours CI',round(x1,3)==l1 and round(y1,3)==h1,f'{x1:.4f}-{y1:.4f} != {l1}-{h1}')
    chk(f'{drug} ak CI',round(x2,3)==l2 and round(y2,3)==h2,f'{x2:.4f}-{y2:.4f} != {l2}-{h2}')
    m,lo,hi,p=paired(b.y.values,b.RvN.values,b.multipass_pa.values)
    chk(f'{drug} delta',round(m,3)==dl,f'{m:+.4f} != {dl}')
    chk(f'{drug} delta CI',round(lo,3)==dlo and round(hi,3)==dhi,f'{lo:+.4f}..{hi:+.4f} != {dlo}..{dhi}')
    if pw is not None: chk(f'{drug} p',round(p,4)==pw,f'{p:.4f} != {pw}')
    else: chk(f'{drug} p<0.001',p<0.001,f'{p:.4f}')
b=base[base.atp_r|base.atp_s].copy(); b['y']=b.atp_r.astype(int)
chk('pooled n_res 140',int(b.y.sum())==140,int(b.y.sum()))
g1=roc_auc_score(b.y,b.multipass_pa); g2=roc_auc_score(b.y,b.RvN)
chk('pooled ours 0.598',round(g1,3)==0.598,f'{g1:.4f}')
chk('pooled ak 0.417',round(g2,3)==0.417,f'{g2:.4f}')
x1,y1=boot(b.y.values,b.multipass_pa.values); x2,y2=boot(b.y.values,b.RvN.values)
chk('pooled ours CI 0.548-0.647',round(x1,3)==0.548 and round(y1,3)==0.647,f'{x1:.4f}-{y1:.4f}')
chk('pooled ak CI 0.369-0.469',round(x2,3)==0.369 and round(y2,3)==0.469,f'{x2:.4f}-{y2:.4f}')
m,lo,hi,p=paired(b.y.values,b.RvN.values,b.multipass_pa.values)
chk('pooled delta -0.180',round(m,3)==-0.180,f'{m:+.4f}')
chk('pooled delta CI',round(lo,3)==-0.246 and round(hi,3)==-0.114,f'{lo:+.4f}..{hi:+.4f}')
chk('pooled p<0.001',p<0.001,f'{p:.4f}')

# ---- caveat scopes ----
SC=[('all',base,1679,0.605,0.668,0.598,0.400,0.435,0.417,-0.180,None),
    ('minusDH',base[~((base.nda_class=='Deactivating')&(base.nda_conf=='High'))],1390,
       0.579,0.646,0.574,0.432,0.455,0.447,-0.128,None),
    ('minusD',base[base.nda_class!='Deactivating'],1021,
       0.555,0.638,0.558,0.472,0.503,0.480,-0.079,0.048),
    ('Aonly',base[base.nda_class=='Activating'],796,
       0.560,0.651,0.563,0.484,0.523,0.500,-0.064,0.1605)]
for nm,sub,n,i1,n1,p1,i2,n2,p2,dl,pw in SC:
    chk(f'{nm} n',len(sub)==n,f'{len(sub)} != {n}')
    for drug,w1,w2 in [('Imatinib',i1,i2),('Nilotinib',n1,n2)]:
        bb=sub[sub[drug].isin(['Resistant','Sensitive'])].copy(); bb['y']=(bb[drug]=='Resistant').astype(int)
        chk(f'{nm} {drug} ours',round(roc_auc_score(bb.y,bb.multipass_pa),3)==w1,
            f'{roc_auc_score(bb.y,bb.multipass_pa):.4f} != {w1}')
        chk(f'{nm} {drug} ak',round(roc_auc_score(bb.y,bb.RvN),3)==w2,
            f'{roc_auc_score(bb.y,bb.RvN):.4f} != {w2}')
    bb=sub[sub.atp_r|sub.atp_s].copy(); bb['y']=bb.atp_r.astype(int)
    chk(f'{nm} pooled ours',round(roc_auc_score(bb.y,bb.multipass_pa),3)==p1,
        f'{roc_auc_score(bb.y,bb.multipass_pa):.4f} != {p1}')
    chk(f'{nm} pooled ak',round(roc_auc_score(bb.y,bb.RvN),3)==p2,
        f'{roc_auc_score(bb.y,bb.RvN):.4f} != {p2}')
    m,lo,hi,p=paired(bb.y.values,bb.RvN.values,bb.multipass_pa.values)
    chk(f'{nm} pooled delta',round(m,3)==dl,f'{m:+.4f} != {dl}')
    if pw is not None: chk(f'{nm} pooled p',round(p,4)==pw,f'{p:.4f} != {pw}')
    else: chk(f'{nm} pooled p<0.001',p<0.001,f'{p:.4f}')

# ---- threshold sweep ----
TH={1.01:(1237,0.598,0.417),0.9:(1216,0.597,0.422),0.8:(1174,0.592,0.423),
    0.7:(1080,0.589,0.430),0.6:(997,0.585,0.445),0.5:(878,0.578,0.493),
    0.4:(725,0.576,0.511)}
# the report's columns pair n/ours with D<thr and ActivArk one column to the LEFT of the
# 0.4 entry; check each cell against its own threshold instead of trusting the layout
real={}
for thr in [1.01,0.9,0.8,0.7,0.6,0.5,0.4,0.3]:
    sub=base[base.D<thr]; bb=sub[sub.atp_r|sub.atp_s].copy(); bb['y']=bb.atp_r.astype(int)
    real[thr]=(len(bb),round(roc_auc_score(bb.y,bb.multipass_pa),3),round(roc_auc_score(bb.y,bb.RvN),3))
for thr,(n,a1,a2) in TH.items():
    chk(f'thr {thr} n',real[thr][0]==n,f'{real[thr][0]} != {n}')
    chk(f'thr {thr} ours',real[thr][1]==a1,f'{real[thr][1]} != {a1}')
chk('thr 0.3 n 557',real[0.3][0]==557,real[0.3][0])
chk('thr 0.3 ours 0.592',real[0.3][1]==0.592,real[0.3][1])
chk('ActivArk monotone to 0.5',
    [real[t][2] for t in [1.01,0.9,0.8,0.7,0.6,0.5]]==sorted([real[t][2] for t in [1.01,0.9,0.8,0.7,0.6,0.5]]),
    [real[t][2] for t in [1.01,0.9,0.8,0.7,0.6,0.5]])

# ---- correlations ----
s=base.dropna(subset=['RvN','D'])
chk('rho RvN-D +0.563',round(stats.spearmanr(s.RvN,s.D).statistic,3)==0.563,
    f'{stats.spearmanr(s.RvN,s.D).statistic:.4f}')
chk('rho RvN-ours +0.028',round(stats.spearmanr(s.RvN,s.multipass_pa).statistic,3)==0.028,
    f'{stats.spearmanr(s.RvN,s.multipass_pa).statistic:.4f}')
chk('AUC of D 0.312',round(roc_auc_score(b.y,b.D),3)==0.312,f'{roc_auc_score(b.y,b.D):.4f}')

# ---- within position ----
sub=D[D.RvN.notna()&(D.atp_r|D.atp_s)].copy(); sub['y']=sub.atp_r.astype(int)
got={}
for nm,col in [('ours','multipass_pa'),('ActivArk','RvN')]:
    c=dd=t=0; npos=0
    for p_,g in sub.groupby('pos'):
        if g.y.nunique()<2: continue
        npos+=1; p1=g[g.y==1][col].values; p0=g[g.y==0][col].values
        c+=(p1[:,None]>p0[None,:]).sum(); dd+=(p1[:,None]<p0[None,:]).sum(); t+=(p1[:,None]==p0[None,:]).sum()
    got[nm]=(round((c+0.5*t)/(c+dd+t),4),npos,int(c+dd+t),int(t),
             round(roc_auc_score(sub.y,sub.groupby('pos')[col].transform('mean')),4))
chk('ours within 0.5000',got['ours'][0]==0.5,got['ours'])
chk('ak within 0.3795',got['ActivArk'][0]==0.3795,got['ActivArk'])
chk('60 positions',got['ours'][1]==60==got['ActivArk'][1],(got['ours'][1],got['ActivArk'][1]))
chk('249 pairs',got['ours'][2]==249==got['ActivArk'][2],(got['ours'][2],got['ActivArk'][2]))
chk('ties 9 and 1',got['ours'][3]==9 and got['ActivArk'][3]==1,(got['ours'][3],got['ActivArk'][3]))
chk('posmean 0.6078 / 0.4614',got['ours'][4]==0.6078 and got['ActivArk'][4]==0.4614,
    (got['ours'][4],got['ActivArk'][4]))

for s_ in ['0.605','0.400','0.668','0.435','0.598','0.417','+0.563','0.312','0.5000','0.3795',
           '1 686','162','212–505','1 679','10.1186/s13073-025-01564-z','activark.russelllab.org',
           '2 505','441','98','0.941','D276G','T315I','0.976','GBC','37 batches of 50','882026044',
           'P00519','5 343','0.752','0.660','0.604','0.560','0.9255','0.1515','0.1605','0.5465','0.2735','+0.028','0.4614','0.6078']:
    inrep(s_)
chk('does not claim beating their published task', 'beat ActivArk at its own published task' in R)
print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
