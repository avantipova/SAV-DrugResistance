"""Independent re-derivation of every number in REPORT_PRIMEEDIT.md from the downloaded files."""
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

import pandas as pd, numpy as np, re, hashlib, json, sys
from scipy import stats
from sklearn.metrics import roc_auc_score

RAW=f'{_SAVDR_ROOT}/external/primeedit'; REP=f'{_SAVDR_ROOT}/reports/REPORT_PRIMEEDIT.md'
R=open(REP).read()
DRUGS=['Imatinib','Nilotinib','Dasatinib','Bosutinib','Ponatinib','Asciminib']
ATP=[d for d in DRUGS if d!='Asciminib']
ok=bad=0
def chk(name,cond,detail=''):
    global ok,bad
    if cond: ok+=1
    else: bad+=1; print(f'  FAIL {name} {detail}')
def inrep(s,name=None):
    chk(name or f'text {s!r}', s in R, '-- not in report')

# ---------- sources ----------
sizes={'Imatinib':276717,'Nilotinib':277335,'Dasatinib':277844,'Bosutinib':278192,
       'Ponatinib':277665,'Asciminib':276345}
import os
for d,s in sizes.items():
    chk(f'size {d}', os.path.getsize(f'{RAW}/RS_{d}.csv')==s, f'{os.path.getsize(f"{RAW}/RS_{d}.csv")} != {s}')
h=hashlib.sha256(open(f'{RAW}/RS_Imatinib.csv','rb').read()).hexdigest()
m=re.search(r'SHA-256 \(`([0-9a-f]+)…([0-9a-f]+)`\)',R)
chk('sha256 prefix/suffix quoted in report', bool(m) and h.startswith(m.group(1)) and h.endswith(m.group(2)),
    f'actual {h}, report {m.group(1) if m else None}…{m.group(2) if m else None}')

T={d: pd.read_csv(f'{RAW}/RS_{d}.csv') for d in DRUGS}
for d in DRUGS:
    chk(f'rows {d}', len(T[d])==2244, len(T[d]))
vc=T['Imatinib'].mut_type.value_counts()
chk('missense 1848', vc['Missense']==1848, vc['Missense'])
chk('synonymous 290', vc['Synonymous']==290, vc['Synonymous'])
chk('nonsense 106', vc['Nonsense']==106, vc['Nonsense'])
M={d: T[d][T[d].mut_type=='Missense'].set_index('AA_var') for d in DRUGS}
chk('identical variant sets', all(set(M[d].index)==set(M['Imatinib'].index) for d in DRUGS))

# ---------- numbering ----------
seq=''.join(l.strip() for l in open(f'{_SAVDR_ROOT}/inputs/fasta/ABL1_ENST00000318560.faa') if not l.startswith('>'))
chk('reference length 1130', len(seq)==1130, len(seq))
ex=pd.Series(M['Imatinib'].index).str.extract(r'^([A-Z])(\d+)([A-Z])$')
wt,pos,mt=ex[0],ex[1].astype(int),ex[2]
hits=sorted(((sum(1 for w,p in zip(wt,pos) if 1<=p+o<=len(seq) and seq[p+o-1]==w),o)
             for o in range(-30,31)),reverse=True)
chk('offset 0 matches 1848', hits[0]==(1848,0), hits[0])
chk('next best offset 147', hits[1][0]==147, hits[1])
chk('T315I present', 'T315I' in M['Imatinib'].index)
chk('T334I absent', 'T334I' not in M['Imatinib'].index)
chk('window 184-505', (pos.min(),pos.max())==(184,505), (pos.min(),pos.max()))

# ---------- our side ----------
ours=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/ABL1_ENST00000318560.csv')
oe=ours.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$')
ours['AA_var']=oe[0]+oe[1]+oe[2]; ours['pos']=oe[1].astype(int)
chk('our ABL1 rows 709', len(ours)==709, len(ours))
chk('our ABL1 positives 101', int(ours.y.sum())==101, int(ours.y.sum()))
chk('all 101 positives inside 184-505',
    int(ours[(ours.pos>=184)&(ours.pos<=505)].y.sum())==101)
def window(p,m,L=25):
    i=p-1; s=i-L//2; e=s+L
    if s<0 or e>len(seq): return ''
    q=list(seq); q[i]=m; return ''.join(q[s:e])
rep=[window(p,m) for p,m in zip(ours.pos,oe[2])]
chk('709/709 windows reproduced', (pd.Series(rep)==ours.pep).all())

pa=pd.read_csv(f'{_SAVDR_ROOT}/results/data/pe_scan_pa.csv').set_index('AA_var')
chk('1848 scored', len(pa)==1848, len(pa))
chk('scored set == scan missense set', set(pa.index)==set(M['Imatinib'].index))
j=ours.set_index('AA_var').join(pa,how='inner',rsuffix='_re')
chk('overlap 305', len(j)==305, len(j))
chk('overlap positives 91', int(j.y.sum())==91, int(j.y.sum()))
chk('Pa identical on overlap', float((j.multipass_pa-j.multipass_pa_re).abs().max())==0.0,
    float((j.multipass_pa-j.multipass_pa_re).abs().max()))

def boot(y,s,n=4000,seed=0):
    rng=np.random.default_rng(seed); y=np.asarray(y); s=np.asarray(s)
    ip,ineg=np.where(y==1)[0],np.where(y==0)[0]
    v=[roc_auc_score(y[i],s[i]) for i in
       (np.r_[rng.choice(ip,len(ip),True),rng.choice(ineg,len(ineg),True)] for _ in range(n))]
    return np.percentile(v,2.5),np.percentile(v,97.5)

# ---------- per-drug table ----------
want={'Bosutinib':(15,1546,0.708,0.563,0.840,0.053),
      'Nilotinib':(81,1517,0.686,0.630,0.740,0.087),
      'Ponatinib':(15,1554,0.671,0.568,0.768,0.051),
      'Imatinib':(137,1449,0.620,0.573,0.668,0.080),
      'Dasatinib':(33,1550,0.591,0.474,0.707,0.025),
      'Asciminib':(244,1338,0.494,0.453,0.535,-0.009)}
for d,(nr,ns,auc,lo,hi,rho) in want.items():
    t=M[d].join(pa)
    b=t[t.Classification.isin(['Resistant','Sensitive'])].copy()
    b['y']=(b.Classification=='Resistant').astype(int)
    chk(f'{d} n_res',int(b.y.sum())==nr,f'{int(b.y.sum())} != {nr}')
    chk(f'{d} n_sens',int((1-b.y).sum())==ns,f'{int((1-b.y).sum())} != {ns}')
    a=roc_auc_score(b.y,b.multipass_pa); l,hh=boot(b.y.values,b.multipass_pa.values)
    chk(f'{d} auc',round(a,3)==auc,f'{a:.4f} != {auc}')
    chk(f'{d} lo',abs(l-lo)<1e-3,f'{l:.4f} != {lo}')
    chk(f'{d} hi',abs(hh-hi)<1e-3,f'{hh:.4f} != {hi}')
    r=stats.spearmanr(t.multipass_pa,t.Resistance_score).statistic
    chk(f'{d} rho',round(r,3)==rho,f'{r:.4f} != {rho}')
    inrep(f'| {auc:.3f}'.replace('| 0','| **0')+'**' if False else f'{auc:.3f}', f'{d} auc in report')

C=pd.DataFrame({d:M[d].Classification for d in DRUGS}).loc[pa.index]
atp_r=(C[ATP]=='Resistant').any(axis=1); atp_s=(C[ATP]=='Sensitive').all(axis=1)
sel=atp_r|atp_s; y=atp_r[sel].astype(int).values; s=pa.multipass_pa[sel[sel].index].values
chk('atp n_res 147',int(y.sum())==147,int(y.sum()))
chk('atp n_sens 1201',int((1-y).sum())==1201,int((1-y).sum()))
a=roc_auc_score(y,s); l,hh=boot(y,s)
chk('atp auc 0.612',abs(a-0.612)<5e-4,f'{a:.4f}')
chk('atp CI 0.565-0.658',abs(l-0.565)<1e-3 and abs(hh-0.658)<1e-3,f'{l:.4f}-{hh:.4f}')
any6=(C=='Resistant').any(axis=1); all6=(C=='Sensitive').all(axis=1)
s6=any6|all6; y6=any6[s6].astype(int).values
chk('6drug n_res 294',int(y6.sum())==294,int(y6.sum()))
chk('6drug n_sens 1038',int((1-y6).sum())==1038,int((1-y6).sum()))
a6=roc_auc_score(y6,pa.multipass_pa[s6[s6].index].values); l6,h6=boot(y6,pa.multipass_pa[s6[s6].index].values)
chk('6drug auc 0.525',abs(a6-0.525)<5e-4,f'{a6:.4f}')
chk('6drug CI 0.487-0.562',abs(l6-0.487)<1e-3 and abs(h6-0.562)<1e-3,f'{l6:.4f}-{h6:.4f}')
chk('asciminib resistant 244',int((C.Asciminib=='Resistant').sum())==244)
chk('ponatinib resistant 15',int((C.Ponatinib=='Resistant').sum())==15)

# threshold robustness
for d,(n,a_) in {'Imatinib':(169,0.593),'Nilotinib':(107,0.658),'Asciminib':(250,0.494)}.items():
    t=M[d].join(pa); yy=(t.Resistance_score>=2.5).astype(int)
    chk(f'thr {d} n',int(yy.sum())==n,int(yy.sum()))
    chk(f'thr {d} auc',round(roc_auc_score(yy,t.multipass_pa),3)==a_,
        f'{roc_auc_score(yy,t.multipass_pa):.4f}')

# ---------- within-position ----------
sub=pa.join(atp_r.rename('r')).join(atp_s.rename('s'))
sub['pos']=[int(''.join(c for c in i if c.isdigit())) for i in sub.index]
sub=sub[sub.r|sub.s].copy(); sub['y']=sub.r.astype(int)
conc=disc=tie=0; npos=0; ws=[]
for p,g in sub.groupby('pos'):
    if g.y.nunique()<2: continue
    npos+=1; p1=g[g.y==1].multipass_pa.values; p0=g[g.y==0].multipass_pa.values
    c=(p1[:,None]>p0[None,:]).sum(); dd=(p1[:,None]<p0[None,:]).sum(); t=(p1[:,None]==p0[None,:]).sum()
    conc+=c; disc+=dd; tie+=t; ws.append((c+0.5*t)/(c+dd+t))
wauc=(conc+0.5*tie)/(conc+disc+tie)
chk('61 mixed positions',npos==61,npos)
chk('250 pairs',conc+disc+tie==250,conc+disc+tie)
chk('9 ties',tie==9,tie)
chk('within auc 0.4980',abs(wauc-0.4980)<5e-5,f'{wauc:.4f}')
chk('within mean 0.5211',abs(np.mean(ws)-0.5211)<5e-5,f'{np.mean(ws):.4f}')
rng=np.random.default_rng(0); null=[]
for _ in range(2000):
    c=dd=t=0
    for p,g in sub.groupby('pos'):
        if g.y.nunique()<2: continue
        yy=rng.permutation(g.y.values); v=g.multipass_pa.values
        a1,a0=v[yy==1],v[yy==0]
        c+=(a1[:,None]>a0[None,:]).sum(); dd+=(a1[:,None]<a0[None,:]).sum(); t+=(a1[:,None]==a0[None,:]).sum()
    null.append((c+0.5*t)/(c+dd+t))
null=np.array(null); pv=(np.sum(null>=wauc)+1)/(len(null)+1)
chk('permutation p 0.53',abs(pv-0.53)<0.005,f'{pv:.4f}')
chk('null mean 0.5007',abs(null.mean()-0.5007)<5e-5,f'{null.mean():.4f}')
chk('null 95% 0.412-0.588',abs(np.percentile(null,2.5)-0.412)<1e-3 and abs(np.percentile(null,97.5)-0.588)<1e-3,
    f'{np.percentile(null,2.5):.4f}-{np.percentile(null,97.5):.4f}')
pm=sub.groupby('pos').multipass_pa.transform('mean')
apm=roc_auc_score(sub.y,pm); lpm,hpm=boot(sub.y.values,pm.values)
chk('pos-mean auc 0.6145',abs(apm-0.6145)<5e-5,f'{apm:.4f}')
chk('pos-mean CI 0.567-0.659',abs(lpm-0.567)<1e-3 and abs(hpm-0.659)<1e-3,f'{lpm:.4f}-{hpm:.4f}')
chk('residual auc 0.5120',abs(roc_auc_score(sub.y,sub.multipass_pa-pm)-0.5120)<5e-5,
    f'{roc_auc_score(sub.y,sub.multipass_pa-pm):.4f}')
rng=np.random.default_rng(0); dd_=[]
ip,ineg=np.where(sub.y.values==1)[0],np.where(sub.y.values==0)[0]
v1,v2,yv=sub.multipass_pa.values,pm.values,sub.y.values
for _ in range(4000):
    i=np.r_[rng.choice(ip,len(ip),True),rng.choice(ineg,len(ineg),True)]
    dd_.append(roc_auc_score(yv[i],v1[i])-roc_auc_score(yv[i],v2[i]))
dd_=np.array(dd_)
chk('delta -0.0023',abs(dd_.mean()+0.0023)<5e-5,f'{dd_.mean():+.4f}')
chk('delta CI -0.0095..+0.0051',abs(np.percentile(dd_,2.5)+0.0095)<1e-3 and abs(np.percentile(dd_,97.5)-0.0051)<1e-3,
    f'{np.percentile(dd_,2.5):+.4f}..{np.percentile(dd_,97.5):+.4f}')
pdl=2*min((dd_<=0).mean(),(dd_>=0).mean())
chk('delta p 0.55',abs(pdl-0.55)<0.006,f'{pdl:.3f}')

# within-position Spearman
for d,(w,wp) in {'Imatinib':(0.019,0.41),'Nilotinib':(0.037,0.12),'Asciminib':(0.023,0.32)}.items():
    t=pa.join(M[d][['Resistance_score']])
    t['pos']=[int(''.join(c for c in i if c.isdigit())) for i in t.index]
    t['a']=t.groupby('pos').multipass_pa.rank(pct=True)-0.5
    t['b']=t.groupby('pos').Resistance_score.rank(pct=True)-0.5
    big=t[t.groupby('pos').pos.transform('size')>=3]
    r=stats.spearmanr(big.a,big.b)
    chk(f'within rho {d}',abs(r.statistic-w)<5e-4,f'{r.statistic:.4f} != {w}')
    chk(f'within rho p {d}',abs(r.pvalue-wp)<5e-3,f'{r.pvalue:.3f} != {wp}')
    chk(f'within rho n {d}',len(big)==1842,len(big))
    chk(f'within rho npos {d}',big.pos.nunique()==318,big.pos.nunique())

# ---------- label validation ----------
o=ours.set_index('AA_var').join(C,how='inner')
P=o[o.y==1]; N=o[o.y==0]
pr=(P[ATP]=='Resistant').any(axis=1); nr=(N[ATP]=='Resistant').any(axis=1)
chk('91 positives in scan',len(P)==91,len(P))
chk('214 negatives in scan',len(N)==214,len(N))
chk('53 experimentally resistant',int(pr.sum())==53,int(pr.sum()))
chk('14 negatives resistant',int(nr.sum())==14,int(nr.sum()))
chk('rate 0.582',abs(pr.mean()-0.582)<5e-4,f'{pr.mean():.4f}')
chk('rate 0.065',abs(nr.mean()-0.065)<5e-4,f'{nr.mean():.4f}')
chk('base 0.080',abs(atp_r.mean()-0.080)<5e-4,f'{atp_r.mean():.4f}')
orr,pf=stats.fisher_exact([[53,91-53],[14,214-14]])
chk('OR 19.9',abs(orr-19.9)<0.05,f'{orr:.3f}')
chk('fisher p 6.6e-22',abs(np.log10(pf)-np.log10(6.6e-22))<0.01,f'{pf:.3g}')
dsc=P[~pr]
chk('38 discordant',len(dsc)==38,len(dsc))
chk('4 asciminib-only',int((dsc.Asciminib=='Resistant').sum())==4,int((dsc.Asciminib=='Resistant').sum()))
chk('13 intermediate',int((dsc[ATP]=='Intermediate').any(axis=1).sum())==13,
    int((dsc[ATP]=='Intermediate').any(axis=1).sum()))
chk('24 fully sensitive',int((dsc[DRUGS]=='Sensitive').all(axis=1).sum())==24,
    int((dsc[DRUGS]=='Sensitive').all(axis=1).sum()))
sp=(P[ATP].isin(['Resistant','Intermediate'])).any(axis=1)
sn=(N[ATP].isin(['Resistant','Intermediate'])).any(axis=1)
sb=(C[ATP].isin(['Resistant','Intermediate'])).any(axis=1)
chk('soft 66/91',int(sp.sum())==66,int(sp.sum()))
chk('soft 0.725',abs(sp.mean()-0.725)<5e-4,f'{sp.mean():.4f}')
chk('soft 75/214',int(sn.sum())==75,int(sn.sum()))
chk('soft 0.350',abs(sn.mean()-0.350)<5e-4,f'{sn.mean():.4f}')
chk('soft base 0.350',abs(sb.mean()-0.350)<5e-4,f'{sb.mean():.4f}')

# ---------- unverified claims actually true ----------
tot=sum(int((M[d].Classification=='Resistant').sum()) for d in DRUGS)
chk('525 resistant pairs',tot==525,tot)
chk('all resistant are missense',
    sum(int((T[d].Classification=='Resistant').sum()) for d in DRUGS)==525)
chk('Intermediate max 5.92 > Resistant min 2.76',
    abs(M['Imatinib'].query("Classification=='Intermediate'").Resistance_score.max()-5.9222)<1e-3 and
    abs(M['Imatinib'].query("Classification=='Resistant'").Resistance_score.min()-2.7581)<1e-3)

# ---------- report text pins ----------
for s in ['1 848','2 802','2 892','T315I','T334I','0.498','0.6145','−0.0023','19.9',
          '6.6 × 10⁻²²','58.2 %','6.5 %','8.0 %','0.612','0.565–0.658','61','250','0.53',
          '72.5 %','35.0 %','525','361','709 of 709','0.000000','184–505','91','214','24',
          '10.1038/s41551-025-01531-4','Goosang-Yu/CML_VUS','K562','MAGeCK']:
    inrep(s)
chk('no claim of beating the comparators', 'outperform' not in R.lower())
print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
