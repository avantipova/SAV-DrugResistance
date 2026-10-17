"""Re-derive every number in REPORT_EGFR_ANOMALY.md."""
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

import pandas as pd, numpy as np, re, sys
from scipy import stats
from sklearn.metrics import roc_auc_score
REP=f'{_SAVDR_ROOT}/reports/REPORT_EGFR_ANOMALY.md'
R=open(REP).read(); N=re.sub(r'\s+',' ',R)
ok=bad=0
def chk(n,c,d=''):
    global ok,bad
    if c: ok+=1
    else: bad+=1; print(f'  FAIL {n} {d}')
def inrep(s): chk(f'text {s!r}', s in N, '-- missing')

MS=pd.read_csv(f'{_SAVDR_ROOT}/multipass/margin_sweep.csv'); MS['short']=MS.target.str.split('_').str[0]
sp={r.short:r.span for _,r in MS[(MS.model=='multipass_full_trained')&(MS.margin=='10')].iterrows()}
NAME={'ABL1':'ABL1_ENST00000318560','ALK':'ALK','EGFR':'EGFR','KIT':'KIT','NT5C2':'NT5C2','SMO':'SMO'}
EXP={'ABL1':('212-495',0.7340,0.6001,0.7982,0.676),'ALK':('1123-1269',0.6878,0.4222,0.7190,0.895),
     'EGFR':('718-797',0.6777,0.7655,0.6607,0.837),'KIT':('578-829',0.6676,0.5987,0.6983,0.692),
     'NT5C2':('26-445',0.5905,0.5690,0.6677,0.218),'SMO':('231-535',0.5921,0.4712,0.6779,0.585)}
print('== step 1: decomposition ==')
for s,full in NAME.items():
    L=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/{full}.csv')
    ex=L.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); L['pos']=ex[1].astype(int)
    span,aw_e,ai_e,ao_e,sh_e=EXP[s]
    chk(f'{s} span',sp[s]==span,f'{sp[s]} != {span}')
    p0,p1=[int(x) for x in span.split('-')]
    L['inreg']=(L.pos>=p0)&(L.pos<=p1)
    pos=L[L.y==1]; neg=L[L.y==0]; nin=neg[neg.inreg]; nout=neg[~neg.inreg]
    def auc_vs(n):
        yy=np.r_[np.ones(len(pos)),np.zeros(len(n))]
        return roc_auc_score(yy,np.r_[pos.multipass_pa.values,n.multipass_pa.values])
    aw=roc_auc_score(L.y,L.multipass_pa); ai=auc_vs(nin); ao=auc_vs(nout); sh=len(nout)/len(neg)
    chk(f'{s} whole',round(aw,4)==aw_e,f'{aw:.4f} != {aw_e}')
    chk(f'{s} vs in-region',round(ai,4)==ai_e,f'{ai:.4f} != {ai_e}')
    chk(f'{s} vs out-region',round(ao,4)==ao_e,f'{ao:.4f} != {ao_e}')
    chk(f'{s} share out',round(sh,3)==sh_e,f'{sh:.4f} != {sh_e}')
    chk(f'{s} decomposition identity',abs((sh*ao+(1-sh)*ai)-aw)<5e-5,
        f'{sh*ao+(1-sh)*ai:.6f} vs {aw:.6f}')
    if s=='EGFR':
        chk('EGFR is the only reversal (in > out)',ai>ao)
        u=stats.mannwhitneyu(nout.multipass_pa,nin.multipass_pa,alternative='greater').pvalue
        chk('EGFR out>in Mann-Whitney 3.3e-8',abs(np.log10(u)-np.log10(3.3e-8))<0.02,f'{u:.2e}')
        chk('EGFR mean Pa in 0.0543',round(float(nin.multipass_pa.mean()),4)==0.0543,
            round(float(nin.multipass_pa.mean()),4))
        chk('EGFR mean Pa out 0.0689',round(float(nout.multipass_pa.mean()),4)==0.0689,
            round(float(nout.multipass_pa.mean()),4))
    else:
        chk(f'{s} in-region negatives score higher than out',
            nin.multipass_pa.mean()>nout.multipass_pa.mean(),
            f'{nin.multipass_pa.mean():.4f} vs {nout.multipass_pa.mean():.4f}')

print('== step 2: the activation-loop block ==')
L=pd.read_csv(f'{_SAVDR_ROOT}/inputs/loto_scores/EGFR.csv')
ex=L.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); L['pos']=ex[1].astype(int)
L['bin']=(L.pos//50)*50
g=L.groupby('bin').multipass_pa.mean()
chk('850-899 is the top window',g.idxmax()==850,g.idxmax())
for b,v in {850:0.1034,800:0.0923,900:0.0908,700:0.0620,750:0.0558}.items():
    chk(f'window {b} mean Pa {v}',round(float(g.loc[b]),4)==v,round(float(g.loc[b]),4))
al=L[(L.pos>=830)&(L.pos<=930)&(L.y==0)]; rr=L[(L.pos>=718)&(L.pos<=797)&(L.y==0)]; po=L[L.y==1]
chk('830-930 n=172',len(al)==172,len(al)); chk('718-797 n=196',len(rr)==196,len(rr))
chk('17 positives',len(po)==17,len(po))
chk('830-930 mean 0.1002',round(float(al.multipass_pa.mean()),4)==0.1002,round(float(al.multipass_pa.mean()),4))
chk('718-797 mean 0.0543',round(float(rr.multipass_pa.mean()),4)==0.0543,round(float(rr.multipass_pa.mean()),4))
chk('positives mean 0.0879',round(float(po.multipass_pa.mean()),4)==0.0879,round(float(po.multipass_pa.mean()),4))
u=stats.mannwhitneyu(al.multipass_pa,rr.multipass_pa,alternative='greater').pvalue
chk('Mann-Whitney 1.2e-33',abs(np.log10(u)-np.log10(1.2e-33))<0.03,f'{u:.2e}')
yy=np.r_[np.ones(len(po)),np.zeros(len(al))]
a_loto=roc_auc_score(yy,np.r_[po.multipass_pa.values,al.multipass_pa.values])
chk('AUC pos vs loop 0.4345 (LOTO)',round(a_loto,4)==0.4345,round(a_loto,4))
chk('below chance',a_loto<0.5)
T=pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/EGFR.csv')
et=T.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); T['pos']=et[1].astype(int)
for p_ in [854,858,861,863,866]:
    s_=T[T.pos==p_]
    chk(f'{p_} present and all negative',len(s_)>0 and int(s_.y.sum())==0,
        (len(s_),int(s_.y.sum())))
posn=sorted(T[T.y==1].pos.unique().tolist())
chk('the 17 positives sit at the nine ATP-site/gatekeeper positions',
    posn==[718,719,747,761,766,790,792,796,797],posn)
print('== step 3: where the prior came from ==')
MED={'KIT':(24,578,829,816),'FLT3':(6,691,835,835),'ABL1_ENST00000318560':(101,212,495,328),
     'ALK':(24,1123,1269,1174),'MET':(8,1181,1249,1247),'BTK_ENST00000308731':(5,316,481,481)}
for nm,(n_,lo_,hi_,md) in MED.items():
    t=pd.read_csv(f'{_SAVDR_ROOT}/inputs/variants/{nm}.csv')
    e=t.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); t['pos']=e[1].astype(int)
    p=t[t.y==1].pos
    chk(f'{nm.split("_")[0]} positives/span/median',
        (len(p),p.min(),p.max(),int(p.median()))==(n_,lo_,hi_,md),
        (len(p),p.min(),p.max(),int(p.median())))
print('== step 4: the within-target control ==')
W=pd.read_csv(f'{_SAVDR_ROOT}/results/data/within_target_EGFR.csv')
ew=W.aa.str.extract(r'^p\.([A-Z])(\d+)([A-Z])$'); W['pos']=ew[1].astype(int)
M=W.merge(L[['aa','multipass_pa']],on='aa'); chk('merge complete',len(M)==len(W),len(M))
alw=M[(M.pos>=830)&(M.pos<=930)&(M.y==0)]; rrw=M[(M.pos>=718)&(M.pos<=797)&(M.y==0)]; pow_=M[M.y==1]
chk('within 830-930 mean 0.0591',round(float(alw.wt_pa.mean()),4)==0.0591,round(float(alw.wt_pa.mean()),4))
chk('within 718-797 mean 0.2114',round(float(rrw.wt_pa.mean()),4)==0.2114,round(float(rrw.wt_pa.mean()),4))
chk('within positives mean 0.3184',round(float(pow_.wt_pa.mean()),4)==0.3184,round(float(pow_.wt_pa.mean()),4))
yy=np.r_[np.ones(len(pow_)),np.zeros(len(alw))]
a_wt=roc_auc_score(yy,np.r_[pow_.wt_pa.values,alw.wt_pa.values])
chk('AUC pos vs loop 0.9485 (within-target)',round(a_wt,4)==0.9485,round(a_wt,4))
whole_wt=roc_auc_score(M.y,M.wt_pa); reg=M[(M.pos>=718)&(M.pos<=797)]
r_wt=roc_auc_score(reg.y,reg.wt_pa); r_lo=roc_auc_score(reg.y,reg.multipass_pa)
whole_lo=roc_auc_score(M.y,M.multipass_pa)
chk('within whole 0.9015',round(whole_wt,4)==0.9015,round(whole_wt,4))
chk('within region 0.7098',round(r_wt,4)==0.7098,round(r_wt,4))
chk('within change -0.192',round(r_wt-whole_wt,3)==-0.192,round(r_wt-whole_wt,4))
chk('LOTO change +0.088',round(r_lo-whole_lo,3)==0.088,round(r_lo-whole_lo,4))
chk('within-target reverses the ordering',
    alw.wt_pa.mean()<rrw.wt_pa.mean()<pow_.wt_pa.mean())
for s_ in ['0.6777 to 0.7655','0.1002','0.0543','0.0879','0.4345','0.9485','0.0591','0.2114',
           '0.3184','0.9015','0.7098','−0.192','+0.088','1.2 × 10⁻³³','3.3 × 10⁻⁸',
           '850–899 (0.1034)','T854, L858, L861, G863 and E866','816','835','172','196',
           'target-specific models are required']:
    inrep(s_)
chk('report states the mechanism is inferred, not read out of the model',
    'no feature attribution was done' in N)
# The region LOTO AUC is exactly 0.7654562: MultiPASS emits Pa at 3 decimals and folds.py parses
# it with float() and no rounding, so the stored per-variant scores are exact.  auc_families.csv
# keeps it at 4 dp as 0.7655, and an earlier draft rounded THAT up again to 0.766 -- a double
# rounding.  At 3 dp the value is 0.765.  Both the report and the manuscript quote 4 dp.
import pandas as _pd
from sklearn.metrics import roc_auc_score as _auc
_d=_pd.read_csv(f'{_SAVDR_ROOT}/results/data/egfr_pa_profile.csv'); _r=_d[(_d.pos>=718)&(_d.pos<=797)]
_v=_auc(_r.y,_r.multipass_pa)
chk('region LOTO recomputes to 0.7654562', round(_v,7)==0.7654562, round(_v,7))
chk('3-dp rounding is 0.765, not 0.766', round(_v,3)==0.765, round(_v,3))
chk('4-dp rounding is the 0.7655 stored in auc_families.csv', round(_v,4)==0.7655, round(_v,4))
# Pa is natively 3 dp, which is what makes the recomputation above exact rather than approximate.
import glob as _g
_tot=_over=0
for _f in _g.glob(f'{_SAVDR_ROOT}/inputs/cv_folds/fold_*/scores.csv'):
    for _x in _pd.read_csv(_f, dtype={'multipass_pa':str}).multipass_pa.dropna():
        _tot+=1
        if '.' in _x and len(_x.split('.')[1])>3: _over+=1
chk('MultiPASS Pa is exactly 3 dp in every stored score', _over==0 and _tot>20000, f'{_over}/{_tot}')
# The old 3-dp pair survives in the report only inside the errata parenthetical that explains it.
# Guard the claim itself, not the explanation: the headline sentence must carry 4 dp, and the
# stale pair must not appear anywhere outside that parenthetical.
_claim = N.split('The answer in one line')[0]
chk('headline anomaly sentence quotes 4 dp', '0.6777 to 0.7655' in _claim, _claim[-200:])
chk('stale 3-dp pair appears only inside the errata parenthetical',
    _claim.count('0.678 to 0.766')==1 and 'double-rounding error' in _claim
    and N.count('0.678 to 0.766')==1, N.count('0.678 to 0.766'))
inrep('0.7654562')

print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
