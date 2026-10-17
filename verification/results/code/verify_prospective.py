"""Re-derive every number in REPORT_PROSPECTIVE.md."""
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
REP=f'{_SAVDR_ROOT}/reports/REPORT_PROSPECTIVE.md'
R=open(REP).read(); N=re.sub(r'\s+',' ',R)
ok=bad=0
def chk(n,c,d=''):
    global ok,bad
    if c: ok+=1
    else: bad+=1; print(f'  FAIL {n} {d}')
def inrep(s): chk(f'text {s!r}', s in N, '-- missing')

C=pd.read_csv(f'{_SAVDR_ROOT}/results/data/prospective_candidates.csv')
chk('82 candidates',len(C)==82,len(C))
chk('32 rejected by the wild-type gate',int((~C.wt_matches_reference).sum())==32,
    int((~C.wt_matches_reference).sum()))
ok1=C[C.wt_matches_reference]
chk('50 pass the wild-type gate',len(ok1)==50,len(ok1))
chk('24 are already training positives',
    int(ok1[ok1.in_our_missense_set].already_a_training_positive.sum())==24,
    int(ok1[ok1.in_our_missense_set].already_a_training_positive.sum()))
for g,v,ref in [('ALK','V600G','L'),('MET','G12C','L'),('BTK','L265P','P')]:
    r=C[(C.gene==g)&(C.variant==v)]
    chk(f'{g}/{v} rejected, reference {ref}',len(r)==1 and r.iloc[0].ref_residue==ref and
        not r.iloc[0].wt_matches_reference, r.ref_residue.tolist())
CU=pd.read_csv(f'{_SAVDR_ROOT}/results/data/prospective_curated.csv')
chk('15 curated',len(CU)==15,len(CU))
chk('7 judged resistance',int((CU.verdict=='resistance').sum())==7,int((CU.verdict=='resistance').sum()))
chk('8 excluded',int((CU.verdict!='resistance').sum())==8,int((CU.verdict!='resistance').sum()))
for v in ['L858R','S768I','V843I','D761N','L747S','V579A','Y167H','R132H']:
    chk(f'{v} is recorded as excluded',
        v in CU[CU.verdict!='resistance'].variant.tolist(),v)

S=pd.read_csv(f'{_SAVDR_ROOT}/results/data/prospective_scored.csv')
P=S[S.kind=='prospective']
chk('7 prospective rows',len(P)==7,len(P))
WANT={('ESR1','L536R'):(0.523,0.007,'resistant',11,12,7),
      ('EGFR','V802F'):(0.095,0.141,'not resistant',37,30,16),
      ('FLT3','Y842C'):(0.086,0.524,'not resistant',7,5,5),
      ('BTK','T474I'):(0.180,0.132,'resistant',9,3,5),
      ('ESR1','L370F'):(0.012,0.900,'not resistant',11,12,7),
      ('BTK','L528W'):(0.103,0.576,'not resistant',9,3,5),
      ('EGFR','L844V'):(0.003,0.892,'not resistant',37,30,16)}
for (g,v),(pa,pi,call,L,M,npos) in WANT.items():
    r=P[(P.gene==g)&(P.variant==v)]
    chk(f'{g}/{v} present',len(r)==1,len(r))
    if len(r)!=1: continue
    r=r.iloc[0]
    chk(f'{g}/{v} Pa',round(r.Pa,3)==pa,f'{r.Pa} != {pa}')
    chk(f'{g}/{v} Pi',round(r.Pi,3)==pi,f'{r.Pi} != {pi}')
    chk(f'{g}/{v} call',r.call==call,f'{r.call} != {call}')
    chk(f'{g}/{v} config',(r.L,r.MNA,r.n_pos)==(L,M,npos),(r.L,r.MNA,r.n_pos))
    chk(f'{g}/{v} call follows Pa>Pi',(r.Pa>r.Pi)==(call=='resistant'))
hits=P[P.call=='resistant']
chk('2 of 7 recalled',len(hits)==2,len(hits))
chk('the two hits are ESR1 L536R and BTK T474I',
    set(zip(hits.gene,hits.variant))=={('ESR1','L536R'),('BTK','T474I')},
    set(zip(hits.gene,hits.variant)))
ci=stats.binomtest(2,7).proportion_ci(method='wilson')
chk('sensitivity 0.286',round(2/7,3)==0.286)
chk('Wilson CI 0.082-0.641',round(ci.low,3)==0.082 and round(ci.high,3)==0.641,
    (round(ci.low,4),round(ci.high,4)))
ctrl=S[S.kind=='training positive']
chk('33 training positives',len(ctrl)==33,len(ctrl))
chk('33/33 recalled',int((ctrl.call=='resistant').sum())==33,int((ctrl.call=='resistant').sum()))
MEAN={'BTK':0.52,'EGFR':0.81,'ESR1':0.92,'FLT3':0.59}
for g,m in MEAN.items():
    got=ctrl[ctrl.gene==g].Pa.mean()
    chk(f'{g} control mean Pa {m}',round(got,2)==m,round(got,4))
DI=pd.read_csv(f'{_SAVDR_ROOT}/results/data/prospective_distance.csv')
EXPD={'L536R':0,'V802F':5,'Y842C':7,'T474I':7,'L370F':10,'L528W':47,'L844V':47}
for v,d in EXPD.items():
    r=DI[DI.variant==v]
    chk(f'{v} distance {d}',int(r.nearest_training_positive.iloc[0])==d,
        int(r.nearest_training_positive.iloc[0]))
chk('L536R sits on a training positive',
    bool(DI[DI.variant=='L536R'].same_position_is_a_training_positive.iloc[0]))
chk('L536H is the nearest to L536R',DI[DI.variant=='L536R'].nearest_variant.iloc[0]=='L536H')
chk('T474I is not on a training positive',
    not bool(DI[DI.variant=='T474I'].same_position_is_a_training_positive.iloc[0]))
chk('both 47-distance variants were missed',
    all(DI[DI.nearest_training_positive==47].call!='resistant'))
chk('a miss at distance 5 and a hit at 7 exist, so distance is not decisive',
    DI[DI.variant=='V802F'].call.iloc[0]!='resistant' and
    DI[DI.variant=='T474I'].call.iloc[0]=='resistant')
chk('T474I margin 0.048',round(0.180-0.132,3)==0.048)
for s_ in ['2 of 7','33 of 33','0.286','0.082–0.641','0.92 (ESR1)','0.523','0.180','0.048','v104, 19 May 2026',
           '82 candidate','32 of 82','24 of 50','8 of 15','L536H is a training positive',
           '10.1182/blood.2024027009','26 %','16 %','not really a prospective success',
           'neither necessary nor sufficient']:
    inrep(s_)
chk('report states the control distinguishes model from pipeline failure',
    'not broken and the pipeline is not misconfigured' in N)
print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
