"""Re-derive REPORT_LABEL_COVERAGE.md and check the superseded notice is in place."""
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

import pandas as pd, re, sys
REP=f'{_SAVDR_ROOT}/reports/REPORT_LABEL_COVERAGE.md'
OLD=f'{_SAVDR_ROOT}/reports/REPORT_PROSPECTIVE.md'
R=open(REP).read(); N=re.sub(r'\s+',' ',R); O=re.sub(r'\s+',' ',open(OLD).read())
ok=bad=0
def chk(n,c,d=''):
    global ok,bad
    if c: ok+=1
    else: bad+=1; print(f'  FAIL {n} {d}')
def inrep(s): chk(f'text {s!r}', s in N, '-- missing')

A=pd.read_csv(f'{_SAVDR_ROOT}/results/data/label_coverage_audit.csv')
chk('7 audited rows',len(A)==7,len(A))
YEAR={'V802F':'2009','L844V':'2014','L536R':'2019','T474I':'2024','L528W':'2024',
      'Y842C':'2024','L370F':'2026'}
for v,y in YEAR.items():
    r=A[A.variant==v]
    chk(f'{v} first reported {y}',str(r.first_reported.iloc[0])==y,r.first_reported.iloc[0])
chk('exactly one variant is from 2026',int((A.first_reported.astype(str)=='2026').sum())==1)
chk('that variant is ESR1 L370F',
    A[A.first_reported.astype(str)=='2026'].variant.iloc[0]=='L370F')
GRADE={'L370F':'sound','L536R':'sound','V802F':'sound','L528W':'indirect',
       'T474I':'excluded','Y842C':'excluded','L844V':'weak'}
for v,g in GRADE.items():
    chk(f'{v} graded {g}',A[A.variant==v].evidence_grade.iloc[0]==g,
        A[A.variant==v].evidence_grade.iloc[0])
PANEL={'L370F':'yes','L536R':'yes','V802F':'yes','L844V':'yes',
       'L528W':'class only','T474I':'no','Y842C':'no'}
for v,p in PANEL.items():
    chk(f'{v} panel status {p}',A[A.variant==v].drug_in_our_panel.iloc[0]==p,
        A[A.variant==v].drug_in_our_panel.iloc[0])
sound=A[A.evidence_grade=='sound']
chk('3 sound variants',len(sound)==3,len(sound))
chk('1 of 3 recalled',int((sound.call=='resistant').sum())==1,int((sound.call=='resistant').sum()))
chk('the sound hit is L536R',sound[sound.call=='resistant'].variant.iloc[0]=='L536R')
chk('L536R sits on a training positive',
    int(sound[sound.variant=='L536R'].nearest_training_positive.iloc[0])==0)
new=A[A.first_reported.astype(str)=='2026']
chk('0 of 1 on the genuinely new variant',int((new.call=='resistant').sum())==0)
chk('L370F Pa 0.012 Pi 0.900',
    round(float(new.Pa.iloc[0]),3)==0.012 and round(float(new.Pi.iloc[0]),3)==0.900,
    (new.Pa.iloc[0],new.Pi.iloc[0]))
chk('T474I was the only hit at an unseen position',
    A[(A.call=='resistant')&(A.nearest_training_positive>0)].variant.tolist()==['T474I'],
    A[(A.call=='resistant')&(A.nearest_training_positive>0)].variant.tolist())
# the control is unchanged
S=pd.read_csv(f'{_SAVDR_ROOT}/results/data/prospective_scored.csv'); ctrl=S[S.kind=='training positive']
chk('33/33 control',len(ctrl)==33 and int((ctrl.call=='resistant').sum())==33,
    (len(ctrl),int((ctrl.call=='resistant').sum())))
for g,m in {'BTK':0.52,'EGFR':0.81,'ESR1':0.92,'FLT3':0.59}.items():
    chk(f'{g} control mean Pa {m}',round(ctrl[ctrl.gene==g].Pa.mean(),2)==m,
        round(ctrl[ctrl.gene==g].Pa.mean(),4))
# every one of the seven is in our missense set but not a positive
C=pd.read_csv(f'{_SAVDR_ROOT}/results/data/prospective_candidates.csv')
for v in YEAR:
    r=C[C.variant==v]
    chk(f'{v} in our missense set, not a positive',
        bool(r.in_our_missense_set.iloc[0]) and not bool(r.already_a_training_positive.iloc[0]),
        (r.in_our_missense_set.iloc[0],r.already_a_training_positive.iloc[0]))
chk('the superseded notice names the replacement','REPORT_LABEL_COVERAGE.md' in O)
chk('the superseded notice states the date error','date filter constrained the search window' in O)
chk('report states it replaces the earlier one',
    'replaces `REPORT_PROSPECTIVE.md`' in R and 'whose framing was wrong' in N)
chk('report does not use the 2 of 7 figure as a result',
    'should not be used' in N)
for s_ in ['1 of 3','0 of 1','33 of 33','2009','2014','2019','2026','0.012','0.900','0.523',
           'L536H','label coverage','pirtobrutinib','compound',
           '10.1182/blood.2024027009','26 %','16 %',
           'the negative class contains known resistance variants']:
    inrep(s_)
print(f'  (part 1: {ok} passed, {bad} failed)')

# ---- addendum: the post-May-2025 search ----
print('\n== addendum ==')
A2=pd.read_csv(f'{_SAVDR_ROOT}/results/data/post_may2025_candidates.csv')
chk('44 candidates',len(A2)==44,len(A2))
chk('one candidate has zero earlier papers',int((A2.papers_before_may2025==0).sum())==1)
chk('it is MET M546L',A2[A2.papers_before_may2025==0].variant.iloc[0]=='M546L')
chk('MET M546L has no papers at all',float(A2[A2.variant=="M546L"].papers_total.iloc[0])==0)
chk('zero candidates qualify',int(A2.qualifies.sum())==0,int(A2.qualifies.sum()))
chk('L861Q has 2068 earlier papers',
    int(A2[A2.variant=='L861Q'].papers_before_may2025.iloc[0])==2068)
chk('L370F has 1 earlier paper',int(A2[A2.variant=='L370F'].papers_before_may2025.iloc[0])==1)
chk('max earlier count is 12030',int(A2.papers_before_may2025.max())==12030)
order=A2.nsmallest(7,'papers_before_may2025').variant.tolist()
chk('the seven fewest are as tabulated',
    order[:7]==['M546L','L370F','G401A','F692L','M766T','V802F','G680R'],order[:7])
for s_ in ['44 candidates','2 068','Zero candidates qualify','existence guard',
           'English word','0.012','0.900','prospective test set does not exist yet',
           '12 to 12 030']:
    inrep(s_)
print(f'\n{ok} checks passed, {bad} failed')
sys.exit(1 if bad else 0)
