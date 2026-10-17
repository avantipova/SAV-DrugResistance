"""Was any resistance substitution in our eleven targets FIRST reported after May 2025?

Design, after two earlier attempts went wrong:
  harvest   keyed on OUR panel drugs, so a candidate is tied from the start to a drug we model
            (Europe PMC, papers first published 2025-05-01 onwards, abstracts plus open-access
            full text); variant tokens kept only when the gene name and a resistance word both
            appear within 250 characters
  novelty   count papers mentioning (variant AND gene) with FIRST_PDATE before 2025-05-01.
            A COUNT is used, not a sort on publication date: sorting mis-ranks records with a
            missing pubYear, and it wrongly reported EGFR L861Q (2068 earlier papers) as first
            seen in 2026.
  guard     a candidate must ALSO have a non-zero unrestricted count. Without this guard the
            novelty test passes anything that does not exist at all -- which is how MET "M546L"
            got through: the context window matched the English word "met", not the gene.
"""
import pandas as pd, json

# (gene, variant, papers mentioning it before 2025-05-01, papers in total)
H = [
 ('ABL1','T315I',1583,None),('ABL1','E255K',398,None),
 ('ALK','I1171N',173,None),('ALK','G1202R',767,None),('ALK','L1196M',815,None),
 ('ALK','F1174L',637,None),('ALK','F1174V',103,None),
 ('BTK','C481S',564,None),('BTK','T474M',19,None),('BTK','T474I',94,128),
 ('BTK','T316A',39,None),('BTK','A428D',39,None),('BTK','L528W',86,124),
 ('EGFR','T790M',12030,None),('EGFR','C797S',1915,None),('EGFR','V802F',9,11),
 ('EGFR','L861Q',2068,2355),('EGFR','M766T',6,None),('EGFR','L718Q',276,None),
 ('EGFR','L844V',58,62),('EGFR','T725M',22,32),('EGFR','V769L',61,77),
 ('FLT3','F691L',258,None),('FLT3','D835Y',539,None),('FLT3','F692L',4,None),
 ('FLT3','D835V',162,None),('FLT3','Y842C',84,99),('FLT3','N676K',102,None),
 ('FLT3','K429E',13,None),
 ('KIT','D816V',1760,None),('KIT','D842V',765,None),('KIT','G680R',10,15),
 ('KIT','T670I',216,None),('KIT','R50C',12,None),('KIT','V560D',140,None),
 ('KIT','N822K',253,286),
 ('MET','D1228N',104,None),('MET','L1213V',23,None),('MET','Y1248C',29,None),
 ('MET','G401A',2,None),('MET','M546L',0,0),
 ('SMO','D473G',21,None),
 ('ESR1','L370F',1,3),('ESR1','L536R',57,75),
]
D = pd.DataFrame(H, columns=['gene','variant','papers_before_may2025','papers_total'])
D['novel_by_count'] = D.papers_before_may2025 == 0
D['exists'] = D.papers_total.fillna(1) > 0
D['qualifies'] = D.novel_by_count & D.exists
D.to_csv(f'{_SAVDR_ROOT}/results/data/post_may2025_candidates.csv', index=False)
print(f'{len(D)} candidates harvested from 2025-05-01 onwards, keyed on our panel drugs\n')
print('candidates with ZERO papers before May 2025:')
z = D[D.novel_by_count]
for _, r in z.iterrows():
    print(f'  {r.gene}/{r.variant}: before {int(r.papers_before_may2025)}, '
          f'total {r.papers_total}  -> '
          + ('QUALIFIES' if r.qualifies else 'REJECTED by the existence guard'))
print(f'\nqualifying: {int(D.qualifies.sum())}')
print(f'\nthe five with the fewest earlier papers (closest to novel):')
print(D.nsmallest(5,'papers_before_may2025')[['gene','variant','papers_before_may2025','papers_total']]
      .to_string(index=False))
json.dump(dict(n_candidates=len(D), n_qualifying=int(D.qualifies.sum()),
    rejected_by_guard=z[~z.qualifies][['gene','variant']].to_dict('records'),
    fewest=D.nsmallest(5,'papers_before_may2025')[['gene','variant','papers_before_may2025']].to_dict('records')),
    open('post_may2025_summary.json','w'), indent=1)
