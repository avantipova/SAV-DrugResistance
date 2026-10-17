"""Is the CIViC / COSMIC disagreement a labelling error or a drug-coverage gap?

Our labels are per (variant, drug) from the Drug column of each *_mutations.csv; CIViC's are
per (variant, therapies). For every CIViC resistance call this asks whether the drug CIViC
names is even present in our label set for that target. A call about a drug we never labelled
cannot be a labelling error on our side -- it is outside our coverage.
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

import pandas as pd, glob, os, re, csv, collections
U = f'{_SAVDR_ROOT}/external/source-data'
TARGET = {'ABL1': 'ABL1_ENST00000318560', 'ALK': 'ALK', 'EGFR': 'EGFR',
          'KIT': 'KIT', 'NT5C2': 'NT5C2', 'SMO': 'SMO'}
FUS = re.compile(r'\b(?:BCR|EML4|NPM1|LMNA|TPM3|ATIC|CLTC|STRN|KIF5B)::(ABL1|ALK)\b')
PAT = re.compile(r'\b(ABL1|ALK|EGFR|KIT|NT5C2|SMO)\s+([A-Z])(\d+)([A-Z])\b')
# antibodies vs small-molecule inhibitors, for reading the EGFR result
ANTIBODY = {'cetuximab', 'panitumumab', 'necitumumab', 'trastuzumab', 'pertuzumab'}
# our COSMIC export labels NT5C2 with the drug CLASS "purine"; CIViC names the members
CLASS_MEMBERS = {'purine': {'mercaptopurine', 'thioguanine', '6-mercaptopurine', '6-thioguanine'},
                 'tyrosine kinase inhibitor - ns': set()}   # a class with no member list given

def norm(s):
    s = s.strip().lower()
    s = re.sub(r'\s+mesylate$', '', s)
    s = re.sub(r'^alk inhibitor\s+', '', s)
    s = re.sub(r'\s+hydrochloride$', '', s)
    return s

ours = {}
for f in sorted(glob.glob(f'{U}/mutations/*_mutations.csv')):
    t = os.path.basename(f).replace('_mutations.csv', '')
    d = pd.read_csv(f)
    ours[t] = {}                                   # variant -> set(drugs)
    for _, r in d.iterrows():
        ours[t].setdefault(r['AA Mutation'], set()).add(norm(str(r['Drug'])))
drugset = {t: set().union(*v.values()) if v else set() for t, v in ours.items()}

civic = collections.defaultdict(lambda: collections.defaultdict(set))
for r in csv.DictReader(open(f'{U}/civic/01-Oct-2026-AcceptedClinicalEvidenceSummaries.tsv',
                             newline='', encoding='utf-8'), delimiter='\t'):
    if r['significance'] != 'Resistance' or r['evidence_type'] != 'Predictive': continue
    if r['evidence_direction'] != 'Supports': continue
    for g, wt, p, mu in PAT.findall(FUS.sub(r'\1', r['molecular_profile'])):
        for th in re.split(r'[;,]', r['therapies'] or ''):
            if th.strip(): civic[g][f'p.{wt}{p}{mu}'].add(norm(th))

rows = []
for g, t in TARGET.items():
    for v, ths in sorted(civic[g].items(), key=lambda x: int(re.search(r'\d+', x[0]).group())):
        our_drugs = ours[t].get(v, set())
        expanded = set(drugset[t])
        for c, mem in CLASS_MEMBERS.items():
            if c in drugset[t]: expanded |= mem
        shared = ths & expanded
        rows.append(dict(gene=g, variant=v,
                         civic_drugs='|'.join(sorted(ths)),
                         civic_drugs_are_antibodies=bool(ths & ANTIBODY) and not (ths - ANTIBODY),
                         our_drugs_for_this_variant='|'.join(sorted(our_drugs)),
                         drug_in_our_label_set='|'.join(sorted(shared)) if shared else '',
                         covered=bool(shared)))
df = pd.DataFrame(rows)

# attach our binary label from the analysed variant set
sc = pd.read_csv(f'{_SAVDR_ROOT}/multipass/scores_all.csv.gz')
lab = {}
for g, t in TARGET.items():
    d = sc[sc.target == t].drop_duplicates('aa')
    lab[g] = dict(zip(d.aa, d.y))
df['our_label'] = [{1: 'positive', 0: 'negative'}.get(lab[r.gene].get(r.variant), 'absent')
                   for r in df.itertuples()]
df.to_csv(f'{_SAVDR_ROOT}/civic/therapy_match.csv', index=False)

print('drugs our labels cover, per target:')
for t, s in drugset.items(): print(f'   {t:24s} {", ".join(sorted(s))}')
print(f'\nCIViC resistance calls examined: {len(df)}')
print('\n=== the discordant ones: CIViC says resistance, our label says negative ===')
d = df[df.our_label == 'negative']
print(f'{"gene":6s} {"variant":10s} {"covered":8s} {"CIViC drug(s)":40s} shared with our label set')
for r in d.itertuples():
    print(f'{r.gene:6s} {r.variant:10s} {str(r.covered):8s} {r.civic_drugs[:40]:40s} {r.drug_in_our_label_set}')
print(f'\n  of {len(d)} discordant calls, {int(d.covered.sum())} name a drug we do label '
      f'and {int((~d.covered).sum())} name a drug outside our coverage')
print(f'  antibody-only calls among them: {int(d.civic_drugs_are_antibodies.sum())}')
print('\n=== agreement on the matched ones: our positives that CIViC also calls resistant ===')
p = df[df.our_label == 'positive']
print(f'  {len(p)} such variants; drug sets overlap for {int(p.covered.sum())} of them')
ex = p[~p.covered]
if len(ex):
    print('  positives where CIViC names only drugs outside our label set:')
    for r in ex.itertuples():
        print(f'     {r.gene:6s} {r.variant:10s} CIViC={r.civic_drugs[:50]:50s} ours={r.our_drugs_for_this_variant}')
