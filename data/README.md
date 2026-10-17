# data/

```
data/
├── fasta/            29 × <TARGET>.faa            protein sequences (UniProt / NCBI Protein)
├── csv/              29 × <TARGET>.csv            NOT INCLUDED — COSMIC export, see below
├── mutations/        29 × <TARGET>_mutations.csv  NOT INCLUDED — COSMIC export, see below
├── prepared_data/    47 × .csv                    NOT INCLUDED — regenerated locally from
│                                                  the two above by the notebook
├── sdf/
│   ├── charged/        peptides with formal charges assigned
│   └── uncharged/      peptides without charges
├── multipass results/
│   └── merged.txt      concatenated MultiPASS run logs, parsed for AUC (IAP)
└── final models/       best model per target and drug, <TARGET>/<DRUG>/
                        (.HST .LST .MBFA .MSAR .MSAS)
```

The SDF files are the input to MultiPASS, driven by `scripts/multipass/*.ps1`;
`multipass results/` and `final models/` hold what comes back. The five model
file types are MultiPASS output: `.MSAR` is the model itself, `.MSAS` its SAS
export, `.MBFA` the descriptor data, `.HST` the full run log and `.LST` the
validation summary.

## The three folders that are not here

COSMIC requires registration and its licence does not permit redistribution, so
no COSMIC export is included in this repository or in its history. Three
folders therefore have to be rebuilt locally; `.gitignore` keeps them out once
they are. Everything else above is ours and is committed.

**`csv/<TARGET>.csv`** — the per-gene variant export from COSMIC, one file per
entry of `data/fasta/`. Columns, in order:

```
Position, CDS Mutation, AA Mutation, Legacy Mutation ID, Count, Type
```

Obtain it from <https://cancer.sanger.ac.uk/cosmic>: open the gene (or the
transcript, for the eight `<GENE>_ENST…` entries), restrict to missense
substitutions, and export the variant table as CSV. Release **v100** reproduces
the published models; **v103/v104** is what the independent validation of
Supplementary Tables S4–S6 used.

**`mutations/<TARGET>_mutations.csv`** — the per-gene drug-resistance export
from the same page, which carries the drug that each resistant substitution was
reported against. Columns, in order:

```
Gene Name, Transcript, Census Tier 1, Sample Name, Sample ID, Drug,
AA Mutation, CDS Mutation, Resistant Mutation, Primary Tissue,
Tissue Subtype 1, Tissue Subtype 2, Histology, Histology Subtype 1,
Histology Subtype 2, Pubmed Id, CGP Study, Somatic Status, Sample Type,
Zygosity, Genomic Co-ordinates
```

Only `Drug`, `AA Mutation` and `Resistant Mutation` are read by the pipeline;
the rest come with the export and are ignored.

**`prepared_data/`** — generated, not downloaded. Running
`notebooks/aas_analysis.ipynb` over the two folders above writes
`<TARGET>_prepared.csv` (one row per substitution, one column per drug, plus
the peptide window and its length) and `<TARGET>_filtered.csv` (the same after
collapsing repeated records of one substitution).

Nothing else in the pipeline reads these three folders, so `data/sdf/`,
`data/multipass results/` and `data/final models/` can be used as committed
without them. The derived tables that the manuscript's numbers are computed
from are in `verification/inputs/` and `verification/results/`, and the full
generated dataset is archived under CC BY 4.0 at
<https://zenodo.org/records/21166107>.
