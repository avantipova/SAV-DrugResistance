# Resistance-prediction methods the manuscript does not cite

Generated 2026-10-09. A literature sweep for variant-level drug-resistance predictors outside
Table S3, done because the Comparison section will be read as a claim of completeness. Everything
below predicts, or supplies data about, **resistance of a specific substitution to a specific
drug** — pathogenicity predictors are out of scope by design (see
`../level2/REPORT_LEVEL2.md` for why they answer the wrong question).

## The one that is in our own category

**DeepMutDTA / SimSiam-MuTF** — Hu *et al.*, *Advanced Science*, 2025,
[10.1002/advs.202514899](https://doi.org/10.1002/advs.202514899). Deep contrastive learning on
**protein sequence plus drug SMILES, structure-independent**: wild-type and mutant sequence
paired with the same drug, pre-trained on ~1.46 M IC50 points from BindingDB and BioLip over
16 202 proteins, fine-tuned per task. Outputs predicted binding affinity and a ΔAffinity
classification where one end "suggests possible resistance". Evaluated on Platinum (981 mutations
over 280 protein–drug complexes), SARS-CoV-2 Mpro, HIV-1 protease, and a **cancer set of 848
mutations over 84 targets and 126 drugs with independent test splits for ABL1, ALK and EGFR**.

This is the closest published neighbour of SAV-DrugResistance — sequence-only, multi-target,
resistance-framed — and it postdates the manuscript, so missing it is not an oversight, but a
reviewer who knows the field will raise it. It belongs in the Comparison section as the current
sequence-based state of the art, with the honest note that no web server or code repository was
found for it (searched 2026-10-09; the paper's own availability statement was not reachable in
the PMC text that could be retrieved), so it cannot be run on our data either.

## Protein–ligand affinity predictors, usable for resistance, needing a 3D complex

These answer "how much does this substitution change binding affinity", which is the same
quantity SUSPECT-ABL's ΔΔG channel reports and which we already know behaves near chance inside
the ABL1 kinase domain. All three need a co-crystal or model complex per (target, drug), so each
adds one comparison per drug pair, not a panel.

| tool | reference | server | input | note |
|---|---|---|---|---|
| **mCSM-lig** | Pires, Blundell & Ascher, *Sci Rep* **6**, 29575 (2016) | `biosig.lab.uq.edu.au/mcsm_lig/` | PDB complex + mutation | framed explicitly for "emergence of drug resistance"; same lab as SUSPECT-ABL |
| **mmCSM-Lig** | Ascher lab | `biosig.lab.uq.edu.au/mmcsm_lig/` | PDB complex + mutations | the multiple-mutation successor |
| **PremPLI** | Sun *et al.*, *Commun Biol* **4**, 1311 (2021) | `lilab.jysw.suda.edu.cn/research/PremPLI/` | PDB complex + mutation | already in Table S3; page live, but **504 on structure preparation incl. its own example** — `REPORT_STRUCTURAL_SERVERS.md` |

mCSM-lig and mmCSM-Lig are **not** in Table S3 and should be, if only to be dismissed on the same
grounds as PremPLI: they require a structure of the complex, so they cannot scale across a panel
of 18 targets and 18 drugs the way a sequence model can. The SUSPECT family's other members
(`suspect_bdq`, and the pyrazinamide/rifampicin tools) are tuberculosis drugs and irrelevant here.

## Databases and label sources, not predictors

Worth knowing about because they are where an independent evaluation set would come from, and
because a reviewer may ask why CIViC was chosen over them.

- **MdrDB** — Lin *et al.*, *Commun Chem* **6**, 123 (2023),
  `quantum.tencent.com/mdrdb/`. **100 537 samples, 240 proteins, 2 503 mutations, 440 drugs**,
  with ΔΔG per sample and a derived binary resistance call at ΔΔG > 1.36 kcal/mol (8 197 samples,
  8.2 %, called resistant). Downloadable as a core set and a full set with structures and
  features. ABL is covered; the full protein list is not in the paper, so the overlap with our
  eleven targets has to be checked on the site. **This is the most promising untried external
  benchmark in this document** — it is quantitative rather than binary, large, and not derived
  from COSMIC the way DRSP and our own labels are.
- **Platinum** — the experimental protein–ligand affinity-change database DeepMutDTA benchmarks
  on, 981 mutations over 280 complexes.
- **KinaseMD**, **OncoKB** (levels R1/R2), **Cancer Genome Interpreter** — curated resistance
  annotations. OncoKB is already named in `SAV-DR_handover.md` §5.4 as the alternative to CIViC;
  it needs an academic token.

## Cited in the manuscript but not comparable, confirmed

- **Huang *et al.*, *Brief Bioinform* 24, bbad033 (2023)** — reference 21, read on 2026-10-09:
  this is a **review of 59 tools**, not a tool. Cited correctly; nothing to run.
- **CellHit** (bulk RNA-seq → patient drug responsiveness) and **scDrugMap** (single-cell
  foundation-model benchmark) — both in Table S3, both take expression data rather than a
  mutation. They cannot be benchmarked on this task at all, and the table should say so in one
  clause rather than leaving the row to imply a comparison was possible.
- **Xie *et al.*, *Molecules* 30, 3556 (2025)** — MD plus ML for BRAF resistance variants. BRAF
  is not in our panel and no server is offered.

## Published after 2024 — the sweep Anna asked for

Searched 2026-10-10, restricted to work that postdates the manuscript. Two items matter, one of
them a great deal.

### ActivArk — a live competitor with a resistance class

**Singh, Schmenger, Gonzalez-Sanchez, Kutkina, Bremec, Diwan, Mozas, Lopez, Siebert, Sotillo &
Russell, "Discriminating activating, deactivating and resistance variants in protein kinases",
*Genome Medicine* 17 (2025)**, [10.1186/s13073-025-01564-z](https://doi.org/10.1186/s13073-025-01564-z).
Server `activark.russelllab.org`, code `github.com/russelllab/kinaseResistance` under GPLv3.

Trained on 2 505 curated kinase variants — 375 activating, 1 028 deactivating, **98 resistance**,
1 004 neutral — across **441 kinases**, on sequence, evolutionary and structural features, with a
reported mean AUC of 0.941 over the four classes. It exposes a dedicated **resistance-vs-neutral**
probability (`RvN`) alongside its activating/deactivating verdicts, takes a bare `GENE/T315I` and
needs no structure from the user.

This is the comparator the manuscript most needs and does not have. It is sequence-level like
ours, it covers **seven of our eleven targets** (ABL1, ALK, BTK, EGFR, FLT3, KIT, MET — the
kinases; ESR1, IDH1, NT5C2 and SMO are out of its scope by construction), it answers the same
question rather than a neighbouring one, and unlike KRDS, D3EGFR, PremPLI and mCSM-lig **it
works**. It is benchmarked head-to-head against SAV-DrugResistance in
`REPORT_ACTIVARK.md`.

One caveat has to travel with any comparison: their 98 resistance training variants were curated
from the literature and public databases, so an overlap with our COSMIC/CIViC positives is
likely. Their output marks a variant's curated class in a `mutType` field, which makes the
overlap measurable rather than merely suspected, and the comparison is reported both with and
without those variants.

### The methodological precedent for validating against a mutational scan

**Schwarz, Giese, Gupta, Dziubanska-Kusibab, Villalba, Cherniack, Yang, Root, Karsli-Uzunbas,
Greulich, Siegel, Kamburov, Nevedomskaya, Christ, Aldeghi & Mortier, "Protein language
model-based fitness estimates facilitate resistance mutation identification", *J Chem Inf Model*
(2026)**, [10.1021/acs.jcim.6c00768](https://doi.org/10.1021/acs.jcim.6c00768).

Free-energy perturbation affinity estimates combined with protein-language-model fitness
estimates, **validated against deep mutational scanning data** — a public ERK2 / SCH772984 set
and an internal EGFR exon-20 set. Their finding is that a fitness filter removes mutations FEP
calls affinity-decreasing which the scan shows are not resistant, because the substitution
destroys the protein instead.

It is not runnable here — FEP needs a licensed engine and per-complex setup, and one of the two
datasets is proprietary — but it is worth citing for two reasons. It is a 2026 statement from
industry that **a mutational scan is the right yardstick for a resistance predictor**, which is
precisely what `REPORT_PRIMEEDIT.md` does; and its "low fitness is not resistance" point explains
part of our own false-positive pattern, since a saturating scan calls a substitution sensitive
when the cell simply fails to grow.

### Looked at and set aside

- **DRESIS 2.0**, *Nucleic Acids Research* 54(D1):D1387 (2026), `idrblab.org/dresis/`. 2 228
  structural alterations relating resistant mutants to wild type across 61 drug targets, 27
  transporters, 62 metabolising enzymes and 45 pro-survival proteins, resistance to 671 drugs,
  downloadable without a login, mutant structures filled in with AlphaFold3. A **database**, not
  a predictor — a possible fourth label source after CIViC, DRSP and MdrDB, not a comparator.
- **"Calibrated stacked learning for source-aware screening of ABL1/BCR-ABL1 mutation-TKI
  resistance"**, *Egyptian Informatics Journal* (2026), PII S1110866526001623. On topic to the
  letter, but the publisher refuses automated access (`ROBOTS_DISALLOWED`) and the paper was not
  read; flagged so it is not mistaken for a gap that was checked.
- **"In-silico predictions of drug resistance in lung cancers with EGFR mutation"**, *PASC '25*
  (2025), [10.1145/3732775.3733588](https://doi.org/10.1145/3732775.3733588). An HPC workflow —
  structure prediction, molecular dynamics, docking, quantum chemistry — not a tool a reader
  could run, and no server offered.
- **A protein-language-model benchmark for HIV drug resistance**, *Bioinformatics* 42(5) (2026).
  Different disease, no overlap with our targets; relevant only as evidence that PLM baselines
  are now the expected comparison in resistance prediction.
- **DeepMutDTA / SimSiam-MuTF** (2025) — already covered above; still no code or server found.

## Status of every resistance-specific comparator, as of 2026-10-09

| method | in Table S3 | reachable | benchmarked against us | outcome |
|---|---|---|---|---|
| SUSPECT-ABL | yes | yes | **yes** | indistinguishable near chance; `REPORT_SUSPECT_ABL.md` |
| DRSP (Lin 2022) | no, cited as ref 22 | no server, **data published** | **yes, via its dataset** | DRSP ahead at a matched operating point; `REPORT_DRSP.md` |
| Dr. Kinase | yes | yes | no — different question | scores 3 % of a kinase domain; `REPORT_DRKINASE.md` |
| D3EGFR | yes | page up, **engine dead** | no | 2 of 59 substitutions answered; `REPORT_D3EGFR.md` |
| KRDS | yes | **no** (ERR_EMPTY_RESPONSE) | no | server gone |
| PremPLI | yes | page up, **prep step 504s** | no | fails on its own example too; `REPORT_STRUCTURAL_SERVERS.md` |
| mCSM-lig | **no** | page up, **job never completes** | no | accepted, still processing after 21 min; input prepared and verified |
| mmCSM-Lig | **no** | not tried | no | same structural requirement as mCSM-lig |
| DeepMutDTA | **no** (postdates) | no code found | no | the real sequence-based rival |
| CellHit, scDrugMap | yes | — | not possible | expression input, not mutation |
| **ActivArk** (2025) | **no** (postdates) | **yes** | **yes** | the one live sequence-level rival; `REPORT_ACTIVARK.md` |
| MdrDB | no | data via its GitHub mirror | **yes, via its dataset** | our labels supported, model adds nothing; `REPORT_MDRDB.md` |
| prime-editing scan (2025) | **no** (postdates) | **data published** | **yes** | experimental ground truth, 1 848 variants; `REPORT_PRIMEEDIT.md` |
| FEP + PLM fitness (2026) | **no** (postdates) | no, licensed FEP | no | cited as precedent for scan-based validation |

## Scope of verification

**Verified.** Each server's reachability by loading it in a browser on 2026-10-08/09 and reading
the page: D3EGFR up, PremPLI up, KRDS empty response, Dr. Kinase up, SUSPECT-ABL up. That
Huang 2023 is a review, from the publication record. DRSP's data contents, from the file itself.
MdrDB's size and licence terms, from its preprint.

**Not verified.** Whether DeepMutDTA has code or weights available somewhere not surfaced by the
searches run. MdrDB's overlap with our eleven targets — the paper names only ABL among them, and
the download was not attempted. Whether PremPLI or mCSM-lig complete jobs on a longer timescale
than the window observed on 2026-10-09, or whether mmCSM-Lig behaves differently from mCSM-lig,
which was not tried. KRDS over HTTPS or from another network: one client, one attempt over HTTP. That this sweep is complete — it was a targeted search, not a systematic
review, and the field publishes faster than a search can settle.
