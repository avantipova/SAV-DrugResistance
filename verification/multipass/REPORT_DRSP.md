# MultiPASS vs DRSP — the comparator the manuscript dismissed, run anyway

Generated 2026-10-08. Third resistance-specific comparator, after `REPORT_SUSPECT_ABL.md` and
`REPORT_DRKINASE.md`.

DRSP (Lin, Y.-F. *et al.*, *Pharmaceuticals* **15**, 136, 2022) is reference 22 of the
manuscript — the citation behind the sentence *"their practical utility is often limited by the
lack of publicly accessible web interfaces"*. That sentence is true: there is no server and no
code. But DRSP is a method, not a review, and **its labelled dataset is published**, so our model
can be scored on exactly its variants and the result placed beside its published figures. No
server needed.

**Conclusion first: DRSP comes out ahead of MultiPASS at its own operating point, and the
manuscript must not claim otherwise.** At DRSP's specificity of 0.9119 its reported sensitivity
is 0.654; ours on its variants is **0.491 (95% CI 0.370–0.611)** — the interval excludes their
value. Two things qualify that, and both belong in the paper next to it: their figure is
cross-validated on the set they fitted, with a feature set chosen by a genetic algorithm on that
same data and the best of eight runs kept by optimised MCC, so it is optimistically biased by an
amount nobody has measured; and the two figures are not on identical variants, since ours is a
302-variant subset of their 610.

## The method

A classifier of single amino acid variants as drug-resistant or not, from 45 features in four
groups: drug–SAV interaction energies from iGEMDOCK, weighted contact number, structural terms
(B-factor, DSSP) and sequence conservation (BLOSUM62, PAM250, PSSM, entropy). It needs a
**drug–protein complex** — a crystal structure where one exists, a docking model otherwise. Two
variants of the model, exterior and interior of an 8 Å shell, are combined.

Published performance, Table 1, five-fold CV on the 610-SAV training set:

| model | accuracy | sensitivity | specificity | MCC | precision | F1 |
|---|---|---|---|---|---|---|
| **DRSPce+DRSPci (crystal, combined)** | **0.8557** | **0.6541** | **0.9119** | **0.5724** | 0.6744 | **0.6641** |
| DRSPse+DRSPsi (docked, combined) | 0.8508 | 0.6391 | 0.9099 | 0.5567 | 0.6641 | 0.6513 |
| DRSPci (crystal, interior) | 0.7819 | 0.7284 | 0.8224 | 0.5536 | 0.7564 | 0.7421 |

The abstract's "86% / 0.57 / 0.66" is the first row rounded; the arithmetic was re-derived from
their own sensitivity and specificity and reproduces all four of accuracy, precision, F1 and MCC
exactly (TP 87, FP 42, FN 46, TN 435). Their independent test set carries **3 positives** and
they report only specificity for it, so it carries no information. No AUC is reported anywhere
in the paper.

## Getting the data, and the numbering trap

`Dataset.xlsx` sits inside `pharmaceuticals-15-00136-s001.zip` (167 711 bytes, sha256
`1ce4bd21…a0d6`), two sheets:

- **Training set**, 610 SAVs over 11 genes. Positives carry the drug name; negatives are
  `NA_deleterious`, i.e. damaging-but-not-resistance. 133 / 477, 21.8 % positive.
- **Testing set**, 115 SAVs over BRAF, MAP2K2 and ROS1 — none of which we model, so only the
  training set overlaps our panel. Nine of its eleven genes are ours (not MAP2K1, PDGFRA).

**DRSP does not number on the same sequences we do, and this nearly wrecked the comparison.**
The wild-type-residue check that every external join in this project goes through came back with
57 mismatches on first run. The offset was then found by brute force, as the integer that
maximises wild-type agreement:

| gene | offset (DRSP = ours + offset) | wt agreement at that offset | at offset 0 |
|---|---|---|---|
| **ABL1** | **+19** | 64/64 | **0/45** |
| **BTK** | **+34** | 31/31 | **0/12** |
| ALK, EGFR, ESR1, FLT3, KIT, MET, SMO | 0 | 100 % each | same |

+19 is ABL1 isoform b (1149 residues) against our ENST00000318560 (1130) — the same offset that
bit the dbNSFP join. +34 is canonical BTK (693) against our ENST00000308731 (659). The shift is
confirmed a third way: under it the ten canonical ABL1 resistance mutations land exactly where
they should — DRSP's M263V is M244V, T334I is **T315I**, Y272H is Y253H, E274K is E255K, M370T is
M351T, H415R is H396R — and all ten are present in our variant table. Had the offset been left at
0, every ABL1 and BTK row would have been dropped or, worse, matched to the wrong residue.

After the shift: **440 of 511** variants in the nine shared genes matched, **all 130 positives
among them**, and the 71 unmatched are all negatives — `NA_deleterious` variants absent from our
COSMIC v100 missense export. Prevalence on the matched set is therefore 29.5 % against their
21.8 %, which is why the comparison below is made at a fixed specificity rather than on accuracy.

## DRSP is not an independent label source

| | our label 0 | our label 1 |
|---|---|---|
| **DRSP label 0** | 309 | 1 |
| **DRSP label 1** | 0 | 130 |

**99.8 % agreement, Cohen κ = 0.995.** DRSP's resistance calls and ours are, to one variant, the
same calls — both trace to the same curated resistance literature. So this is *not* a CIViC-style
independent test; it is the same labels on a shared subset of variants, which makes it a clean
comparison of *methods* and no test of our labels at all. The one disagreement is **ESR1 E380Q**,
which we label resistant (endocrine therapy) and DRSP labels `NA_deleterious`; E380Q is a
well-described ESR1 activating mutation, so this looks like an omission on their side rather than
on ours.

## Result

Two score panels, because our coverage differs. Nothing is refitted in either.

- **B — within-protein CV**, the `full` setting with each target's own best configuration,
  out-of-fold under position-grouped CV, rank-averaged over the 5 fold-assignment seeds. Exists
  for ABL1, ALK, EGFR, KIT, SMO: 302 of the 440 variants.
- **A — LOTO**, leave-one-target-out at the common L = 25 / MNA 9, where the held-out protein was
  never in training. Exists for all nine genes, all 440 variants. The most conservative reading
  there is.

| panel | n | pos / neg | AUC | 95% CI | AUPRC | MCC at prevalence-matched cut | best MCC |
|---|---|---|---|---|---|---|---|
| **B — within-protein CV, 5 genes** | 302 | 108 / 194 | **0.805** | 0.755–0.849 | 0.722 | 0.409 | **0.483** |
| A — LOTO, all 9 genes | 440 | 130 / 310 | 0.608 | 0.552–0.669 | 0.384 | 0.192 | 0.212 |
| A5 — LOTO, the same 5 genes | 302 | 108 / 194 | 0.597 | 0.534–0.669 | 0.442 | 0.173 | 0.208 |

### The comparison that counts: a matched operating point

Accuracy cannot be compared across the two sets because prevalence differs (29.5 % vs 21.8 %);
MCC moves with prevalence too. Sensitivity at a fixed specificity does not.

| score | n | specificity | sensitivity | 95% CI |
|---|---|---|---|---|
| **DRSP crystal combined (published)** | 610 | 0.9119 | **0.654** | not reported |
| **MultiPASS, within-protein CV** | 302 | 0.9124 | **0.491** | 0.370–0.611 |
| MultiPASS, LOTO, same 5 genes | 302 | 0.9124 | 0.167 | 0.046–0.296 |
| MultiPASS, LOTO, all 9 genes | 440 | 0.9129 | 0.200 | 0.077–0.323 |

The CI on our 0.491 stops at 0.611, below their 0.654. On the face of it DRSP detects about a
third more of these resistance variants than we do at the same false-positive rate.

### Per gene, LOTO panel

| gene | n | pos | AUC | 95% CI |
|---|---|---|---|---|
| BTK | 27 | 4 | 0.957 | 0.870–1.000 |
| MET | 35 | 7 | 0.954 | 0.857–1.000 |
| FLT3 | 49 | 5 | 0.882 | 0.786–0.968 |
| KIT | 62 | 21 | 0.637 | 0.465–0.796 |
| EGFR | 65 | 15 | 0.617 | 0.447–0.757 |
| ABL1 | 59 | 31 | 0.608 | 0.475–0.750 |
| ESR1 | 27 | 6 | 0.587 | 0.381–0.782 |
| SMO | 51 | 17 | 0.522 | 0.349–0.682 |
| ALK | 65 | 24 | 0.507 | 0.363–0.653 |

The three apparently excellent targets are the three with 4–7 positives; the five with 15–31
positives sit between 0.507 and 0.637. That ordering is what small-sample AUC noise looks like,
not a biological finding, and the report does not offer it as one.

## What this means for the manuscript

The Comparison section currently disposes of DRSP with one clause about web interfaces. That is
no longer defensible now that the comparison has been made. Three statements are supportable:

1. **On accuracy we do not win.** At DRSP's operating point MultiPASS recovers 0.491 of these
   resistance variants against their reported 0.654, and our CI excludes their figure. Any
   sentence implying our method is more accurate than DRSP has to go.
2. **Both figures carry selection optimism, and only ours is measured.** DRSP chose its feature
   set by a genetic algorithm on the training set and kept the best of eight runs by MCC; our
   per-target configuration was chosen by maximum CV over a grid of up to 504 configurations.
   We measured ours (ABL1 0.019, NT5C2 0.134 — `REPORT_ADDITIONS.md`); theirs is unquantified.
   This is a point of methodological strength to make explicitly, not a defence of our number.
3. **The real difference is the input requirement, and that is where the claim belongs.** DRSP
   needs a drug–protein complex per pair: a crystal structure, or a docking model plus iGEMDOCK
   energies. That is why it covers 11 genes and 13 drug pairs and cannot be extended by a user.
   MultiPASS needs a sequence and a substitution. Plus the obvious: DRSP has no server and no
   code, so its numbers cannot be reproduced or applied by anyone — the manuscript's original
   observation, now backed by having tried.

Taken with the other two tier-3 comparators the picture is consistent: on a level playing field
this task is hard for everybody, and where a comparator has a structural or training advantage it
tends to use it. SUSPECT-ABL and MultiPASS were indistinguishable near chance in the ABL1 kinase
domain; Dr. Kinase could not be made into a ranker; DRSP, with a 3D complex per drug and its own
data to fit, is ahead of us at a matched operating point.

## Scope of verification

`verify_drsp.py`, **88 checks, all passing**.

**Verified.** The archive's PK signature, its 167 711 bytes against the server's
content-length, its two-file contents and `Dataset.xlsx`'s 25 750 bytes, with the sha256
recorded. Both sheets' sizes and label composition against the paper's Tables 3 and 4 (610 =
133 + 477; 115 with 3 positives), no duplicate (gene, SAV) rows. The numbering offsets
re-derived independently by brute force, with total wild-type agreement at the chosen offset,
zero agreement at offset 0 for ABL1 and BTK, and the ten canonical ABL1 mutations landing on
their textbook positions and existing in our table. That every one of the 130 positives matched
and all 71 unmatched are negatives. That every matched variant is present in our own variant
table. The label cross-tabulation and the single ESR1 E380Q disagreement. Both AUCs recomputed
from the Mann–Whitney U statistic rather than from `sklearn`. The matched-operating-point
sensitivity recomputed by a brute-force scan over every threshold, with the achieved specificity
checked to be at or above the target and within 0.01 of it, and bootstrap CIs from 2000
stratified resamples. All four of DRSP's published metrics re-derived from their own sensitivity
and specificity, and their headline 0.86/0.57/0.66 confirmed as those values rounded. That our
`cv_rank` scores equal the five-seed rank average the main pipeline produces, for all five CV
targets.

**Not verified.** Whether DRSP's 610 training variants overlap COSMIC v100 — given 99.8 % label
agreement they almost certainly derive from the same annotations, but their sourcing is not
stated variant by variant. The size of their selection optimism, which cannot be estimated
without their code. Whether their 0.654 sensitivity would hold on our 302-variant subset rather
than their full 610 — they publish no per-variant predictions, so the two sensitivities are on
overlapping but not identical sets, and that gap is not quantifiable from the paper. Their
feature computation, iGEMDOCK energies and docking models, none of which was reproduced. The
`NA_deleterious` negatives' provenance: how "deleterious but not resistant" was established for
each is not stated in the file, and 71 of them are not in COSMIC's missense export at all.

## Files

| file | contents |
|---|---|
| `drsp_overlap.csv` | the 440 matched variants: DRSP's SAV and ours, the offset applied, both labels, both scores |
| `drsp_compare.csv` | the three score panels with AUC, CI, AUPRC and threshold metrics |
| `drsp_operating_point.csv` | sensitivity at DRSP's specificity, with bootstrap CIs |
| `drsp_pergene.csv` | the per-gene LOTO panel |
| `drsp_compare.py` | the comparison, including offset detection |
| `verify_drsp.py` | the 88 checks |
