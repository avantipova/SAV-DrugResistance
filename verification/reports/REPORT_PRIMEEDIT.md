# The ABL1 saturation prime-editing scan — experimental ground truth, 1848 substitutions

Generated 2026-10-10. This is not another comparator server. It is the first **experimental**
test set for SAV-DrugResistance, and it answers the two questions no amount of comparison with
other predictors could answer: are our resistance labels real, and is our prediction anything
more than a position lookup.

**Source.** Jung, Yu, Oh, Lee, Jeon, Bae, Sim & Kim, *Comprehensive resistance profiling of
chronic myeloid leukaemia associated ABL1 variants against five tyrosine kinase inhibitors
using prime editing*, **Nature Biomedical Engineering**, October 2025,
[10.1038/s41551-025-01531-4](https://doi.org/10.1038/s41551-025-01531-4). Prime editing
installed 2 802 of the 2 892 possible SNVs in the ABL1 kinase domain in K562 cells, covering
about 98 % of possible single amino acid variants, and each variant pool was grown under
imatinib, nilotinib, dasatinib, bosutinib, ponatinib and asciminib. The per-variant results are
published as code and data at `github.com/Goosang-Yu/CML_VUS`.

**Data used.** `SupplementaryCode4/data/resistance_score/ResistanceScore_K562PE4K_<drug>.csv`,
six files, downloaded 2026-10-10. Each holds 2 244 rows — **1 848 missense**, 290 synonymous,
106 nonsense — with two replicate log fold changes, a normalised LFC per replicate, a
`Resistance_score`, and a three-level `Classification` (Sensitive / Intermediate / Resistant).
All six files carry exactly the same variant set.

## Conclusions first

1. **Our labels are right.** Of our 91 ABL1 resistance positives that the scan covers, **53
   (58.2 %) are experimentally resistant** to at least one ATP-site TKI, against **6.5 %** of
   our 214 negatives and an **8.0 %** base rate over the whole scan. Fisher odds ratio
   **19.9**, p = 6.6 × 10⁻²². This is the strongest independent validation of the COSMIC/CIViC
   labelling in the whole project — stronger than MdrDB's κ = 0.28 — and it is a direct answer
   to the passenger objection.
2. **Our prediction is above chance, and far below 0.92.** The leave-one-target-out model,
   which has never seen an ABL1 variant, reaches AUC **0.620** on imatinib, **0.686** on
   nilotinib, **0.708** on bosutinib, **0.671** on ponatinib, **0.591** on dasatinib, and
   **0.494** — chance — on asciminib. Pooled over the five ATP-site TKIs, **0.612
   (95 % CI 0.565–0.658)**.
3. **All of that signal is positional.** Among the 61 positions where the scan contains both a
   resistant and a sensitive substitution, the model's within-position AUC is **0.498**
   (permutation p = 0.53). A predictor that knows only the position — the mean Pa of its
   substitutions — scores **0.6145**, i.e. **slightly better** than the full model
   (ΔAUC −0.0023, 95 % CI −0.0095 to +0.0051, p = 0.55). The model identifies *where* in the
   kinase domain resistance lives; it does not identify *which* substitution at that position
   confers it.
4. **Asciminib is a mechanistic blind spot.** It is the one allosteric inhibitor in the panel
   (myristoyl pocket, not the ATP site), it has by far the most resistant variants (244 of
   1 848, against 15 for ponatinib), and the model is exactly at chance on it. Resistance to it
   sits in a different part of the protein than our training signal.

## Numbering, checked before anything else

The scan's substitution labels are in **canonical ABL1 numbering** — `T315I` is present,
`T334I` is not — the same numbering as our `ABL1_ENST00000318560` (P00519, 1 130 aa). Brute
force over offsets −30 … +30 against **our reference FASTA**: offset 0 matches the wild-type
residue of **1 848 of 1 848** missense variants; the next best offset matches 147. So no offset
correction is applied, and this was established from the sequence, not from the variant table.

## How the model was applied

The base used is `loto11_L25_mna9_ABL1_ENST00000318560/B.MSAS` — the leave-one-target-out base
trained on the other ten targets, with ABL1 wholly held out. Every one of the 1 848 missense
variants was turned into its 25-residue mutant window from the same FASTA and scored in one
`OLMPASSofSDF` pass; no variant failed to produce a result file.

Two checks on the mechanics. The 25-residue windows of all **709** of our own ABL1 variants
were rebuilt from the FASTA and compared with the `pep` column actually used in the LOTO run:
**709 of 709 identical**. And for the **305** variants that appear both in our set and in the
scan, the Pa values from this run and from the LOTO run agree to **max |difference| 0.000000**,
0 differing rows. So the numbers below come from the same model, the same window and the same
executable as Table 1.

Coverage: our ABL1 set has 709 missense variants and 101 resistance positives. All 101 fall
inside the scan's window (positions 184–505); **91** are substitutions the scan actually made,
and 305 of our 709 variants overlap in total. The 1 848-variant evaluation below is **not**
restricted to that overlap — it uses the whole scan, which is the point: a saturation library is
not enriched for variants anyone has observed, so it cannot reward a model for recognising
famous mutations.

## Per-drug performance on the full scan

Resistant vs Sensitive, Intermediate excluded as the source file defines it. Bootstrap CIs,
4 000 resamples stratified on the label, seed 0. ρ is Spearman of our Pa against their
continuous `Resistance_score` over all 1 848 variants.

| drug | resistant | sensitive | AUC | 95 % CI | ρ | p(ρ) |
|---|---|---|---|---|---|---|
| Bosutinib | 15 | 1 546 | **0.708** | 0.563–0.840 | +0.053 | 0.022 |
| Nilotinib | 81 | 1 517 | **0.686** | 0.630–0.740 | +0.087 | 1.7e-4 |
| Ponatinib | 15 | 1 554 | **0.671** | 0.568–0.768 | +0.051 | 0.029 |
| Imatinib | 137 | 1 449 | **0.620** | 0.573–0.668 | +0.080 | 5.4e-4 |
| Dasatinib | 33 | 1 550 | 0.591 | **0.474–0.707** | +0.025 | 0.29 |
| Asciminib | 244 | 1 338 | 0.494 | **0.453–0.535** | −0.009 | 0.71 |

Four of the six confidence intervals exclude 0.5; dasatinib's and asciminib's do not. Pooling:

| label | n resistant | n sensitive | AUC | 95 % CI |
|---|---|---|---|---|
| resistant to ≥1 of the five ATP-site TKIs | 147 | 1 201 | 0.612 | 0.565–0.658 |
| resistant to ≥1 of all six drugs | 294 | 1 038 | 0.525 | 0.487–0.562 |

The six-drug pool is at chance **because** asciminib contributes 244 of its 294 positives. That
is a result about the panel, not about the pooling.

The conclusion does not depend on their three-level convention. Scoring `Resistance_score ≥ 2.5`
as resistant instead gives AUC 0.593 for imatinib (169 positives), 0.658 for nilotinib (107) and
0.494 for asciminib (250).

## The position-controlled test, which is the important one

A saturation scan makes a test possible that no observational dataset allows: hold the position
fixed and ask whether the model ranks the resistant substitution above the sensitive one at that
same position. Positional information is then worth nothing by construction.

| quantity | value |
|---|---|
| positions holding both a resistant and a sensitive substitution (ATP-site pool) | 61 |
| discordant/concordant pairs available | 250 (9 ties) |
| **pooled within-position AUC** | **0.4980** |
| unweighted mean over the 61 positions | 0.5211 |
| permutation p, labels shuffled inside each position (2 000 draws) | **0.53** |
| permutation null mean (95 %) | 0.5007 (0.412–0.588) |
| position-mean Pa alone, one value per position | **0.6145** (CI 0.567–0.659) |
| Pa minus its own position mean | 0.5120 |
| AUC(full Pa) − AUC(position-mean Pa) | **−0.0023** (CI −0.0095…+0.0051), p = 0.55 |

Within-position Spearman against the continuous score tells the same story: imatinib +0.019
(p = 0.41), nilotinib +0.037 (p = 0.12), asciminib +0.023 (p = 0.32), over 318 positions each,
against overall ρ of +0.080, +0.087 and −0.009.

This is the same conclusion the internal `pos_only` and `nearest_pos` baselines pointed at and
that MdrDB reached from the other direction (MultiPASS AUC 0.716 vs the binary label alone
0.712) — but here it is settled against experimental data, with the confound removed by the
design of the experiment rather than by an argument.

## Where the 38 discordant positives went

Of our 91 positives in the scan, 53 are experimentally resistant to an ATP-site TKI. The other
38 break down as: **4** resistant to asciminib only, **13** called Intermediate rather than
Resistant for at least one ATP-site TKI, and **24** sensitive to all six drugs in K562 cells.

Softening the experimental call to "Resistant or Intermediate" raises our positives to
**66/91 = 72.5 %** against **75/214 = 35.0 %** for our negatives — which is exactly the scan's
own base rate of 35.0 %, so under the softer definition our negatives are indistinguishable
from a random substitution while our positives remain enriched twofold. The hard definition is
the more informative one and is what the headline numbers use.

Twenty-four fully sensitive positives is a real false-positive rate in the labels, and the
honest reading is that CIViC/COSMIC evidence of clinical association does not always survive a
cell-line growth assay — K562 is one genetic background, one readout, one drug concentration per
arm.

## What this means for the manuscript

This dataset should go into the paper, and it changes what the paper can claim.

- It supplies the external experimental validation the desk rejections asked for, on 1 848
  substitutions, from a 2025 *Nature Biomedical Engineering* paper, with the model blinded to
  the target by leave-one-target-out.
- The label-validation result (OR 19.9, p = 6.6 × 10⁻²²) belongs in the Results as its own
  paragraph. It is the answer to "your positives may be passengers".
- The honest performance claim for ABL1 against experiment is **0.59–0.71 per drug**, not 0.92.
  The 0.92 figure is a within-COSMIC cross-validation number and should never be set beside an
  experimental AUC without that distinction being made explicitly.
- The positional finding should be stated, not hidden. "SAV-DrugResistance localises resistance
  to residues rather than discriminating substitutions at a residue" is a defensible and useful
  claim for a web tool; letting a reviewer discover it in a saturation scan is not.
- Asciminib should be named as a limitation. The model is at chance on the one allosteric
  inhibitor tested, which bounds the claim of drug-agnostic prediction.

## Scope of verification

**Verified.** The six source files by byte size against the repository tree listing, and
`RS_Imatinib.csv` by SHA-256 (`c0771bf1f347…f5148e31`). That all six carry an identical 1 848-variant
missense set. That the numbering is canonical, by 1 848/1 848 wild-type agreement with our
reference FASTA at offset 0 against 147 for the next-best offset, and by `T315I` present /
`T334I` absent. That the 25-residue windows are built identically to the LOTO run, by rebuilding
all 709 of our ABL1 windows from the FASTA (709/709 identical). That the scoring pipeline is the
same one, by Pa agreement to 0.000000 on all 305 overlapping variants. That every one of the
1 848 variants produced a result file. All AUCs, CIs, Spearman coefficients, the permutation
test, the Fisher test and the paired bootstrap are recomputed from the downloaded files by
`verify_primeedit.py`.

**Not verified.** The paper's headline count of **361 resistance-conferring variant–TKI pairs**
could not be reproduced from these files: the `Classification` column marks **525** resistant
pairs (all of them missense), and the three-level call is not a plain threshold on
`Resistance_score` — Intermediate reaches 5.92 while Resistant starts at 2.76 — so the published
361 evidently applies a further criterion (replicate agreement, or a MAGeCK test) that these
files do not expose. Everything above therefore uses the repository's own `Classification`
column, and the threshold-based sensitivity check is reported for that reason. Also not
verified: that the repository files are the same ones behind the paper's figures, which was not
checkable without the Supplementary Tables; whether the 90 SNVs the library failed to install
are non-random with respect to resistance; and whether a model trained *with* ABL1 would do
better on the scan — only the leave-one-target-out model was applied, which is the conservative
choice.

## Files

| file | contents |
|---|---|
| `pe_scan_pa.csv` | all 1 848 scan missense variants with our LOTO Pa |
| `pe_pooled.csv` | per-variant any-drug-resistant / all-drug-sensitive flags |
| `pe_eval.json`, `pe_eval2.json`, `pe_eval3.json`, `pe_eval4.json` | every number in this report |
