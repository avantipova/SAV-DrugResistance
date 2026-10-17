# MultiPASS against seven external variant-effect predictors

Generated 2026-10-04, **recomputed after the FASTA fix described in `../verification/VERIFICATION.md`**. Script `level2_eval.py`, full numbers `level2_eval.csv`, log `level2_eval.log`.

> **Numbers corrected 2026-10-08.** The two **pooled** tables in this file had been written from
> an earlier run and did not match `level2_eval.csv` / `level2_eval.log`: the pooled negative
> count (4365 instead of 4347), several AUPRCs and confidence-interval bounds, the test count
> (112 instead of 140), and — the one that matters — the MultiPASS-vs-`pos_only` raw p in the
> `full` setting, printed as 0.058 when the run gives **0.027**. The six **per-target** rows were
> checked cell by cell and were already correct, as was the AlphaMissense pathogenic-share table
> (recomputed from the AlphaFold DB files on 08-Oct: 90/80, 88/72, 94/73, 92/64, 100/56, 89/64,
> 0 variants missing a score). Everything below is now read from the CSV.
>
> **What this changes in substance:** nothing about the external predictors, but the headline
> caveat has to be worded more carefully. MultiPASS is nominally ahead of the position-only
> model at raw p = 0.027, not p = 0.058. That still does not survive Holm correction over the
> 140 tests in this file (Holm = 1.00), so the conclusion "not separated from a model that sees
> only the residue number" stands — but it now rests on the multiplicity correction rather than
> on the raw p being above 0.05, and the paper should say so in those words.

## What is compared

**Trained on our labels**, scored out-of-fold on position-grouped folds and rank-averaged over
the 5 fold-assignment seeds, exactly as in the main table — nothing refitted here:
MultiPASS, `pos_only`, `subst_lr`, `full_gbm`.

**External, never trained on drug-resistance labels and not refitted at all:**

| predictor | source | direction |
|---|---|---|
| AlphaMissense | `am_pathogenicity`, per-protein CSV from AlphaFold DB (`amAnnotationsUrl`), CC BY-NC-SA 4.0 | higher = more pathogenic |
| REVEL | dbNSFP via MyVariant.info | higher = more deleterious |
| CADD | phred, via MyVariant.info | higher = more deleterious |
| PolyPhen-2 | HDIV score, via MyVariant.info | higher = more damaging |
| SIFT | via MyVariant.info; the score used is **−SIFT** | lower = more deleterious |
| ESM-2 650M | `esm2_t33_650M_UR50D`, wt-marginal log-ratio log p(mut) − log p(wt) over the full canonical protein; the score used is **−LLR** | lower = more deleterious |
| ESM-2 150M | the same from `esm2_t30_150M_UR50D`, kept alongside to show the model-size effect | lower = more deleterious |

**ESM-2 is now the full 650M model** (`esm2_t33_650M_UR50D`, 652.4 M parameters), with the
150M checkpoint kept beside it. The weights could not be moved into this environment — the
bridge transfer ceiling sat between 50 and 150 MB while the file is 2.5 GB — so the 650M run
was executed on the author's own machine from a locally downloaded copy whose SHA-256 matches
the one HuggingFace publishes (`a08adabb…`). Scoring is wt-marginals (Meier et al., ESM-1v):
one forward pass over the unmasked wild-type sequence, reading the masked-LM head at each
position, on the FULL protein rather than the peptide window the MNA base uses — the window is
a MultiPASS implementation detail, not a property of the variant. The three proteins longer
than ESM-2's 1022-token limit (ABL1, ALK, EGFR) are scored in overlapping windows, each
position taken from the window where it sits furthest from an edge. Join checks: 119 396 rows,
**100 % coverage** of our variants, Spearman against AlphaMissense −0.54 to −0.79 (recomputed
08-Oct; an earlier version of this line said −0.55, rounding KIT's −0.539 the wrong way) and against
REVEL −0.46 to −0.77 (negative by construction, since ESM-2's LLR runs the other way), and the
canonical ABL1 resistance mutations score strongly negative (T315I −7.61, M351T −8.22,
E255K −8.97).

Every model is scored on the **same** variants — the intersection where all eight have a
value — so each paired DeLong test runs on identical variants. Variants dropped for a missing
external score: 16 ABL1, 10 ALK, 34 EGFR, 21 KIT, 0 NT5C2, 3 SMO (dbNSFP holds only SNVs, so
substitutions needing two nucleotide changes are absent). dbNSFP gives per-nucleotide variants,
so several genomic variants can produce one substitution; scores are aggregated per
substitution as max REVEL / max CADD / max PolyPhen-2 / min SIFT.

**Join verification** (not assumed — checked): AlphaMissense's `protein_variant` keys are the
canonical UniProt sequences whose CRC64 we already matched, and each per-protein CSV has
exactly length × 19 rows (**21 470**/19 = 1130 for ABL1, **30 780**/19 = 1620 for ALK, and so on;
119 396 rows over the six files, re-counted 08-Oct — an earlier version of this line quoted
21 471 and 30 781, which are the line counts including the CSV header and are not divisible by
19). Coverage of our variants by AlphaMissense is 4623/4623 (re-checked 08-Oct, per target
712/1124/1214/960/142/471), and the per-protein CSVs were cross-checked value-for-value against AlphaFold DB's separate genome-keyed file (41 238 values, 0 mismatches). Spearman AlphaMissense~REVEL is
+0.75 to +0.85 per target and AlphaMissense~CADD +0.64 to +0.78, which a scrambled join would
not produce. dbNSFP's `aa.pos` arrives as an array over isoforms (T315I appears as
`[334, 315]`, the +19 ABL1 offset), so a row was emitted for **every** element and the join
lands on our numbering without assuming which element is ours. Spot check: the eight canonical
imatinib-resistance mutations all carry high AlphaMissense scores (T315I 0.9991, Y253H 0.9999,
M351T 0.9962, E255K 0.990, G250E 0.9919, F359V 0.9666, M244V 0.929, H396R 0.7922) at FASTA
positions whose wild-type residue matches.

## Result 1 — against protein-wide negatives, MultiPASS beats all seven

Pooled over the six targets (192 pos / 4347 neg), `full` setting:

| model | kind | AUC | 95% CI | AUPRC | DeLong p | Holm |
|---|---|---|---|---|---|---|
| **MultiPASS** | trained here | **0.847** | 0.824–0.869 | 0.239 | — | — |
| full_gbm | trained here | 0.844 | 0.820–0.866 | 0.214 | 0.771 | 1.00 |
| pos_only | trained here | 0.828 | 0.803–0.853 | 0.172 | 0.027 | 1.00 |
| REVEL | external | 0.775 | 0.751–0.798 | 0.097 | <0.0001 | **<0.001** |
| AlphaMissense | external | 0.759 | 0.735–0.784 | 0.091 | <0.0001 | **<0.001** |
| CADD | external | 0.692 | 0.661–0.718 | 0.068 | <0.0001 | **<0.001** |
| PolyPhen-2 | external | 0.670 | 0.635–0.704 | 0.075 | <0.0001 | **<0.001** |
| SIFT (−) | external | 0.656 | 0.619–0.692 | 0.070 | <0.0001 | **<0.001** |
| ESM-2 650M | external | 0.614 | 0.578–0.651 | 0.056 | <0.0001 | **<0.001** |
| ESM-2 150M | external | 0.601 | 0.562–0.641 | 0.053 | <0.0001 | **<0.001** |
| subst_lr | trained here | 0.596 | 0.549–0.648 | 0.066 | <0.0001 | **<0.001** |

**This is the answer to the reviewer's question.** On the project's own data, MultiPASS
significantly outperforms every general-purpose missense predictor, all seven, after Holm
correction across the 140 DeLong tests in this file. The margin on AUPRC is larger than on AUC
(0.239 against 0.053–0.097 for the seven externals), which matters at this class imbalance.

The honest other half stands: **MultiPASS does not separate from a model that sees only the
residue number** (0.847 vs 0.828, raw p = **0.027**, Holm 1.00). Note the raw p is below 0.05
here — the non-separation claim rests on the multiplicity correction over 140 tests, not on the
raw p-value, and should be stated that way rather than as "p > 0.05".

Per target, `full` (AUC; full table in the CSV):

| target | n (pos/neg) | MultiPASS | pos_only | full_gbm | REVEL | AlphaMissense | CADD | PolyPhen-2 | −SIFT | −ESM-2 650M |
|---|---|---|---|---|---|---|---|---|---|---|
| ABL1 | 696 (93/603) | 0.862 | 0.839 | **0.873** | 0.800 | 0.741 | 0.691 | 0.622 | 0.617 | 0.547 |
| ALK | 1114 (24/1090) | 0.930 | **0.951** | 0.920 | 0.892 | 0.812 | 0.638 | 0.714 | 0.748 | 0.730 |
| EGFR | 1180 (16/1164) | **0.927** | 0.924 | 0.889 | 0.728 | 0.768 | 0.698 | 0.761 | 0.728 | 0.716 |
| KIT | 939 (23/916) | **0.908** | 0.869 | 0.831 | 0.804 | 0.816 | 0.779 | 0.730 | 0.656 | 0.576 |
| NT5C2 | 142 (18/124) | 0.696 | 0.660 | 0.741 | 0.659 | **0.774** | 0.750 | 0.650 | 0.704 | 0.720 |
| SMO | 468 (18/450) | 0.823 | 0.785 | **0.844** | 0.746 | 0.801 | 0.692 | 0.673 | 0.622 | 0.685 |

## Result 2 — inside the resistance region, everything collapses, external tools worst

Pooled, `region` setting (192 pos / 1206 neg):

| model | kind | AUC | 95% CI | DeLong p | Holm |
|---|---|---|---|---|---|
| full_gbm | trained here | 0.653 | 0.609–0.696 | 0.3605 | 1.00 |
| **MultiPASS** | trained here | **0.635** | 0.591–0.680 | — | — |
| pos_only | trained here | 0.607 | 0.566–0.649 | 0.084 | 1.00 |
| REVEL | external | 0.571 | 0.530–0.612 | 0.027 | 1.00 |
| subst_lr | trained here | 0.566 | 0.521–0.618 | 0.029 | 1.00 |
| AlphaMissense | external | 0.566 | 0.528–0.602 | 0.016 | 1.00 |
| PolyPhen-2 | external | 0.553 | 0.508–0.597 | 0.0056 | 0.510 |
| SIFT (−) | external | 0.548 | 0.504–0.586 | 0.0064 | 0.576 |
| CADD | external | 0.533 | 0.485–0.576 | 0.0007 | 0.073 |
| ESM-2 650M | external | **0.486** | 0.443–0.524 | <0.0001 | **<0.001** |
| ESM-2 150M | external | **0.481** | 0.443–0.519 | <0.0001 | **<0.001** |

**Almost nothing is significant after Holm — with one exception: ESM-2.** At 0.486 pooled
(650M; 0.481 at 150M) it is the only predictor MultiPASS significantly beats inside the region
(Holm p < 0.001), and the only one whose pooled AUC is below 0.5. Per target the 650M model is
below chance on ABL1 (0.402, Holm 0.044), ALK (0.439) and KIT (0.484). On ABL1 and ALK the other external predictors also go *below*
chance inside the region — AlphaMissense 0.471 and CADD 0.471 on ABL1, CADD 0.377 and
AlphaMissense 0.461 on ALK — i.e. they actively anti-rank.

The mechanism is visible directly. Share of variants AlphaMissense calls pathogenic
(`am_pathogenicity` > 0.564), region setting:

| target | resistance positives | region negatives |
|---|---|---|
| ABL1 | 90 % | 80 % |
| ALK | 88 % | 72 % |
| EGFR | 94 % | 73 % |
| KIT | 92 % | 64 % |
| NT5C2 | 100 % | 56 % |
| SMO | 89 % | 64 % |

AlphaMissense flags almost every substitution in a kinase domain as pathogenic, resistance or
not, so inside the region it has almost no discriminating power left. This is expected rather
than a failure of the tool: resistance mutations must *preserve* kinase function to be
selected, so "generally damaging" is close to the wrong question. It is the cleanest argument
in the paper for why a resistance-specific model is worth having at all — and it is an argument
about the *task*, not about MNA descriptors, since the position-only model shares it.

## Caveats

- **The comparison is not symmetric.** MultiPASS and the baselines are fitted to these
  resistance labels (out-of-fold); the five external tools have never seen a resistance label
  and are applied off the shelf. They are disadvantaged on task specificity and advantaged in
  not being fitted to this data at all. The fair reading of Result 1 is "a model trained on
  resistance beats general pathogenicity predictors on resistance", not "MultiPASS is a better
  predictor".
- **Result 1's negatives are protein-wide**, so part of what every model exploits is domain
  location. Result 2 removes that and the gap to the external tools closes.
- **AlphaMissense and REVEL were trained on data that may include ClinVar**; none of these
  targets' resistance variants are a designated training set, but circularity cannot be fully
  excluded for the pathogenicity tools.
- **n is small.** 192 positives in total, 16–93 per target. Several per-target comparisons have
  CIs spanning 0.2 AUC.
- **EGFR `full` is now 5 seeds for MultiPASS too**, so no cell in this table rests on a single
  fold assignment any more. Its per-seed AUCs spread 0.882–0.935.
- **dbNSFP was read through MyVariant.info, not from the release files**, so the version is
  whatever MyVariant currently serves; the release version is not pinned here.
- **Model size is not what holds ESM-2 back.** Going from 150M to 650M moves the pooled figure
  from 0.601 to 0.614 in the `full` setting and from 0.481 to 0.486 in `region` — it stays the
  weakest external predictor and stays below chance inside the region. Whatever ESM-2 is missing
  here is not capacity.
- **The therapy-matched analysis is now done** and lives in `../civic/REPORT_CIVIC.md`.
- **The windowing confound is now partly resolved.** At 150M the two best targets were the two
  that fit in a single context (NT5C2, SMO), which confounded "weak model" with "split
  context". At 650M that pattern breaks: ALK (0.730) and EGFR (0.716) are windowed and score
  well, while KIT (0.576) and ABL1 (0.547) are its worst despite KIT fitting in one window. So
  the weakness does not track windowing.
- **Robustness of the headline:** the claim does not depend on pooling. Macro-averaging the six
  per-target AUCs equally gives MultiPASS **0.858** against AlphaMissense 0.785, REVEL 0.772,
  CADD 0.708, PolyPhen-2 0.692, SIFT 0.679, ESM-2 650M 0.662 and ESM-2 150M 0.642, and
  MultiPASS is ahead in 5 or 6 of the 6 targets against every external (6/6 against REVEL and
  PolyPhen-2; the one target it loses to the others is always NT5C2). Re-running EGFR at 5 seeds
  moved the pooled MultiPASS figure by 0.001. Against `pos_only` it is ahead in 5/6 with a median
  ΔAUC of **+0.0295**, and against `full_gbm` in 3/6 with a median of **−0.0005** — i.e. the two
  are indistinguishable by this measure as well. (These macro figures are averages of the
  3-decimal per-target values in `level2_eval.csv`.)

## Data provenance

| file | source | licence |
|---|---|---|
| `data/alphamissense/AF-*-F1-aa-substitutions.csv` | AlphaFold DB `amAnnotationsUrl`, 6 accessions | CC BY-NC-SA 4.0 |
| `data/dbnsfp/dbnsfp_myvariant_6genes.tsv` | MyVariant.info, `dbnsfp.genename:<gene>` scrolled in full, 71 749 substitution rows | dbNSFP terms |
| `data/civic/01-Oct-2026-*.tsv` | CIViC release 01-Oct-2026 | CC0 |

The full AlphaMissense release (`AlphaMissense_aa_substitutions.tsv.gz`, 1151 MB,
[10.5281/zenodo.10813168](https://zenodo.org/records/10813168), md5
`b9ccb339e0de6cb0a8d1973ad2026576`) was **not** needed: AlphaFold DB serves the same
predictions per protein, 2.25 MB for all six. dbNSFP's own release was not needed either.
