# Independent re-verification of the whole comparison

2026-10-04. Runnable audit: `verify_all.py` (prints PASS/FAIL per check and exits with a
summary), plus `loss_audit.py`. Every check below was written to try to break the result, not
to confirm it.

> **Re-audited and corrected 2026-10-08.** Four further audit scripts now exist and all pass:
> `verify_new.py` (86 checks, the power / calibration / margin-sweep / LOTO / nested-search
> additions), `verify_external.py` (49 checks, SUSPECT-ABL and Dr. Kinase),
> `verify_civic.py` (159 checks, `../civic/REPORT_CIVIC.md` against its CSV and logs) and
> `verify_level2.py` (265 checks, `../level2/REPORT_LEVEL2.md` against its CSV, the AlphaMissense
> files and the ESM-2 score tables). Together they found stale numbers in two report files,
> both now fixed and both documented in place: 12 per-target cells plus both pooled tables in
> `REPORT_CIVIC.md`, and both pooled tables plus the test count in `REPORT_LEVEL2.md`. **This
> file was stale in one place too** — the robustness table below predates the EGFR 5-seed re-run
> and is corrected, with the superseded figures shown beside the current ones. The "What changed"
> table that follows is a historical record of the FASTA fix and is left as it was, with one
> annotation where its "after" values have since moved again.

## One real defect was found, fixed, and everything recomputed

**Multi-record FASTA files were read as one concatenated sequence.** Every script read a FASTA
with `''.join(line for line in open(f) if not line.startswith('>'))`. For a single-record file
that is correct. **KIT has 8 records, SMO has 3** (FLT3 has 11, but FLT3 is not analysed), so
for those two targets the "sequence" was every isoform glued end to end: KIT 7792 residues
instead of 976, SMO 2101 instead of 787.

Consequences, measured rather than assumed:

- A peptide window near the end of the canonical sequence ran **into the next record's
  N-terminus** instead of being rejected for running off the end. Affected: **7 of 967 KIT
  variants (positions 971–975) and 11 of 482 SMO variants (773–787)**. All 18 are **negatives**;
  no positive was touched. Zero variants affected in the region setting, and the four
  single-record targets (ABL1, ALK, EGFR, NT5C2) were not affected at all.
- `full_gbm`'s ±7-residue window-composition features and the position normaliser used the
  concatenated length for KIT and SMO.

Fixed by a shared `first_record()` helper in `folds.py`, `summarize.py`, `build_sdf.py`,
`compare.py`, `loto.py` and `verify_all.py`. With the fix those 18 variants are correctly
rejected, so the corrected sets are KIT full 960 and SMO full 471. **MultiPASS was re-run from
scratch for KIT full and SMO full, 5 seeds each**, and `summarize.py`, `level2_eval.py`,
`civic_eval.py` and `civic_overlap.py` were all recomputed.

### What changed

| quantity | before | after |
|---|---|---|
| KIT full, MultiPASS AUC | 0.869 | **0.890** |
| KIT full, pos_only | 0.868 | 0.863 |
| KIT full, full_gbm | 0.791 | 0.812 |
| SMO full, MultiPASS AUC | 0.801 | **0.805** |
| SMO full, subst_lr | 0.469 | 0.459 |
| Level-2 pooled `full`, MultiPASS | 0.844 | **0.848** |
| Level-2 pooled `full`, every external | — | moved by ≤ 0.001 |
| Level-2 pooled `region`, all models | — | moved by ≤ 0.001 |
| CIViC set A pooled, MultiPASS | 0.828 | 0.826 → since **0.821** ‡ |
| CIViC set B pooled, MultiPASS | 0.685 | 0.674 (CI lower bound 0.555 → **0.528**) → since **0.626** (CI lower bound **0.484**) ‡ |
| Holm-significant tests in the main table | 4 | **5** (SMO full `subst_lr` crossed over) |

‡ These two "after" values were themselves superseded on 05-Oct, when the CIViC evaluation was
re-run following the KIT/SMO FASTA re-run and EGFR's fifth seed. The current figures are in
`../civic/REPORT_CIVIC.md`; the consequence is recorded there, and it is the one place in the
project where a conclusion did change — set B's CI now includes 0.5, so MultiPASS is no longer
demonstrably above chance on the discordant positives.

**No conclusion changed by the FASTA fix itself.** MultiPASS still beats all seven external
predictors in the `full` setting (Holm p < 0.001 for each), still does not separate from the
position-only model, and still loses all significance inside the resistance region.

The pre-fix fold directories are kept under `stale_preFASTAfix/` rather than deleted.

## Checks that passed

### Labels and the variant set
- Every row of every `*_mutations.csv` has `Resistant Mutation = Yes`, so treating each listed
  AA mutation as a positive is right. **Checked explicitly — this was a plausible place for a
  silent error and there is none.**
- `scores.csv` in **all 56 runs** reproduces exactly, variant for variant and label for label,
  from an independent rebuild straight out of the source CSV and FASTA.
- Positive-loss fully accounted for: COSMIC lists **220** resistance mutations across the six
  targets and **202** become positives. All 18 exclusions are listed individually in
  `positives_excluded.csv`: 13 are not single substitutions (indels, duplications,
  `p.(P627_G664)Xfs*?`, the nonsense `p.Q523*`), 4 are ambiguous/stop residue codes
  (`p.E279Z`, `p.E355X`, `p.E450X`, `p.F359X`), 1 is `p.?`. **One oddity in the source data:
  `p.L340L` is listed as an ABL1 resistance mutation but is synonymous**, so it is absent from
  the COSMIC missense export and correctly excluded. NT5C2 loses the most (5 of 23) and KIT 4
  of 28 — worth stating in the paper, since a reader comparing n_pos = 18 against COSMIC's 23
  will otherwise wonder.

### Fold integrity — the load-bearing claim
- **280 training SDFs read and checked.** Not one contains a variant from its own test fold.
- Each training SDF's molecule names equal its split's training set exactly — no silent
  substitution.
- **No residue position appears on both sides of any split**, in any of the 56 runs.
- The five test folds partition each variant set exactly once.

### No label leakage into the baselines
- Line-by-line, `y` enters `baselines()` only through `fit()` and the splitter, never through a
  feature.
- Labels permuted 12 times and the baselines refitted: mean AUC **0.486 / 0.501 / 0.490** for
  `pos_only` / `subst_lr` / `full_gbm`. A leak would show up here as an AUC above 0.5.

### Statistics
- The DeLong implementation's AUCs equal `sklearn.roc_auc_score` to **2.2 × 10⁻¹⁶**.
- `p(s vs s) = 1` exactly.
- Null calibration over 3000 simulations: type-I error **0.053** at a nominal 0.05, and the
  p-values are indistinguishable from Uniform(0,1) (KS p = 0.23).
- On real ABL1 data the DeLong p-values match an **independent paired stratified bootstrap**
  (4000 resamples): 0.0744 vs 0.0730, 0.1460 vs 0.1505, 0.0000 vs 0.0000.
- The Holm implementation matches `statsmodels.stats.multitest.multipletests(method='holm')`
  **exactly** (max difference 0.0) over 500 random p-vectors.
- Every one of the 226 reported confidence intervals brackets its own point estimate.
- Rank-averaging across seeds moves AUC by at most 0.029 relative to the mean of the per-seed
  AUCs, and no row by more than 0.03.

### The NaN → 0 choice is conservative, and by how much
MultiPASS declared "Not Predictable" for some NT5C2 and EGFR variants, scored 0 so the paired
tests stay on identical variants. Measured cost to MultiPASS of that choice:

| | NaN → 0 | dropping them | 
|---|---|---|
| NT5C2 full | 0.670 | 0.694 |
| NT5C2 region | 0.620 | 0.642 |

So the convention **understates** MultiPASS by about 0.02 on the only target where it matters.
It is not flattering the method.

### External data joins, verified by a second independent route
- **AlphaMissense:** each per-protein CSV has exactly length × 19 rows for all six targets, and
  covers **4623/4623** of our variants. Cross-checked against AlphaFold DB's *separate*
  genome-keyed `hg38` file: **41 238 values compared, 0 mismatches.** That file also reports the
  transcript behind each accession, which independently confirms the Table S7 correction —
  P00519 → `ENST00000318560.6`, our ABL1 target, not the 1149-residue row.
- **dbNSFP via MyVariant.info:** 18 variants spread across the six genes re-queried live with a
  completely different query shape (per-variant rather than the bulk scroll). **All 18 match
  the stored table exactly** on REVEL, SIFT, PolyPhen-2 and CADD — which validates the scroll
  paging and the per-substitution aggregation, the most bug-prone part.
- AlphaMissense~REVEL Spearman is **+0.75 to +0.85** per target; a scrambled join would be ~0.
- The 0.564 threshold used for the "AlphaMissense calls everything pathogenic" table is not a
  guess: it is exactly AlphaMissense's own class boundary (`LBen` ≤ 0.3400, `Amb` 0.3401–0.5637,
  `LPath` ≥ 0.5640), and recomputing that table from the `am_class` labels gives identical
  percentages.

### CIViC parsing
- Re-parsed with two independent extractors — one substituting fusion prefixes then scanning
  globally, one splitting on ` AND ` and matching each clause anchored. **Both return the same
  90 substitutions, symmetric difference 0.**
- The 90 / 70 / 14 / 6 split reproduces exactly.
- All Resistance evidence is `evidence_type = Predictive`, confirmed rather than assumed.

### Table 1 reproduction
Re-confirmed from MultiPASS's own `BaseValidation` log (`ABL1E17.LST`), not from our
reimplementation: **Imatinib IAP 0.8865** at N = 83 (published 0.8887, N = 87),
**TKI-NS 0.8420** (published 0.8435), **Nilotinib 0.8694** (published 0.8710). The N = 83 vs 87
gap is exactly the four ambiguous-code records. Every per-activity count in that log was
re-derived from the SDF's `Drug` fields and matches: Not resistant 611, Imatinib 83, TKI-NS 48,
Dasatinib 10, Nilotinib 9, Bosutinib 5.

### Cross-report consistency
Originally: 14 numbers quoted in the three reports were looked up in the CSVs they come from.
**Redone far more thoroughly on 08-Oct.** Every numeric cell of every table in
`REPORT.md`, `REPORT_ADDITIONS.md`, `REPORT_11TARGETS.md`, `REPORT_SUSPECT_ABL.md`,
`REPORT_DRKINASE.md`, `civic/REPORT_CIVIC.md` and `level2/REPORT_LEVEL2.md` was compared against
the CSV or log it comes from — 229 numbers in the first four reports (0 mismatches), then 159
checks on the CIViC report and 265 on the level-2 report. **That pass is what found the stale
cells listed in the notice at the top of this file.** The two dedicated scripts
(`verify_civic.py`, `verify_level2.py`) are committed so the check is repeatable rather than a
one-off; both now pass with 0 failures.

## Robustness of the new positive claim

"MultiPASS beats the external predictors" does not depend on how targets are pooled:

Corrected 08-Oct to the current `level2_eval.csv`. The earlier version of this table used the
pre-EGFR-5-seed pooled MultiPASS figure of 0.848 and a macro-average of 0.859, and quoted
`full_gbm` as 0.840 pooled / 0.841 macro; the corrected values are 0.847, 0.858, 0.844 and
0.850. No row changes which model is ahead.

| | pooled (variant-level) | macro-average over the 6 targets | targets where MultiPASS is ahead |
|---|---|---|---|
| vs AlphaMissense | 0.847 vs 0.759 | 0.858 vs 0.785 | 5/6 |
| vs REVEL | 0.847 vs 0.775 | 0.858 vs 0.772 | 6/6 |
| vs CADD | 0.847 vs 0.692 | 0.858 vs 0.708 | 5/6 |
| vs PolyPhen-2 | 0.847 vs 0.670 | 0.858 vs 0.692 | 6/6 |
| vs SIFT | 0.847 vs 0.656 | 0.858 vs 0.679 | 5/6 |
| vs ESM-2 650M | 0.847 vs 0.614 | 0.858 vs 0.662 | 5/6 |
| vs ESM-2 150M | 0.847 vs 0.601 | 0.858 vs 0.642 | 5/6 |
| vs pos_only | 0.847 vs 0.828 | 0.858 vs 0.838 | 5/6 (median ΔAUC +0.0295) |
| vs full_gbm | 0.847 vs 0.844 | 0.858 vs 0.850 | 3/6 (median ΔAUC −0.0005) |

In every row except `pos_only` and `full_gbm`, the one target MultiPASS loses is NT5C2 (18
positives).

## Added after the audit: ESM-2 650M

The 650M checkpoint was run on the author's own machine (the weights are 2.5 GB and the bridge
transfer ceiling was 50–150 MB). Integrity was established the same way as every other external
source rather than assumed: the downloaded file's SHA-256 equals the one HuggingFace publishes
(`a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0`); the run reported 652.4 M
parameters and emitted exactly **119 396** rows, which is 6284 positions × 19 substitutions and
matches the AlphaMissense row count for the same six proteins; coverage of our variants is
**100 %**; Spearman against AlphaMissense is −0.54 to −0.79 (corrected 08-Oct from −0.55) and against REVEL −0.46 to −0.77,
negative by construction since ESM-2's log-ratio runs the other way; and the canonical ABL1
resistance mutations score strongly negative (T315I −7.61, M351T −8.22, E255K −8.97).

It did not change any conclusion. Pooled, ESM-2 goes from 0.601 (150M) to 0.614 (650M) in the
`full` setting and 0.481 to 0.486 in `region` — still the weakest external predictor, still
below chance inside the resistance region, still significantly behind MultiPASS after Holm
across the 140 DeLong tests now in that file.

## What this audit did not verify

- **MultiPASS's internals.** It is a closed binary; the audit verifies what goes in, that
  nothing from a test fold goes in, and what comes out — not the MNA computation itself.
- **EGFR `full` is no longer a single seed** — it was re-run at 5 seeds after this audit
  (0.930, 0.935, 0.923, 0.882, 0.935; mean 0.921, sd 0.020), and the whole chain was
  recomputed. Nothing else in this document changed as a result; the pooled level-2 MultiPASS
  figure moved from 0.848 to 0.847.
- **The trained MSAS/MSAR bases and training SDFs of completed runs were deleted** to free
  disk space for those re-runs. The fold-integrity check above was run and passed against all
  280 training SDFs while they existed, but re-running `verify_all.py` now reports
  "0 SDFs read" for that check rather than re-proving it. `folds.py` regenerates the SDFs
  deterministically, so the check can be reproduced by re-running any fold.
- **The published Table 1 could only be compared against the three activities whose values we
  have.** The full published table was not re-derived.
- **dbNSFP's release version is not pinned** — the scores are whatever MyVariant.info served on
  2026-10-04.
- **Ensembl transcript ids were not checked** (both `rest.ensembl.org` and `www.ensembl.org`
  are unreachable from here), so the retirement of `ENST00000454757` behind the inactive
  `E9PFD7` entry remains unconfirmed.
- **The region definition** (positives' span ±10 residues) is a design choice, not a verified
  fact; results inside it are sensitive to that margin and no sensitivity sweep over the margin
  was run.
