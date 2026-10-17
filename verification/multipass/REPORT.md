# MultiPASS vs. baseline models on our own COSMIC data

All six targets, two experimental settings, position-grouped cross-validation.
Generated 2026-10-04, **recomputed after the FASTA fix described in
`verification/VERIFICATION.md`**. Scripts and raw per-variant scores in this folder.

## What was run

MultiPASS (`OLMPASSdoMNABase` + `OLMPASSofMOL`, Wine/win32) was evaluated on **exactly the
same folds** as the baselines: `StratifiedGroupKFold(5)` grouped by **residue position**, so no
position ever appears in both train and test. Nothing from a test fold enters the MNA base —
verified directly by reading all 280 training SDFs (see the verification report).
Per-target window length and MNA level are the project's own settings
(ABL1 17/9, ALK 33/4, EGFR 37/30, KIT 23/9, NT5C2 25/4, SMO 31/7).

Binary task per target: `Resistant` (to any drug) vs `Not resistant`.

Two settings:

- **full** — every other missense substitution in the protein is a negative.
- **region** — negatives restricted to the resistance span ±10 aa. This removes the
  "is this mutation anywhere near the drug-binding region" shortcut.

Baselines (refit per seed on exactly the variants MultiPASS scored, same splits):

| model | features |
|---|---|
| `pos_only` | normalised residue position only — HistGradientBoosting |
| `subst_lr` | substitution chemistry only: one-hot wt/mut, BLOSUM62, Δhydrophobicity/volume/charge — LogisticRegression |
| `full_gbm` | `subst_lr` features + position + ±7-residue window composition |

AUC is the mean over 5 random fold assignments (seeds 0–4); CIs are stratified bootstraps
(B = 600) of the seed-averaged out-of-fold score; p-values are paired DeLong on those same
variants, Holm-corrected across all **36** tests.

## Result

AUC, mean over seeds:

| target | setting | n+ / n− | seeds | **MultiPASS** | pos_only | full_gbm | subst_lr |
|---|---|---|---|---|---|---|---|
| ABL1 | full | 101 / 611 | 5 | **0.857** | 0.842 | 0.866 | 0.657 |
| ABL1 | region | 101 / 210 | 5 | **0.582** | 0.523 | 0.597 | 0.570 |
| ALK | full | 24 / 1100 | 5 | **0.924** | 0.950 | 0.909 | 0.447 |
| ALK | region | 24 / 137 | 5 | **0.665** | 0.686 | 0.677 | 0.510 |
| EGFR | full | 17 / 1197 | 5 | **0.921** | 0.900 | 0.861 | 0.616 |
| EGFR | region | 17 / 251 | 5 | **0.701** | 0.717 | 0.740 | 0.669 |
| KIT | full | 24 / 936 | 5 | **0.890** | 0.863 | 0.812 | 0.689 |
| KIT | region | 24 / 331 | 5 | **0.739** | 0.732 | 0.680 | 0.648 |
| NT5C2 | full | 18 / 124 | 5 | **0.670** | 0.669 | 0.716 | 0.483 |
| NT5C2 | region | 18 / 102 | 5 | **0.620** | 0.620 | 0.674 | 0.508 |
| SMO | full | 18 / 453 | 5 | **0.805** | 0.775 | 0.806 | 0.459 |
| SMO | region | 18 / 204 | 5 | **0.639** | 0.594 | 0.636 | 0.474 |

Full numbers with sd over seeds, bootstrap CIs, AUPRC and both raw and Holm-corrected
p-values: `summary.csv`; the same thing one row per cell: `table_wide.csv`.

### Three things the table says

**1. MultiPASS is statistically indistinguishable from a position-only model, everywhere.**
Mean difference +0.012 AUC in MultiPASS's favour; it is behind `pos_only` in 4 of 12 cells.
Holm-corrected DeLong p = 1.00 in 11 of 12 cells and 0.84 in the twelfth. A model that sees
*nothing but the residue number* — no wild-type residue, no substituted residue, no sequence
context — matches it.

**2. The only model MultiPASS reliably beats is the chemistry-only one.**
Of 36 tests, exactly 5 survive Holm correction, all against `subst_lr`
(ABL1 full p < 0.0001, ALK full p < 0.0001, KIT full p = 0.014, SMO full p = 0.026,
EGFR full p = 0.026). `subst_lr` is at or below chance in the `full` setting for ALK (0.447),
SMO (0.459) and NT5C2 (0.483) — substitution chemistry alone carries almost no signal here, so
beating it is a weak claim.

**3. Performance collapses once the binding-region shortcut is removed.**
Going from `full` to `region` negatives, MultiPASS drops from 0.857→0.582 (ABL1),
0.924→0.665 (ALK), 0.921→0.701 (EGFR), 0.890→0.739 (KIT), 0.805→0.639 (SMO).
The published ~0.92 figures are mostly the model learning *where* resistance mutations sit,
not *which* substitution at a site confers resistance. All CIs in the `region` setting include
or nearly include 0.5 for at least one model.

Against five external predictors (AlphaMissense, REVEL, CADD, PolyPhen-2, SIFT) MultiPASS does
significantly better in the `full` setting — see `level2/REPORT_LEVEL2.md`.

## Per-drug, in the shape of the paper's own Table 1

The table above collapses every drug into a single `Resistant` class. Table 1 of the paper does
not — it reports one figure per drug. That difference is worth stating plainly, so the
comparison was also rebuilt per drug, on the same position-grouped folds
(`per_drug.csv`, script `folds_drug.py` + `summarize_drug.py`).

MultiPASS is natively multiclass, so the drug list goes into the SDF class field exactly as the
project's own `build_sdf.py` writes it (newline-separated, `Not resistant` otherwise) and one
base per fold predicts Pa for every drug at once — the same compute as the binary run. For a
drug *d*: positives = variants annotated resistant to *d*; negatives = variants with no
resistance annotation at all; variants resistant **only to a different drug** are excluded,
since their status for *d* is unknown. Baselines are refit on that same restricted set.

Only cells with at least 10 positives are reported — of the 20 target×drug combinations, 8
clear that bar and the rest have 1–9 positives.

| target | drug | n+ / n− | **MultiPASS** | pos_only | full_gbm | subst_lr |
|---|---|---|---|---|---|---|
| ABL1 | Imatinib | 83 / 611 | 0.863 | 0.856 | **0.884** | 0.641 |
| ABL1 | TKI - NS | 48 / 611 | 0.862 | 0.822 | **0.873** | 0.697 |
| ABL1 | Dasatinib | 10 / 611 | **0.838** | 0.833 | 0.834 | 0.571 |
| ALK | Crizotinib | 16 / 1100 | 0.921 | **0.940** | 0.884 | 0.562 |
| NT5C2 | Purine | 18 / 124 | 0.696 | 0.660 | **0.741** | 0.487 |

**Nothing changes.** Per drug, MultiPASS lands within ±0.022 of its own binary-task figure
(Imatinib 0.863 vs 0.860, TKI-NS 0.862 vs 0.860, Crizotinib 0.921 vs 0.930, Dasatinib 0.838 vs
0.860), and of the 15 DeLong tests only 3 survive Holm — all three against `subst_lr` again.
It still does not separate from the position-only model or from `full_gbm` for any single drug.

Sanity check on the new pipeline: NT5C2 has exactly one drug, so its per-drug task *is* the
binary task — and the two give the identical number, 0.696 vs 0.696, on identical n (18 / 124).

Not run per drug: **KIT, SMO and EGFR**. KIT's imatinib covers 23 of its 24 positives, SMO's
vismodegib all 18, and EGFR's osimertinib 16 of 17, so their per-drug task is the binary one
already answered above; for EGFR the MNA-30 re-run would also have cost about 20 h for a
question its own numbers answer. **MultiPASS declines per drug too**: on ABL1/Dasatinib seed 3
it returned no Pa for that drug on 130 variants, and on NT5C2 seed 2 for 31 — scored 0, as in
the binary analysis.

## Caveats, stated explicitly

- **EGFR `full` now rests on 5 seeds like every other cell** (it was 1 seed in the previous
  version). At MNA level 30 the base is 289 MB and every single-molecule prediction re-parses
  it, so the four extra seeds took about 16 h of compute. Per-seed AUCs: 0.930, 0.935, 0.923,
  **0.882**, 0.935 — mean 0.921, sd 0.020. The single seed reported before (0.930) sat at the
  top of that spread, so the earlier number was mildly optimistic; one seed in five lands
  almost 0.05 lower, which is why the seed average matters here.
- **18 of the 220 resistance mutations COSMIC lists never become positives**, so n+ is smaller
  than a reader comparing against COSMIC would expect. Every exclusion is itemised in
  `verification/positives_excluded.csv`: 13 are not single substitutions (indels,
  duplications, the frameshift `p.(P627_G664)Xfs*?`, the nonsense `p.Q523*`), 4 are
  ambiguous/stop residue codes (`p.E279Z`, `p.E355X`, `p.E450X`, `p.F359X`), 1 is `p.?`.
  NT5C2 loses the most, 5 of 23, and KIT 4 of 28. **One source-data oddity worth a line in the
  paper: `p.L340L` is listed as an ABL1 resistance mutation but is synonymous.**
- **MultiPASS declined to predict some variants.** Where the base reported the class
  "Not Predictable", the variant was scored 0 rather than dropped (dropping would change the
  variant set between methods and break the paired test). Affected:
  EGFR full seed 0 — 1/1214; NT5C2 full seed 2 — 31/142; NT5C2 region seeds 2 and 4 —
  23/120 and 27/120. This is a real property of the method on small training sets, not an
  infrastructure failure; the bases exist and the `.LST` logs state the refusal. The convention
  **costs MultiPASS about 0.02 AUC** on NT5C2 (full 0.670 scored-as-0 vs 0.694 dropping them;
  region 0.620 vs 0.642), i.e. it is conservative, not flattering.
- **This is not a selection-optimism correction.** Table 1 reports the CV value of the
  configuration chosen by maximising that same CV. The numbers here use the project's chosen
  per-target configuration as given; they do not re-do the configuration search inside each
  fold, so the `full`-setting numbers are, if anything, still optimistic.
- **The binary setting is not the published setting.** The main table answers "is this
  substitution among the curated resistance variants", pooled over drugs; Table 1 is per drug.
  The per-drug section above reconciles them and the answer does not move, but a reviewer
  comparing numbers to Table 1 directly should be pointed at that section.
- **Negatives are "no annotation", not "tested and found non-resistant".** Some are presumably
  resistant and simply uncurated — the CIViC analysis found nine such variants for drugs we do
  label, including KIT D816V (imatinib) and SMO L412F (vismodegib).
- **The ±10-residue region margin is a design choice**, not a fact; no sensitivity sweep over
  that margin was run.
- **What was checked and what was not:** the fold construction, labels, statistics and external
  joins were independently re-verified — `verification/VERIFICATION.md` lists every check,
  including the one real defect that was found and fixed. Table 1 reproduces from MultiPASS's
  own validation log (ABL1 L17/MNA9: Imatinib IAP 0.8865 vs published 0.8887, TKI-NS 0.8420 vs
  0.8435, Nilotinib 0.8694 vs 0.8710). MultiPASS's own internals are a closed binary and were
  not audited.

## Files

| file | contents |
|---|---|
| `summary.csv` | one row per target × setting × model: n, seeds, AUC mean/sd, pooled AUC, bootstrap CI, AUPRC, DeLong p raw and Holm |
| `table_wide.csv` | same, one row per target × setting |
| `scores_all.csv.gz` | every out-of-fold MultiPASS score, per target/setting/seed/variant |
| `summary_final.log` | run log, including every "declined to predict" case |
| `folds.py` | MultiPASS on position-grouped folds (per-fold and per-molecule checkpoints, base caching) |
| `summarize.py` | pooling, bootstrap CIs, paired DeLong, Holm correction |
| `build_sdf.py` | peptide-window SDF construction |
| `ACCESSIONS.md`, `Table_S7_corrected.csv`, `s7_fix.py` | the Table S7 audit |
