# Errata in the submitted manuscript, found by checking its own tables

Generated 2026-10-13. Source: Google Drive, `submitted version/Manuscript.docx`, read
2026-10-13. Everything below is an arithmetic or attribution check of the manuscript **against
itself** and against the variant tables the models were trained on. Four defects, all fixable
before resubmission, and one of them is the claim the Comparison section is built on.

## 1. The abstract's headline AUCs do not match Table 1

Table 1 lists 30 target–drug models with an LOOCV and a 20-fold AUC each, and prints its own
average row: **0.920** and **0.919**. Recomputing the mean of the 30 printed values gives
**0.9196** and **0.9191** — so Table 1's average row is right.

The abstract says **0.923** and **0.921**. Neither is the mean of Table 1.

| | abstract | Table 1's average row | recomputed from the 30 rows |
|---|---|---|---|
| mean AUC, LOOCV | 0.923 | 0.920 | **0.9196** |
| mean AUC, 20-fold CV | 0.921 | 0.919 | **0.9191** |

Fix: the abstract should read **0.920** and **0.919**. The same 0.923 appears twice more in the
Comparison section ("sustained a robust mean AUC of 0.923") and must change with it.

## 2. The SUSPECT-ABL comparison compares two different things

The Comparison section reads: *"the structure-guided SUSPECT-ABL predictor, which is exclusively
tailored to the Abelson 1 kinase (ABL1), reports a cross-validated AUC of 0.84. In contrast, our
optimized ABL1-specific SSPR model achieved a significantly higher AUC of 0.923 under LOOCV
validation."*

**No ABL1 model in Table 1 has an AUC of 0.923.** The seven ABL1 rows are 0.934, 0.915, 0.894,
0.890, 0.871, 0.869 and 0.839. The 0.923 is the panel-wide mean across all 11 targets — so the
sentence sets SUSPECT-ABL's ABL1-specific figure against our average over everything, which is a
category error rather than an optimistic rounding.

Two further problems in the same sentence. **"Significantly higher" was never tested** — no
significance test against SUSPECT-ABL appears anywhere in the manuscript. And when SUSPECT-ABL is
actually run on our variants it is **indistinguishable from chance** on them
(`REPORT_SUSPECT_ABL.md`), so the comparison cannot be repaired by swapping in the right number.

Fix: drop the head-to-head claim. What the data supports instead is in `REPORT_ACTIVARK.md` — a
comparison against a method that works, decided on experimental data.

## 3. Table 2, first row: the mean is wrong

| peptide length | MNA | AUC imatinib | AUC TKI class | AUC Not resistant | printed mean | correct mean |
|---|---|---|---|---|---|---|
| 15 | 6 | 0.928 | 0.871 | 0.938 | **0.928** | **0.9123** |
| 33 | 5 | 0.934 | 0.835 | 0.937 | 0.902 | 0.9020 ✓ |
| 11 | 5 | 0.900 | 0.890 | 0.914 | 0.901 | 0.9013 ✓ |

The printed 0.928 is the imatinib cell repeated, not the row mean. The other two rows are
correct. This matters beyond arithmetic: with the mean corrected to 0.912, the row L = 15,
MNA = 6 is no longer the best row of the table — L = 33, MNA = 5 and L = 11, MNA = 5 sit at 0.902
and 0.901, so the ordering is 0.912 > 0.902 > 0.901 and the first row still leads, but by 0.010
rather than the 0.026 the printed value implies.

## 4. The NR column counts variants the models were not trained on

Table 1's `NR` column ("number of drug resistant amino acid substitutions in the training set")
was checked cell by cell against the variant tables the models were built from. **17 of the 19
checkable cells match exactly.** Two do not:

| row | manuscript NR | actually trained on |
|---|---|---|
| ABL1_ENST00000318560 / Imatinib | 87 | **83** |
| FLT3 / Sorafenib | 4 | **3** |

Both are explained completely, and by the same cause. Tracing the ABL1 cell: COSMIC holds **90**
distinct imatinib-labelled AA mutations for this transcript; three of them are not missense
substitutions at all — `p.?`, `p.K294delinsRGG` and `p.L340L` — and are correctly excluded,
giving **87**, exactly the printed NR. Four more have an **unspecified mutant residue** —
`p.E279Z`, `p.E355X`, `p.F359X`, `p.E450X` — which a sequence model cannot score, so the pipeline
drops them, leaving **83**. FLT3 is the same story in miniature: 4 sorafenib-labelled mutations
minus `p.D835X` leaves **3**.

So the NR column is the count *before* that filter, while the AUCs beside it come from models
trained *after* it. Fix: print the trained-on counts, and state the filter in the Methods —
`REPORT_11TARGETS.md` quantifies it at 14 records in 11 836, of which 5 carry a resistance label,
and shows that no residue is lost because every one of the five sits at a position that keeps
other resistance-labelled substitutions.

**One consequence is worth noticing.** The Methods set the inclusion threshold at "at least three
experimentally or clinically validated resistance-conferring substitutions". FLT3 / sorafenib was
admitted on 4 and actually trained on **3** — exactly at the threshold, with no margin. It still
qualifies, but the manuscript should say 3 rather than 4 so a reviewer checking the threshold
reaches the same conclusion.

**A warning about fixing item 1 by find-and-replace.** The string `0.923` occurs four times in
the manuscript, and **one of them is correct**: Table 1's **FLT3 / Quizartinib** row has a
genuine LOOCV AUC of 0.923. The three to change are the abstract ("0.923 and 0.921"), the
Comparison section's "significantly higher AUC of 0.923", and the Comparison section's "sustained
a robust mean AUC of 0.923". A blind replace will corrupt the table.

## Smaller things to reconcile while the file is open

- **"11 core targets" vs "18 target-transcript variations" vs "11 targets and 8 transcripts."**
  All three are defensible but they are never reconciled. Table 1 holds 18 distinct
  target-transcript entries over 11 genes: ten genes appear under their bare name and eight
  appear as transcripts, and IDH1 appears *only* as `IDH1_ENST00000415913`. One sentence in the
  Methods naming that asymmetry removes the apparent contradiction.
- **"More than 80 000 unique SSPR model configurations" (Results) vs "more than 1700 distinct
  SSPR model configurations ... for each protein-drug pair" (Methods.)** The two are consistent
  only if the reader multiplies 1 700 by the number of pairs; say so.
- **Three rows have an LOOCV/20-fold gap above 0.05, and two of them go the wrong way.**
  The median absolute gap over the other **27** rows is **0.002**. In **ALK / Lorlatinib**
  (0.909 → 0.983, **+0.074**) and **ALK_ENST00000618119 / Alectinib** (0.913 → 0.979,
  **+0.066**) the 20-fold value *exceeds* LOOCV, which is the surprising direction: 20-fold
  holds out 5 % at a time and on small, imbalanced sets usually scores at or below LOOCV. In the
  third, **ABL1_ENST00000318560 / Bosutinib** (0.839 → 0.782, **−0.057**), LOOCV exceeds 20-fold,
  which is the expected direction on five positives and needs no explanation. Neither ALK row is
  an error on its face, but both will draw a question and are worth re-deriving.
- The Discussion's FLT3 **D835** paragraph is **sound** and survives the filter above: it names
  D835Y and D835H for sorafenib and D835F, D835V and D835I for quizartinib, and all five are
  present and resistance-labelled in the training data. Only the unspecified `p.D835X` was
  dropped, so the argument stands as written.

## Scope of verification

**Verified.** The 30 rows of Table 1 and the 3 rows of Table 2 were transcribed from the docx and
every mean recomputed; Table 1's printed average row matches its own columns to three decimals
and the abstract's figures do not. That no ABL1 row equals 0.923. That the Comparison section's
0.923 is the panel mean, by comparing it with the ABL1 rows and with the table average. Table 2's
first-row mean, by arithmetic. The NR column against the variant tables for all 19 cells whose
(target, drug) pair exists in them, and the two mismatches traced to named COSMIC records. That
`p.?`, `p.K294delinsRGG` and `p.L340L` are the three non-missense imatinib records. The
LOOCV/20-fold gap for the two ALK rows against the distribution over the other 28 rows. All of it
is recomputed by `check_manuscript.py`.

**Not verified.** The 11 cells of Table 1 whose (target, drug) pair has no counterpart in our
tables — the bare-gene rows for transcripts we ingested under one identifier only — so the NR
check covers 19 of 30 rows, not all 30. The AUC values themselves: these are MultiPASS LOOCV and
20-fold figures from the original pipeline, and this report checks their arithmetic and
attribution, not whether re-running MultiPASS reproduces each cell. The two ALK rows were flagged
from the LOOCV/20-fold gap alone; no attempt was made to re-derive them. The Supplementary
Information and Supplementary Tables 4–6 were not checked here. Whether the journal's submitted
PDF matches this docx.
