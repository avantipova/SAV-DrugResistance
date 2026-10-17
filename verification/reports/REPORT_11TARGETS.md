# Eleven proteins, and the provenance of Table 1

Generated 2026-10-08. Adds five genes (BTK, ESR1, FLT3, IDH1, MET) to the six already
analysed, and audits where the published Table 1 figures come from, using the project's own
notebook (`aas_2.ipynb`) and selection spreadsheet (`total_table_avg_iap_final_2.xlsx`).

**The one result here that points the other way from everything else in `REPORT_ADDITIONS.md`
is in §3.** It is suggestive, not established, and the reasons are stated with it.

---

## 1. Data intake

The uploaded `csv/` and `mutations/` exports contain 13 genes beyond the original six. All 12
files for the six already analysed are **byte-identical** to the ones every earlier result was
computed from, so nothing downstream shifts.

A target is a stem with all three of `fasta/<stem>.faa`, `csv/<stem>.csv` and
`mutations/<stem>_mutations.csv`. `ABL1` and `ABL1_ENST00000318560` are *different* targets with
different numbering and are never mixed — an earlier version of the validator grouped by gene
name and paired one transcript's labels against another's missense set, which surfaced as ABL1
reporting 0 positives instead of 105.

Which sequence the numbering refers to is decided from the data, by agreement of each
`p.X###Y` wild-type letter with the reference residue. It is decisive: ABL1 on isoform b agrees
**7.23%**, on isoform a **100.00%**. All eleven usable targets reach 100% except FLT3 (99.881%).

One transcript per gene, the one with the most positives (`ingest_check.csv`):

| gene | target | positives (ingest) | positives at L = 25 | positions with both classes | drugs |
|---|---|---|---|---|---|
| ABL1 | ABL1_ENST00000318560 | 105 | 101 | 34 | 5 |
| ALK | ALK | 24 | 24 | 9 | 5 |
| KIT | KIT | 24 | 24 | 14 | 3 |
| SMO | SMO | 18 | 18 | 3 | 1 |
| NT5C2 | NT5C2 | 18 | 18 | 3 | 1 |
| EGFR | EGFR | 17 | 17 | 8 | 5 |
| MET | MET | 8 | 8 | 1 | 4 |
| FLT3 | FLT3 | 7 | 6 | 2 | 3 |
| ESR1 | ESR1_ENST00000206249 | 7 | 7 | 2 | 1 |
| IDH1 | IDH1_ENST00000415913 | 6 | 6 | 3 | 1 |
| BTK | BTK_ENST00000308731 | 5 | 5 | 0 | 1 |

The two positive columns differ for ABL1 (105 → 101) and FLT3 (7 → 6), and for no other target.
The cause is not the window length but **an unspecified mutant residue**: COSMIC records some
substitutions with `X` or `Z` in the mutant position, and a sequence model cannot score a peptide
whose mutant residue is unknown, so `Chem.MolFromSequence` refuses it and the variant is dropped.
Across the whole eleven-target ingest this happens **14 times in 11 836 missense records**, and
**5 of the 14 carry a resistance label**:

| variant | other substitutions at that position | of them resistance-labelled |
|---|---|---|
| ABL1 `p.E279Z` | 3 | 3 — E279K, E279A, E279Y |
| ABL1 `p.E355X` | 2 | 2 — E355G, E355A |
| ABL1 `p.F359X` | 5 | 4 — F359I, F359V, F359C, F359A |
| ABL1 `p.E450X` | 4 | 3 — E450A, E450G, E450K |
| FLT3 `p.D835X` | 9 | 5 — D835Y, D835H, D835F, D835V, D835I |

**No residue is lost.** Every one of the five sits at a position that keeps at least two other
resistance-labelled substitutions, so the filter costs the per-variant task five labels and costs
the residue-level task nothing — which matters, because the residue-level claim is the one the
manuscript makes. The filter is worth one sentence in the Methods; it was previously implicit in
the pipeline and unstated.

**The expansion adds no new per-target row.** Only the original six clear 15 positives. The
within-position test gains 8 positions (71 → 79), which will not rescue it. BRAF, FGFR2, JAK1,
JAK2, MAP2K1, MAP2K2 and PDGFRA have 0–3 positives each and no FASTA, so their sequences are not
worth chasing. What the expansion does buy is §3: ten training proteins instead of five.

Two pipeline changes, both no-ops on the original six:

- **a wild-type guard** — any variant whose wild-type letter disagrees with the reference is
  dropped. It fires exactly once across 7299 variants: FLT3 `p.V571I`, whose reference residue
  is R at 571 and which is not a resistance label. No other FLT3 isoform explains it (the other
  ten records of that multi-record file agree 3–20%).
- **a common window length** (L = 25) for all eleven, since the per-target lengths came from a
  hand-curated shortlist the new genes have no entry in, and a cross-protein base should be
  homogeneous anyway. On ABL1 this gives the same 101 positives as L = 17 did.

Totals: 7299 variants, 234 positives, 11 targets.

---

## 2. Where Table 1's numbers come from

Reconstructed from `aas_2.ipynb` cells 29–34 and `total_table_avg_iap_final_2.xlsx`, whose
`total_table_avg_iap` sheet holds 12 261 rows: peptide lengths 3–37 odd (18 values) × MNA levels
3–30 (28 values) × two charge states — **up to 504 configurations per protein**, 324–885 rows
each. `avg_iap` is MultiPASS's own `Average IAP` from its training log, i.e. internal
validation, not a held-out estimate. My recomputation of the sheet's rule reproduces its own
`avg_iap` on **8589 of 8595 rows** (median |difference| 0.000000), so the metric is read
correctly.

**The selection was not the maximum.** `best_comb` is a hand-curated shortlist with `+` marks
and a column reading "sure" / "not sure" / "script is running". The published configuration is
the shortlist top for EGFR, KIT, NT5C2 and SMO, but **3rd of 5 for ABL1** (L17/MNA9 at 0.85915,
behind L23/MNA7 at 0.87060) and **4th of 6 for ALK** (L33/MNA4 at 0.930, behind L21/MNA5 at
0.94177). The notebook's own two selection runs also disagree with each other and with the
paper — ABL1 at L23/MNA7 in one and L17/MNA10 in the other, ESR1 at MNA 12 and MNA 30 — while
`avg_iap` moves by under 0.004. A flat landscape with an undocumented, partly manual criterion
is harder to defend than either taking the maximum or nesting the search.

Three defects in the metric, all verified against the sheet's own cells (`table1_recomputed.csv`,
`table1_all_configs.csv`):

**(a) Declined classes are dropped from the mean, not scored.** A configuration therefore scores
better for failing on a hard drug. Across the full grid, however, the cost to the *published
value* is small, because with ~670 configurations the search finds another one that does not
decline: at most **0.0119** AUC, mean 0.0070 over the four multi-class targets. What it changes
is *which* configuration wins — in 5 of 13 targets. (An earlier estimate of 0.13 from the
five-row shortlist overstated this; the full-grid figure above supersedes it.) Declines are
nonetheless pervasive: bosutinib yields any value in **25 of ABL1's 674 configurations** —
declined in 96% — nilotinib in 561 of 674, and 2312 of 8597 configurations carry at least one
decline.

**(b) `Not resistant` is averaged in as though it were a drug**, adding the large easy negative
class to a number presented as per-drug performance. Cost to the published value: at most 0.0167
(FLT3), mean 0.0115.

**(c) For 9 of 13 targets the published figure averages exactly one drug class.** EGFR's 0.9602
is osimertinib alone; BTK, ESR1, IDH1, MET, NT5C2 likewise. Only ABL1_ENST (5 classes),
ALK_ENST (3), ABL1 (2) and FLT3 (2) are genuinely multi-drug. Table 1 is not a per-drug average
for most of the table.

**One row is corrupted.** Row 2835, `ABL1ENST00000318560`, L15/MNA3, carries three
ALK-inhibitor values — Crizotinib 0.7746, Alectinib 0.9281, Lorlatinib 0.8868 — for a protein
with no ALK-inhibitor labels, and averages all nine numbers:
(0.8943+0.8707+0.8900+0.8289+0.8391+0.7746+0.9281+0.8868+0.8117)/9 = 0.858244…, which is exactly
why that row's `avg_iap` does not reproduce from its own ABL1 columns. It is a paste
misalignment. It also contaminated my first recomputation, until the attempted-class set was
intersected with the drugs each label file actually names.

So the published figures are **not** materially inflated by (a) or (b). What inflates them, if
anything, is selection over ~670 configurations — which `REPORT_ADDITIONS.md` §7 measures at
0.019 (ABL1, 101 positives) and 0.134 (NT5C2, 18 positives) on a six-point grid, both lower
bounds for a grid this size.

---

## 3. Eleven-protein leave-one-target-out — the one result pointing the other way

One MNA base per held-out protein, trained on all variants of the other ten (6078–7157
variants), scoring every variant of the held-out one. L = 25, MNA 9. No folds, no seeds: the
split is by protein. `loto11.csv`.

| held out | family | n+ / n− | **MultiPASS** | 95% CI | full_gbm | subst_lr | pos_only |
|---|---|---|---|---|---|---|---|
| FLT3 | kinase | 6 / 811 | **0.970** | 0.955–0.983 | 0.903 | 0.855 | 0.854 |
| MET | kinase | 8 / 862 | **0.964** | 0.911–0.994 | 0.917 | 0.773 | 0.363 |
| BTK | kinase | 5 / 343 | 0.792 | 0.438–0.977 | 0.836 | 0.778 | 0.665 |
| ABL1 | kinase | 101 / 608 | 0.734 | 0.678–0.784 | 0.438 | 0.548 | 0.433 |
| ALK | kinase | 24 / 1104 | 0.688 | 0.569–0.790 | 0.606 | 0.627 | 0.555 |
| EGFR | kinase | 17 / 1204 | 0.678 | 0.525–0.814 | 0.825 | 0.685 | 0.580 |
| KIT | kinase | 24 / 935 | 0.668 | 0.534–0.797 | 0.588 | 0.641 | 0.620 |
| SMO | GPCR/Frizzled | 18 / 458 | 0.592 | 0.488–0.691 | 0.743 | 0.565 | 0.614 |
| NT5C2 | 5′-nucleotidase | 18 / 124 | 0.591 | 0.471–0.700 | 0.623 | 0.475 | 0.602 |
| ESR1 | nuclear receptor | 7 / 389 | 0.574 | 0.496–0.661 | 0.352 | 0.789 | 0.128 |
| IDH1 | dehydrogenase | 6 / 227 | 0.455 | 0.173–0.730 | 0.684 | 0.472 | 0.687 |
| **mean** | | | **0.701** | | 0.683 | 0.655 | 0.555 |

Two things follow.

**Cross-protein transfer improves with more training proteins.** Mean 0.604 (five training
proteins) → 0.701 (ten). Five of the six original targets rose: ABL1 0.635→0.734,
ALK 0.599→0.688, EGFR 0.613→0.678, KIT 0.504→0.668, NT5C2 0.573→0.591; SMO fell 0.700→0.592.
**This comparison is confounded**: the six-protein run used each target's own window length, this
one uses L = 25 for all. Separating "more proteins" from "different window" needs a six-protein
run at L = 25, which has not been done. Until then the +0.097 is not attributable.

**Every kinase transfers; no non-kinase does.** All seven kinases land at 0.668–0.970, all four
non-kinases at 0.455–0.592 — **complete separation**, means 0.785 vs 0.553, difference +0.232.
Mann–Whitney one-sided p = 0.0030. It survives removing the two outliers (5 v 4, separation
holds, p = 0.0079) and restricting to targets with ≥15 positives (4 v 2, separation holds,
p = 0.0667).

Four reasons to treat this as suggestive rather than established:

1. **The grouping is post-hoc.** It was chosen after seeing the AUCs. The p-value describes the
   separation; it does not test a pre-registered hypothesis.
2. **In every subset the p-value equals the minimum attainable for those group sizes** (0.0030
   for 7 v 4, 0.0079 for 5 v 4, 0.0667 for 4 v 2). So the evidence is exactly "complete
   separation" and nothing more; with eleven proteins that is as strong as the test can get,
   which is not very strong.
3. **The CIs overlap heavily.** BTK is 0.438–0.977 on five positives; IDH1 is 0.173–0.730 on six.
   Five of eleven CIs include 0.5.
4. **The baselines trend the same way without reaching separation** — full_gbm +0.130
   (p = 0.21), subst_lr +0.126 (p = 0.12), pos_only +0.074 (p = 0.54). The clean separation is
   specific to MultiPASS, but at n = 11 the baselines failing to reach it is weak evidence of a
   real difference.

Two confounds that do **not** explain it:

- **Not sample size.** Spearman ρ between positives and AUC across the eleven targets is
  **−0.08**. FLT3 and MET are not high because they have few positives.
- **Not domain location.** Normalised position does transfer somewhat among kinases (pos_only
  reaches 0.854 on FLT3), but it collapses where MultiPASS is strongest: MET pos_only **0.363**
  against MultiPASS 0.964, ABL1 0.433 against 0.734. Averaged, pos_only's kinase advantage is
  +0.074 against MultiPASS's +0.232.

**And MultiPASS still does not beat the baselines overall**: mean 0.701 against full_gbm 0.683,
ahead in 6 of 11 targets, with only 2 of 33 DeLong tests surviving Holm. The kinase pattern is a
statement about *where* the method transfers, not that it transfers better than a gradient
booster on substitution chemistry.

If it holds up, the claim it supports is narrow and testable: MNA descriptors capture something
about resistance in conserved kinase contexts that does not generalise to other protein
families. The test is more kinases — and the exports already contain FLT3-adjacent candidates
(BRAF, FGFR2, JAK1, JAK2, MAP2K1, MAP2K2, PDGFRA) that need only their reference sequences and
have 0–3 positives each, so they would enter as training data rather than as held-out rows.

---

## Scope of verification

**Verified.** Byte-identity of all 12 re-uploaded files for the original six. Wild-type
agreement for every candidate sequence of every gene. That the wild-type guard fires exactly
once and on a non-label variant. That ABL1/ALK/EGFR/KIT/NT5C2/SMO reproduce their original
positive counts at L = 25. That the recomputed selection metric reproduces the spreadsheet's own
`avg_iap` on 8589 of 8595 rows. The arithmetic of row 2835's nine-value mean. That the LOTO
variant order matches the built tables for all eleven targets (asserted in the script). The
kinase separation under three subsettings, with the minimum attainable p-value reported for each.

**Not verified.** The 0.604 → 0.701 improvement is confounded with window length and is not
attributable until a six-protein run at L = 25 exists. The kinase grouping is post-hoc. MNA 4
has not been run on the eleven-protein set, so the MNA sensitivity shown for six proteins is not
established here. Nested configuration selection has not been run for the five new genes, nor
for ALK, EGFR, KIT or SMO. MultiPASS's internals remain a closed binary. Whether FLT3's and
MET's positives are specifically activation-loop residues homologous to training-set positions
was not checked structurally — the kinase/non-kinase grouping is by protein family, not by
residue homology.

## Files

| file | contents |
|---|---|
| `ingest_check.py`, `ingest_check.csv` | per-target validation: sequence choice, numbering, label hygiene, feasibility |
| `build11.py` | variant tables for the eleven targets, with the wild-type guard |
| `loto11.py`, `loto11_summary.py`, `loto11.csv` | the eleven-protein LOTO and its summary |
| `recompute_table1.py`, `table1_recomputed.csv`, `table1_all_configs.csv` | the Table 1 provenance audit |
