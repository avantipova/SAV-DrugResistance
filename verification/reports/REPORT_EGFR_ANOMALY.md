# Why EGFR was the one target whose AUC rose inside the resistance region

Generated 2026-10-13. `REPORT_ACTIVARK_3TARGETS.md` left this unexplained. It now has a complete
answer, and the answer is a limitation worth putting in the manuscript.

**The anomaly.** For every target, restricting the evaluation from the whole protein to the ±10
resistance region lowers the leave-one-target-out AUC — except EGFR, where it **rises**, from
0.6777 to 0.7655. (Earlier drafts of this report wrote "0.678 to 0.766". The 0.766 was a
double-rounding error: the region AUC is exactly 0.7654562 — MultiPASS emits Pa at three
decimals and the parser applies no rounding of its own, verified across 21,608 stored values, so
the per-variant scores are exact and so is the AUC computed from them — which this folder's
`auc_families.csv` stores at 4 dp as 0.7655, and that 0.7655 was then rounded up a second time.
At 3 dp the value is 0.765. It is quoted at 4 dp here and in the manuscript to keep the two
roundings from being confused again.)

**The answer in one line.** The target-naive model transfers an *activation-loop* resistance prior
learned from KIT and FLT3 onto EGFR, whose labelled resistance sits in the ATP site instead. The
mis-scored residues all lie outside the resistance region, so removing them raises the in-region
AUC. The within-target model does not make the error, and EGFR then behaves like every other
target.

## Step 1: the whole-protein AUC decomposes exactly, and EGFR is the only reversal

Splitting the negatives by location gives an exact identity — the whole-protein AUC is the
share-weighted average of the AUC against in-region and against out-of-region negatives. The
identity was checked and reproduces all six values to four decimals.

| target | ±10 span | AUC whole | vs in-region negatives | vs out-of-region negatives | share out |
|---|---|---|---|---|---|
| ABL1 | 212–495 | 0.7340 | 0.6001 | 0.7982 | 0.676 |
| ALK | 1 123–1 269 | 0.6878 | 0.4222 | 0.7190 | 0.895 |
| **EGFR** | **718–797** | **0.6777** | **0.7655** | **0.6607** | 0.837 |
| KIT | 578–829 | 0.6676 | 0.5987 | 0.6983 | 0.692 |
| NT5C2 | 26–445 | 0.5905 | 0.5690 | 0.6677 | 0.218 |
| SMO | 231–535 | 0.5921 | 0.4712 | 0.6779 | 0.585 |

Everywhere except EGFR the out-of-region negatives are the *easy* ones. On EGFR they are the hard
ones, which is the whole anomaly. The mean Pa says the same thing: on five targets the in-region
negatives score higher than the out-of-region ones — nearby negatives resemble positives, the
positional effect — while on EGFR it inverts, 0.0543 in-region against 0.0689 out
(Mann–Whitney, out > in, p = 3.3 × 10⁻⁸).

## Step 2: the high-scoring residues are the activation loop

Mean leave-one-target-out Pa along EGFR in 50-residue windows peaks not at the labelled
resistance region but at **850–899 (0.1034)**, with 800–849 at 0.0923 and 900–949 at 0.0908,
against 0.0620 and 0.0558 for the two windows covering 718–797.

Taking 830–930 as the activation-loop and C-lobe block:

| | n | mean Pa | median Pa |
|---|---|---|---|
| negatives in 830–930 | 172 | **0.1002** | 0.1000 |
| negatives in 718–797 (the labelled region) | 196 | 0.0543 | 0.0500 |
| **our 17 EGFR positives** | 17 | **0.0879** | — |

Mann–Whitney, 830–930 above 718–797: **p = 1.2 × 10⁻³³**. And the activation-loop negatives
outrank **our own positives**: AUC(positives versus 830–930 negatives) = **0.4345**, below chance.

The individual residues identify themselves. The highest-scoring out-of-region substitutions
cluster at positions 811, 826, 829, 832, 848, 854, 860, 861, 863, 864, 865, 866, 914, 915, 919 —
that is **T854, L858, L861, G863, E866**, all present in our EGFR set and all labelled
**non-resistant**, because L858R and L861Q are the canonical *sensitising* mutations. Our 17
positives are instead L718, G719, L747, D761, M766, T790, L792, G796 and C797 — the ATP site and
the gatekeeper, almost all osimertinib.

## Step 3: where the prior came from

The model scoring EGFR was trained on the other ten targets, and two of them keep their
resistance positives in exactly this part of the fold:

| training target | positives | span | median position |
|---|---|---|---|
| KIT | 24 | 578–829 | **816** |
| FLT3 | 6 | 691–835 | **835** |
| ABL1 | 101 | 212–495 | 328 |
| ALK | 24 | 1 123–1 269 | 1 174 |
| MET | 8 | 1 181–1 249 | 1 247 |
| BTK | 5 | 316–481 | 481 |

KIT's D816 and FLT3's D835 are activation-loop hotspots. The model learns that this structural
neighbourhood carries resistance, which is correct for those kinases, and applies it to EGFR —
where the equivalent residues are activating rather than resistance-conferring in our labels.

## Step 4: the decisive test — the within-target model does not make the error

If the explanation is right, a model trained on EGFR itself should not elevate the activation
loop. It does not:

| | 830–930 negatives | 718–797 negatives | our positives | AUC(positives vs 830–930 negatives) | whole protein | inside 718–797 | change |
|---|---|---|---|---|---|---|---|
| leave-one-target-out | **0.1002** | 0.0543 | 0.0879 | **0.4345** | 0.6777 | 0.7655 | **+0.088** |
| within-target | 0.0591 | **0.2114** | **0.3184** | **0.9485** | 0.9015 | 0.7098 | **−0.192** |

The within-target model reverses every part of the pattern: the activation-loop negatives drop to
the bottom, the labelled region and the positives rise to the top, the positives separate from the
loop at 0.9485, and the whole-protein-to-region change flips to **−0.192**, the normal direction
and the largest drop of any target. So the anomaly is **not a property of EGFR** — it is a
property of scoring EGFR with a model that has never seen it.

## What this means for the manuscript

- The anomaly is resolved and should be reported as resolved, not dropped. It is one sentence in
  the Results and a paragraph in the Discussion.
- **The real limitation it exposes:** the positional prior these models rely on is
  *target-specific*, and does not transfer between kinases whose resistance mechanism differs.
  ATP-site resistance (EGFR, osimertinib) and activation-loop resistance (KIT D816, FLT3 D835) are
  different structural stories, and a model trained on one assigns confident scores to the wrong
  part of the other — confident enough to outrank the true positives.
- This supports, rather than undermines, the platform's architecture. The manuscript's own
  argument is that a **single global model is inadequate and target-specific models are required**
  (the hyperparameter-heterogeneity section). This is the sharpest piece of evidence for that
  claim in the whole project, and it is quantitative: 0.4345 against 0.9485 on the same variants,
  the only difference being whether the model saw the target.
- It also bounds the leave-one-target-out numbers reported elsewhere: they are a floor for a
  *new* target, and on a target whose mechanism is unlike the training set they can be worse than
  chance. The platform does not operate in that regime, and the manuscript should say so rather
  than letting a reviewer infer that leave-one-target-out is the deployment number.

## Scope of verification

**Verified.** The decomposition identity, by reproducing all six whole-protein AUCs from the
share-weighted in-region and out-of-region AUCs to four decimals. The per-window Pa profile along
EGFR. The 830–930 versus 718–797 contrast under both models, with Mann–Whitney p-values. That
AUC(positives versus activation-loop negatives) is 0.4345 under leave-one-target-out and 0.9485
within-target, on the same 17 positives and 172 negatives. That T854, L858, L861, G863 and E866
are in our EGFR set and all labelled non-resistant, and that our 17 positives are the nine
ATP-site/gatekeeper positions listed. The training targets' positive spans and medians, from the
variant tables. All of it is recomputed by `verify_egfr_anomaly.py`.

**Not verified.** That the model's internal representation is "activation loop" rather than some
correlated sequence feature — MultiPASS is a closed binary and no feature attribution was done,
so the mechanism is inferred from *where* the scores land, not read out of the model. Whether
L858R and L861Q would be labelled resistant under a different curation (they are sensitising to
first-generation inhibitors and the literature does report acquired resistance at L858R-containing
compound genotypes; our COSMIC extraction labels them non-resistant and that labelling was taken
as given). Whether the same transfer failure occurs on the other targets in the reverse direction
— only EGFR was traced. The within-target column uses one fold seed, so its figures carry
fold-assignment variance that is not quantified here.
