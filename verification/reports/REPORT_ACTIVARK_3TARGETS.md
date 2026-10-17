# ActivArk on ALK, KIT and EGFR — a tie, not a win and not a loss

Generated 2026-10-13. Extension of `REPORT_ACTIVARK.md` to the other kinase targets with enough
positives to support an AUC. **The ABL1 advantage does not reproduce here, but neither does a
deficit: with the models matched, ALK and KIT are ties.** Read the correction notice below
before the tables — the first version of this report got ALK wrong.

**What was run.** All 3 308 of our ALK, KIT and EGFR variants were submitted to
`activark.russelllab.org` (algorithm GBC) in 67 batches of 50 on 2026-10-13; none failed.
BTK (5 positives), FLT3 (6) and MET (8) were excluded beforehand on power grounds — a bootstrap
CI on an AUC with fewer than about fifteen positives is too wide to decide anything — and that
exclusion was decided before any result was seen.

**The arena is different from ABL1's, and weaker.** On ABL1 the referee was the prime-editing
saturation scan, i.e. experimental ground truth neither method had seen. Here the only labels
are **our own** COSMIC/CIViC-derived ones. Both methods were built on curated literature and
database annotations of this kind, so neither arrives naive to the label-generating process, and
no claim of a neutral referee is made. This is a replication check on whether the ABL1 ordering
holds — it does not, in either direction.

> **Correction resolved 2026-10-13.** The first version of this report compared our
> **leave-one-target-out** scores — a model that never saw the target — against ActivArk, which
> was trained on **441 kinases including ALK, KIT and EGFR**, and reported that ActivArk beat us
> on ALK by 0.229 (p = 0.0005). That was a comparison between a target-naive model and a
> target-aware one, and it was unfair to us rather than to ActivArk. It was also the wrong model:
> the web platform serves target-specific models trained on the target itself. The comparison has
> been recomputed with within-target position-grouped cross-validation (`within_target_cv.py`,
> same L = 25, MNA = 9, same windows, same executable — the only variable changed is what the
> model was trained on). **The ALK result does not survive: matched model against matched model
> it is a tie.** The corrected numbers are in the table below; the superseded LOTO figures are
> kept beside them, because the gap between the two columns is itself the finding.

## Coverage, and why EGFR cannot be judged at all

| | ALK | KIT | EGFR |
|---|---|---|---|
| our variants | 1 128 | 959 | 1 221 |
| our positives | 24 | 24 | 17 |
| ActivArk returned an `RvN` | **254 (23 %)** | **493 (51 %)** | **557 (46 %)** |
| of our positives, scored | **24 of 24** | **24 of 24** | **17 of 17** |
| positions scored | 1 087–1 412 | 559–952 | 685–992 |
| **in ActivArk's own resistance training set** | 7 | 16 | **16 of our 17** |

ActivArk scores only its kinase-domain window — **1 304 of the 3 308 variants submitted, 39 %** —
but every one of our positives falls inside it on all three targets, so coverage never costs us a
positive.

**EGFR is not comparable.** Sixteen of our seventeen EGFR positives are in ActivArk's own curated
resistance training set, leaving at most one positive outside it; no honest AUC can be computed.
That is worth reporting as a result in its own right: our EGFR positives are, almost to the
variant, the canonical literature set. It supports the labels and it says there is little novel
EGFR signal for any method to find.

## The comparison, ActivArk's training variants removed

Our out-of-fold Pa against ActivArk's `RvN`, both judged on our binary label, over the
ActivArk-scored variants that are not in its training set. Δ is ActivArk minus our within-target
model, so positive Δ favours ActivArk. Bootstrap, 4 000 resamples, seed 0.

| target | n | positives | ours, within-target | 95 % CI | ours, LOTO | ActivArk | 95 % CI | Δ | 95 % CI | p |
|---|---|---|---|---|---|---|---|---|---|---|
| ALK | 247 | 17 | **0.739** | 0.608–0.851 | 0.476 | 0.707 | 0.571–0.830 | **−0.034** | −0.182…+0.125 | 0.6575 |
| KIT | 477 | 8 | 0.726 | 0.499–0.907 | 0.706 | 0.783 | 0.592–0.948 | +0.060 | −0.133…+0.325 | 0.6895 |
| EGFR | — | ≤1 | — | — | — | — | — | — | — | — |

**Both targets are ties.** On ALK the within-target model reaches 0.739 against ActivArk's 0.707,
and the paired interval straddles zero comfortably (p = 0.66); on KIT, 0.726 against 0.783, also
straddling zero (p = 0.69). With 17 and 8 positives the intervals are 0.24 and 0.41 wide, so
"tie" here means "this data cannot tell them apart", not "they are equal".

The LOTO column is kept because the ALK contrast is instructive: 0.476 target-naive against
0.739 target-aware on the same 247 variants. **Target awareness is worth 0.26 of AUC on ALK** —
more than any difference between the two methods — which is the quantitative reason the first
version of this report reached the wrong conclusion.

## Why our ALK number collapses, and why that is the real finding

The ALK figure is not an artefact of the subset — it is the positional effect, measured directly:

| target | scope | n | positives | our LOTO AUC | 95 % CI |
|---|---|---|---|---|---|
| ALK | whole protein | 1 128 | 24 | **0.6878** | 0.564–0.791 |
| ALK | ActivArk's window, 1 087–1 412 | 254 | 24 | **0.5094** | 0.400–0.617 |
| KIT | whole protein | 959 | 24 | 0.6676 | 0.537–0.791 |
| KIT | ActivArk's window, 559–952 | 493 | 24 | 0.6326 | 0.500–0.762 |
| EGFR | whole protein | 1 221 | 17 | 0.6777 | 0.525–0.811 |
| EGFR | ActivArk's window, 685–992 | 557 | 17 | 0.6231 | 0.473–0.767 |

All 24 ALK positives sit inside 1 087–1 412. Across the whole protein the model separates them
from the many distant negatives easily; **inside the kinase domain that separation is gone.** The
drop is 0.18 for ALK, 0.035 for KIT and 0.055 for EGFR — largest where the domain is the smallest
share of the protein. This is the same conclusion as the prime-editing scan and the internal
`pos_only` baselines, reached a third way, and it is the quantitative statement the hotspot
framing needs: **our whole-protein AUCs are inflated by the easy part of the problem.**

ActivArk, built on kinase-domain alignments and structural features, does not depend on that easy
part — which is exactly why it overtakes us once the easy part is removed.

The two scores remain nearly independent: Spearman(`RvN`, our Pa) is +0.180 on ALK and +0.080 on
KIT. And ActivArk's `RvN` again tracks its own deactivating probability (+0.561 on ALK, +0.608 on
KIT), the same contamination its help page warns about.

## What this means for the manuscript

- **The ABL1 claim still cannot be generalised, but neither can a loss.** On the ABL1 saturation
  scan, judged by experiment, SAV-DrugResistance discriminates imatinib and nilotinib resistance
  better than ActivArk. On ALK and KIT, judged by our own labels inside the kinase domain and with
  the models matched, the two are indistinguishable. The defensible sentence is that
  SAV-DrugResistance is **competitive with the current sequence-level method where both apply**,
  and that is a claim the manuscript can make.
- **Report the whole-protein and in-domain AUCs side by side for every target.** A reviewer will
  restrict to the domain sooner or later; doing it ourselves is the difference between a limitation
  and a retraction risk.
- The honest positioning is unchanged but sharpened: the tool localises resistance to residues
  across proteins **without structures or kinase alignments**, which is where ActivArk cannot go
  at all — it covers 23–51 % of these three proteins and nothing outside a kinase domain, and
  four of our eleven targets (ESR1, IDH1, NT5C2, SMO) are outside its scope entirely. That is the
  defensible claim: broader reach, not better discrimination.
- **The EGFR training-overlap finding should be stated.** It is a clean, checkable fact about the
  field's shared evidence base, and it explains why EGFR cannot referee any comparison.

## The three families of AUC in this project, reconciled

A reader moving between reports will meet three different ALK AUCs and could reasonably call
them contradictory. They are three different quantities, and this is the table that says so.
All three use MultiPASS; they differ in what the model was trained on.

**Whole protein**

| target | n | positives | leave-one-target-out | within-target CV (`REPORT.md`) | protein-wide trained (`margin_sweep`) |
|---|---|---|---|---|---|
| ABL1 | 709 | 101 | 0.734 | 0.857 | 0.860 |
| ALK | 1 128 | 24 | 0.688 | 0.924 | 0.930 |
| EGFR | 1 221 | 17 | 0.678 | 0.921 | 0.922 |
| KIT | 959 | 24 | 0.668 | 0.890 | 0.903 |
| NT5C2 | 142 | 18 | 0.591 | 0.670 | 0.696 |
| SMO | 476 | 18 | 0.592 | 0.805 | 0.824 |

**Inside the ±10 resistance region, same span for all three columns**

| target | span | n | positives | leave-one-target-out | within-target CV | protein-wide trained |
|---|---|---|---|---|---|---|
| ABL1 | 212–495 | 298 | 101 | 0.600 | 0.582 | 0.599 |
| ALK | 1 123–1 269 | 140 | 24 | **0.422** | 0.665 | 0.725 |
| EGFR | 718–797 | 213 | 17 | **0.765** | 0.701 | 0.679 |
| KIT | 578–829 | 312 | 24 | 0.599 | 0.739 | 0.770 |
| NT5C2 | 26–445 | 115 | 18 | 0.569 | 0.620 | 0.664 |
| SMO | 231–535 | 208 | 18 | **0.471** | 0.639 | 0.668 |

Three things follow. **The whole-protein-to-region collapse is real in all three families** —
mean drop 0.087 for leave-one-target-out, 0.187 within-target, 0.172 protein-wide trained — so it
is a property of the task, not of one evaluation. **Leave-one-target-out costs 0.08–0.24 on the
whole protein**, which is the price of a target the model has never seen, and it is why the
figures in this report sit so far below `REPORT.md`'s. And **inside the region the two families
stop agreeing per target**: Spearman is only +0.09 across the six, with ALK and SMO below chance
under leave-one-target-out while EGFR is *higher* in-region than whole-protein, the one target
where the collapse reverses. **The EGFR reversal is now explained** — it is a cross-target
transfer artefact of the target-naive model, which puts its highest EGFR scores on the activation
loop it learned from KIT and FLT3 rather than on EGFR's own ATP-site resistance region; the
within-target model does not make the error and EGFR then drops like every other target
(−0.192). `REPORT_EGFR_ANOMALY.md` has the full trace.

For the manuscript the consequence is simple: **every AUC must say which of the three it is.**
The platform serves target-specific models, so the deployment-relevant figures are the
within-target column; the leave-one-target-out column answers "what if the target is new", which
the manuscript does not claim and should not be read as claiming.

## Scope of verification

**Verified.** That all 3 308 variants were submitted and all 3 308 accounted for, with no batch
failing. That the extracted values are aligned to the right variants: the submitted order was
regenerated independently in the container and checked against a per-gene rolling checksum
computed in the browser — KIT 159113742, EGFR 372431080, ALK 388601878, all three matched — and
the scored counts (493, 557, 254) were reproduced from the parsed stream. That all of our
positives on each target fall inside ActivArk's scored window. Every AUC, CI and paired bootstrap
is recomputed by `ak3_fair.py` (matched models) and `ak3_compare.py` (the superseded
leave-one-target-out version), and the window effect by the table above. The bootstrap is drawn
after a deterministic sort on the variant name, because resampling by positional index made the
third decimal depend on merge order — that was found and fixed while verifying this report.

**Not verified.** Whether the within-target figures would hold under more fold seeds — one seed (0) was used
here, against five in `REPORT.md`, so the within-target column carries fold-assignment variance
this report does not quantify; for scale, the whole-protein values it produces (ALK 0.9014,
KIT 0.8836, EGFR 0.9015) sit 0.006–0.023 below `REPORT.md`'s five-seed means. BTK, FLT3 and MET, excluded on power grounds and not submitted. The `RF`
algorithm option; only GBC was run. Whether ActivArk's curated training set overlaps our ALK and
KIT positives beyond the 7 and 16 variants its own `mutType` field flags. Whether an ensemble of
the two scores beats either. Why ActivArk declines to score 49–77 % of these proteins — the
window boundaries were read off the returned positions, not from its documentation. And the
central caveat above: on these three targets there is no experimental referee, so every number
here is measured against our own labels.
