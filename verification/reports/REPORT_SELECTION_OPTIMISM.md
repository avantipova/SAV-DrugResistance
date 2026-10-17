# How much of Table 1 is configuration-selection optimism

Generated 2026-10-10. Closes §5.2 of `SAV-DR_handover.md` — the one item on that plan that had
never been computed — for all eleven targets.

**The question.** Table 1 reports, for each target, the cross-validated accuracy of the
configuration (peptide window length L, MNA level, charge) that was chosen *by maximising that
same cross-validated accuracy*, over a grid of up to 674 configurations per protein
(`table1_recomputed.csv`). The number published is therefore a maximum over many correlated
estimates, and a maximum is biased upwards. The honest quantity re-does the configuration search
*inside* each training fold, so the configuration is never picked using the variants it is then
scored on. The difference is the optimism.

**The answer, in one line.** Where it can be measured it is **0.0010 to 0.1337 of AUC**, and on a
grid of only six configurations — a hundredth of the real one — the spread between the best and
worst configuration is already **0.04 at 101 positives and 0.11–0.36 at 7–24 positives**.

> **Superseded 2026-10-16 — see the addendum at the end of this file.** This line originally read
> "0.019 to 0.134", said the optimism "is larger where resistant variants are fewer", and counted
> three targets corrected and four unmeasurable. All six targets that can be run are now done:
> the range is 0.0010 to 0.1337 over **six** targets, **three** are provably not computable, and
> the dependence on the number of resistant variants is gone — Spearman -0.12, exact p = 0.84,
> against -0.63 when only five targets were measured. The clause about positives has been deleted
> from the line above rather than rewritten, because it is the second claim in this report to die
> on more data.

## Method

`nested_config.py`. Grid L ∈ {17, 25, 33} × MNA ∈ {4, 9}, six configurations. Outer 5 folds,
inner 3 folds, both `StratifiedGroupKFold` grouped by residue position, one outer seed. Every
configuration is scored on the **intersection of the variant sets across L**, because a longer
window drops variants near the protein's ends and configurations must otherwise be compared on
different data.

- **flat** — one full outer CV per configuration; `flat_max` is the quantity Table 1 reports.
- **nested** — per outer fold, the configuration is chosen by an inner CV on that fold's training
  variants only, then applied to the held-out fold. The configurations differ between folds; that
  is the point.
- **optimism** = flat_max − nested.

Every (configuration, fold) block is checkpointed, because no single run finishes in one call.

**A note on cost, and a 8× speedup found on the way.** The pipeline scored one molecule per
`OLMPASSofMOL.exe` launch, and the bottleneck was Wine process startup, not scoring. MultiPASS
also ships `OLMPASSofSDF.exe`, which takes a whole SDF in one launch. Switching to it cut 40
molecules from ~40 s to **4.9 s**, and the values are identical: 40/40 molecules agreed exactly
on a trained SMO base, and a full re-run of one stored block reproduced all **93** of its Pa
values with **0** differences. Without this the eleven-target sweep would not have been possible
at all; base *training* (120 bases per target) is now the floor.

## Result

| target | n | positives | configs trainable | grid min | median | max | spread | nested | **optimism** | stage |
|---|---|---|---|---|---|---|---|---|---|---|
| ABL1 (ENST00000318560) | 709 | 101 | 6/6 | 0.803 | 0.817 | **0.844** | 0.041 | 0.824 | **0.019** | full |
| ALK | 1128 | 24 | 6/6 | 0.784 | 0.912 | **0.929** | 0.145 | 0.928 | **0.001** | full |
| KIT | 959 | 24 | 6/6 | 0.774 | 0.847 | **0.898** | 0.124 | 0.892 | **0.006** | full |
| SMO | 476 | 18 | 6/6 | 0.679 | 0.740 | **0.857** | 0.178 | 0.821 | **0.036** | full |
| NT5C2 | 142 | 18 | 6/6 | 0.563 | 0.624 | **0.710** | 0.147 | 0.576 | **0.134** | full |
| EGFR | 1221 | 17 | 6/6 | 0.814 | 0.827 | **0.923** | 0.110 | — | — | flat only |
| MET | 870 | 8 | 6/6 | 0.885 | 0.896 | **0.934** | 0.049 | — | — | flat only |
| ESR1 (ENST00000206249) | 396 | 7 | 6/6 | **0.497** | 0.678 | **0.861** | **0.364** | — | not computable | — |
| IDH1 (ENST00000415913) | 233 | 6 | 4/6 | 0.499 | 0.500 | 0.500 | 0.002 | — | not computable | — |
| FLT3 | 817 | 6 | **0/6** | — | — | — | — | — | — | nothing trainable |
| BTK (ENST00000308731) | 348 | 5 | **0/6** | — | — | — | — | — | — | nothing trainable |

### The five measured values

| target | positives | flat max | nested | optimism |
|---|---|---|---|---|
| ALK | 24 | 0.9288 | 0.9277 | **0.0010** |
| KIT | 24 | 0.8977 | 0.8917 | **0.0060** |
| ABL1 | 101 | 0.8436 | 0.8242 | **0.0194** |
| SMO | 18 | 0.8572 | 0.8209 | **0.0363** |
| NT5C2 | 18 | 0.7099 | 0.5762 | **0.1337** |

The measured range is **0.001 to 0.134**, a factor of 130, and nothing about a target's size
anticipates where in that range it falls. NT5C2 loses 0.134, which takes it from 0.710 to
**0.576** — from "modest" to "barely above chance". ALK, with the same six-configuration grid and
a *wider* spread than NT5C2 (0.145 against 0.147), loses 0.001. Two targets at 18 positives differ
by a factor of nearly four. It has to be measured per target, which is the whole point.

**A claim in an earlier version of this report was wrong and is withdrawn.** It said "where the
picks are stable the optimism is small", and that the stability of the configuration choice
explains NT5C2's larger loss. With five targets measured it does not hold: SMO's choice is the
*most* stable of the five — only 2 distinct configurations across the 5 outer folds — and its
optimism is 0.036, while ALK chose 3 different configurations and lost 0.001. Over the five
targets the number of distinct configurations chosen has Spearman +0.22 with optimism (exact
p = 1.00). What distinguishes ALK from SMO is not stability but **what the choices were worth**:
ALK's three picks are worth 0.9275, 0.9288 and 0.9275 by the outer CV — a span of 0.001, so the
choice cannot hurt — whereas SMO's two picks are worth 0.857 and 0.803, a span of 0.054.

| target | distinct configs chosen | true AUC of the picks | span | optimism |
|---|---|---|---|---|
| ALK | 3 | 0.9275–0.9288 | 0.001 | 0.0010 |
| KIT | 3 | 0.8647–0.8977 | 0.033 | 0.0060 |
| ABL1 | 3 | 0.8123–0.8436 | 0.031 | 0.0194 |
| SMO | 2 | 0.8029–0.8572 | 0.054 | 0.0363 |
| NT5C2 | 4 | 0.5706–0.7099 | 0.139 | 0.1337 |

### Does anything predict the optimism? Honestly, no

Seven candidate statistics were tested against optimism over the five measured targets, with
**exact permutation p-values** — scipy's asymptotic p is badly wrong at n = 5 (it returns 0.037
for a statistic whose exact p is 0.083). At n = 5 the smallest attainable two-sided p is
**0.0167**, so this analysis can only ever be suggestive.

| statistic | ρ | exact p | status |
|---|---|---|---|
| flat max − mean true AUC of the picks | +1.000 | 0.0167 | **tautological** — this *is* optimism, by construction |
| span of the picks' true AUCs | +0.900 | 0.0833 | near-tautological, and not significant |
| mean AUC over the whole grid | **−1.000** | 0.0167 | independent of the choice; perfect rank order |
| number of resistant variants | −0.632 | 0.2667 | not significant |
| grid spread (best − worst) | +0.500 | 0.4500 | not significant |
| spread among the top three configurations | +0.500 | 0.4500 | not significant |
| number of distinct configurations chosen | +0.224 | 1.0000 | not significant |

Two of the three "significant" rows are circular and must not be reported as findings: the first
is optimism rewritten, and the second is one step downstream of the same choice. The only
non-circular candidate is the **mean AUC over the grid**, where the rank order is perfect and
inverse — targets whose configurations all work well have little selection optimism — but it
reaches exactly the floor p of 0.0167 because a perfect order is the *only* way to get any
significance at n = 5, and seven statistics were tested, so a Holm correction would demand
p < 0.0024, which is unattainable here. **Nothing survives correction.** The defensible statement
for the manuscript is the measured range and the per-target values, not a rule.

This also withdraws a second thing said earlier in this project: that optimism grows as the
number of positives falls. Over the five measured targets that correlation is −0.63 with an exact
p of 0.27, and ALK (24 positives, 0.001) against NT5C2 (18 positives, 0.134) is the counterexample.

### The spread, which is the cheap half of the argument

For EGFR and MET the nested stage has not finished — only the flat
grid is in (30 blocks instead of 120); ALK and KIT have since completed and are in the table
above. That still answers a question a reviewer will ask
immediately: how much does the answer depend on which configuration you pick?

On **six** configurations the spread between best and worst is 0.041 for ABL1 (101 positives) and
0.110–0.178 for every target with 17–24 positives. ESR1, with 7 positives, spans **0.497 to
0.861** — three of its six configurations are at chance and three are near 0.86, and the
published number is the maximum of that. Table 1's real grid is ~674 configurations, so these
spreads are **lower bounds** on how much room the selection had.

MET is the apparent exception, spread 0.049 on 8 positives; with that few positives a narrow
spread is not evidence of stability, and MET's nested value was not computed.

### Where the question cannot be answered

- **ESR1 and IDH1**: nested selection is not computable. MultiPASS declines to train on the inner
  training folds — `B.HST` reports every activity type "Failed", 0 selected, average IAP 0.0000 —
  so no configuration survives the inner comparison and there is nothing to choose between.
  IDH1 additionally has 2 of 6 configurations untrainable even in the outer CV, and its four
  trainable ones all sit at 0.499–0.500, i.e. at chance.
- **FLT3 and BTK**: **no configuration is trainable at all** under position-grouped 5-fold CV.
  With 6 and 5 resistant variants, grouping by residue position leaves training folds with 0–1
  positives and MultiPASS refuses outright, six refusals out of six in each case. These two
  targets cannot be evaluated by this protocol, let alone corrected.

A cheap probe (`nested_probe.py`) was used to predict computability before paying for full runs.
It trains bases on the inner folds of outer fold 0 only, with no scoring. It was right about BTK
and FLT3 and wrong about ESR1 and MET: ESR1's outer fold 0 trains fine and outer fold **2** is
the one that fails, so the probe is a necessary but not sufficient condition and is reported as
such. MET passed the probe via MNA 4 only — MNA 9 declined on an inner fold with 6 positives,
which shows the refusal depends on the MNA level and not only on the positive count.

## What this means for the manuscript

1. **Table 1 needs a second column.** For the three targets where it is measurable, report the
   nested value beside the maximum: ABL1 0.844 → 0.824, SMO 0.857 → 0.821, NT5C2 0.710 → 0.576.
   A sentence of the form "the configuration was chosen by maximising cross-validated accuracy;
   re-selecting inside each training fold lowers it by 0.019–0.134" converts the single biggest
   methodological objection into a disclosed, quantified limitation.
2. **Say plainly that four rows cannot be corrected, and two cannot be evaluated.** ESR1, IDH1,
   FLT3 and BTK have 5–7 resistant variants each. That is the honest reading of what happens to a
   target-specific model at that sample size, and it is far better to state it than to let a
   reviewer discover that a 0.995 in Table S2 rests on five variants.
3. **The grid spread belongs in the supplement for all eleven targets.** It costs nothing to
   present, it is measured here, and it pre-empts the "why this configuration?" question.
4. **These optimism figures are lower bounds**, on six configurations against the real ~674, with
   one outer seed. Do not present them as the full correction.

## Scope of verification

`verify_new.py`, section 8 — the whole file now runs **179 checks, 0 failures**.

**Verified.** That all eleven targets have a result file, and that each one's n and n_pos match
the variant table. That trainable plus untrainable configurations equals six for every target.
That each target's grid minimum, median and maximum in `nested_summary.csv` match its own JSON,
and that `flat_best_config` really is the argmax of the grid. That optimism equals flat_max minus
nested for all three measured targets, that every configuration the nested procedure picked is in
the grid, and that all five outer folds picked one. That BTK and FLT3 have zero trainable
configurations, six listed as untrainable, and exactly six recorded refusals in their logs. That
the target with the most positives has the narrowest grid spread. That exactly three targets have
a measured optimism and that those three are ABL1, SMO and NT5C2 with values 0.0194, 0.0363 and
0.1337. Separately, that the batched scoring route reproduces the per-molecule route exactly: 40
molecules on a trained base, and all 93 Pa values of a stored block, with 0 differences. The
nested AUC for ABL1, NT5C2 and SMO is also recomputed from the per-fold checkpoints, with the
five held-out folds checked to be disjoint and to cover the common variant set.

**Verified, as of 2026-10-17.** Three of the sentences above describe the verifier as it stood on
2026-10-10, not as it stands now: "all three measured targets", "exactly three targets" and "those
three are ABL1, SMO and NT5C2 with values 0.0194, 0.0363 and 0.1337". `nested_summary.csv` was
maintained by hand and was not updated when ALK, EGFR and KIT completed -- it still recorded them
as "flat only" -- so those three checks were asserting a stale file rather than catching it. The
table is now derived from the per-target checkpoints by `multipass/make_nested_summary.py`, and
`verify_new.py` checks that the committed file is byte-for-byte what the checkpoints produce. It
also checks that six targets carry a measured optimism -- ABL1 0.0194, ALK 0.0010, EGFR 0.0019,
KIT 0.0060, SMO 0.0363, NT5C2 0.1337 -- that the three with no inner-fold fit are MET, ESR1 and
IDH1, and that the two where nothing trains at all are BTK and FLT3. The identity "optimism equals
flat_max minus nested" now runs on all six rather than three, with a tolerance of 2e-4: that is
the size of the double-rounding artefact in ALK, where 0.9287500 - 0.9277462 = 0.0010038 rounds to
0.0010 while the two printed columns give 0.9288 - 0.9277 = 0.0011. The tolerance is no larger
than the artefact, so a real disagreement in the fourth decimal still fails. `verify_new.py` makes
191 assertions, up from 179; the no-download suite is 824, up from 812.

**Not verified** *(as of 2026-10-10; EGFR and MET were completed on 2026-10-16, see the
addendum)*. The nested value for ALK, KIT, EGFR and MET — only the flat grid was run, so
their optimism is unknown and the table says so rather than guessing it. Whether the optimism on
the real ~674-configuration grid is larger than on this six-point grid; it must be, since a
maximum over more correlated estimates is more biased, but by how much is not measured. Any
dependence on the outer fold assignment: one seed, because the main table already reports
fold-assignment variance as a standard deviation. Why MultiPASS refuses to train on some folds —
the refusal is observed and recorded from `B.HST`, not explained; it is a closed binary. Whether
a different grid (other L values, MNA levels, or the charge flag that Table 1 also varies) would
give different optimism; charge was not varied here at all.

## Files

| file | contents |
|---|---|
| `nested_summary.csv` | the eleven-target table above |
| `nested_<target>.json` | full runs: grid, nested picks per fold, nested AUC, optimism |
| `nested_<target>_flat.json` | flat-only runs: grid and its spread |
| `nested_probe.json` | the trainability probe for the five small-positive targets |
| `nested_config.py` | the nested-selection pipeline, with batched scoring and `FLAT_ONLY` |
| `nested_probe.py` | the cheap computability probe |
| `nested_run.sh`, `nested_flat.sh` | resume-driven drivers |
| `../nested_checkpoints.tar.gz` | every (configuration, fold) checkpoint, so the four unfinished targets can resume |


---

# Addendum, 2026-10-16: EGFR and MET completed, and the limit of the procedure

EGFR finished all 120 blocks; MET returned a verdict rather than a number. With those two the
nested stage is closed: there is nothing further to compute on this grid.

## The six targets where the nested estimate exists

| target | resistant variants | grid best | nested estimate | optimism | grid spread |
|---|---|---|---|---|---|
| ALK | 24 | 0.9288 | 0.9277 | **0.0010** | 0.1452 |
| EGFR | 17 | 0.9233 | 0.9214 | **0.0019** | 0.1097 |
| KIT | 24 | 0.8977 | 0.8917 | **0.0060** | 0.1236 |
| ABL1 (ENST00000318560) | 101 | 0.8436 | 0.8242 | **0.0194** | 0.0411 |
| SMO | 18 | 0.8572 | 0.8209 | **0.0363** | 0.1781 |
| NT5C2 | 18 | 0.7099 | 0.5762 | **0.1337** | 0.1472 |

EGFR picked `L33_MNA9` on three outer folds and `L25_MNA9` on two, with inner AUCs of 0.84–0.89
throughout, so an imperfect choice costs almost nothing — the same mechanism that keeps ALK at
0.0010.

**ALK's row does not subtract, and 0.0010 is right.** 0.9288 − 0.9277 = 0.0011, but both AUC
columns are rounded independently and for ALK the two moves go in opposite directions
(0.928750 up to 0.9288, 0.927746 down to 0.9277, together 0.0000962). Recomputed from the
per-fold checkpoints at full precision the AUCs are 0.9287500 and 0.9277462 and the optimism is
**0.0010038**. The `optimism` field in
`nested_ALK.json` is computed before rounding and is the authoritative one; the earlier
verification tolerated this because it compared the difference with a 6e-4 tolerance, which is
wider than the 1e-4 artefact.

## The dependence on dataset size is gone

| correlation with optimism | five targets (10 Oct) | six targets (16 Oct) |
|---|---|---|
| number of resistant variants | -0.63, exact p = 0.27 | **-0.12, exact p = 0.84** |
| grid spread | +0.50, exact p = 0.45 | **+0.43, exact p = 0.42** |

Exact two-sided permutation p over all 720 orderings; the floor at n = 6 is 0.0014.

EGFR is the single observation that does it: **17 resistant variants and optimism 0.0019, against
NT5C2's 18 and 0.1337** — a factor of seventy at effectively the same dataset size. Any account
of the optimism that runs through sample size is now refuted by a matched pair, not merely
unsupported by a weak correlation. This is the third claim about *what predicts* the optimism to
fail in this report; the honest position is that the only thing that predicts it is what the
candidate configurations are worth on the outer folds, which is itself only measurable after the
fact.

## Three targets where the nested procedure is not computable

Not "not finished" — not computable. On each of these, **every one of the six candidate
configurations declines to train on the inner folds**, so no outer fold can make a selection:

| target | resistant variants | negatives | configurations failing on the inner folds | flat grid |
|---|---|---|---|---|
| MET | 8 | 855 | 6 of 6 | trains at all six, 0.885–0.934 |
| ESR1 (ENST00000206249) | 7 | 387 | 6 of 6 | trains at all six, but three sit at ~0.50 |
| IDH1 (ENST00000415913) | 6 | 226 | 6 of 6 | two do not train at all; the four that do are all ~0.50 |

In each case `B.HST` reports every activity type as Failed, 0 selected, average IAP 0.0000. The
cause is arithmetic: a three-fold split of an outer training set holding seven or eight positives
leaves as few as **two** resistant variants in an inner training fold, and MultiPASS will not build
a base on two. MET's log records it directly — "declined to train on 459 variants (2 resistant)".

The flat grid trains normally on all three, so this is a limit of the nested *procedure*, not of
the models.

**BTK and FLT3 fail one step earlier.** At 5 and 6 resistant variants, not one of the six
configurations trains even on the *flat* grid — `flat_auc` is empty and `flat_not_trainable`
lists all six for both — so they have no grid best, and therefore nothing for the nested stage to
be optimistic about. An earlier draft of this addendum called them "never attempted", which was
wrong: they were attempted and nothing trained.

**So the procedure has a floor.** It runs at 17 resistant variants and above, and fails at 8 and
below; nothing in our panel sits between 9 and 16, so the floor is bracketed but not located. The
consequence for the paper is worth stating plainly: the targets whose published accuracy is most
likely to be inflated by selection are exactly the ones on which the inflation cannot be measured.

## IDH1 deserves a separate note

IDH1's flat grid is at chance — 0.5, 0.4985, 0.5, 0.5 on the four configurations that train.
Table 1 reports 0.816 for IDH1 with ivosidenib under LOOCV. These are different quantities (per
target-drug best-of-grid LOOCV against drug-pooled position-grouped CV on a six-point grid) and
with 6 resistant variants neither is worth much, but the gap is large enough that the manuscript's
existing statement that these small-n models are provisional is doing real work and should not be
softened.

## Verified

`verify_optimism.py` (formerly `verify_optimism5.py`) re-derives every number above from the
`nested_*.json` summaries and the per-fold checkpoints: the six optimism values, the three
not-computable verdicts with their failure counts, both Spearman correlations with exact
permutation p-values, and ALK's full-precision 0.0010038 recomputed from its checkpoints.
