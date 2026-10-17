# Seven additional analyses

Generated 2026-10-07, on the post-FASTA-fix data. Companion to `REPORT.md`; nothing here
changes a number in that report's main table. Scripts and CSVs are in this folder.

These close the gaps the first round listed as not done, plus one defect found along the way.
Independent re-verification of everything below: `verify_new.py` / `verify_new.log`,
**86 checks, 0 failures**, with the un-verified scope stated at the end of that log and
repeated at the end of this file.

---

## 1. How large a difference could these comparisons have detected?

Every headline claim in `REPORT.md` is a null — MultiPASS is indistinguishable from a
position-only model, Holm-corrected DeLong p = 1.00 in 11 of 12 cells. That reads as "n is too
small to see anything" unless the sensitivity is quantified. `power.py` → `power.csv`.

For MultiPASS vs `pos_only`, using the DeLong standard error of the difference (which carries
the real correlation between the two models' scores, so it is the paired sensitivity, not two
independent ones):

| target | setting | Δ observed | SE | min. detectable Δ (80% power) | **95% upper bound on Δ** |
|---|---|---|---|---|---|
| ABL1 | full | +0.019 | 0.0105 | 0.029 | **0.036** |
| ABL1 | region | +0.059 | 0.0271 | 0.076 | 0.104 |
| ALK | full | −0.020 | 0.0178 | 0.050 | **0.009** |
| ALK | region | −0.021 | 0.0404 | 0.113 | 0.045 |
| EGFR | full | +0.001 | 0.0195 | 0.055 | **0.033** |
| EGFR | region | −0.024 | 0.0318 | 0.089 | 0.029 |
| KIT | full | +0.035 | 0.0207 | 0.058 | 0.069 |
| KIT | region | +0.008 | 0.0355 | 0.099 | 0.066 |
| NT5C2 | full | +0.036 | 0.0443 | 0.124 | 0.108 |
| NT5C2 | region | +0.038 | 0.0743 | 0.208 | 0.161 |
| SMO | full | +0.038 | 0.0473 | 0.133 | 0.116 |
| SMO | region | +0.051 | 0.0446 | 0.125 | 0.124 |

Median detectable difference 0.094 AUC (range 0.029–0.208). **In none of the 12 cells does the
observed |Δ| reach the detectable threshold**, so the nulls are nulls rather than missed signal.

**The column that should go in the paper is the last one.** On ABL1, ALK and EGFR in the `full`
setting the data *exclude* MultiPASS being better than a position-only model by more than
0.01–0.04 AUC. That is a bounded positive claim, not an absence of evidence. On NT5C2 and SMO,
and throughout the `region` setting, the bound is 0.10–0.16 and the null genuinely is
uninformative — stated as such rather than pooled with the rest.

**One check failed and was replaced.** A binormal Monte-Carlo power simulation gave 0.26–0.98
power at the analytic MDE (mean 0.70, not 0.80). Its type-I rate at Δ = 0 was 0.057, so the
DeLong test is calibrated and the *generator* was wrong: a binormal model cannot reproduce the
tie structure here, where Pa is coarsely quantised and `pos_only` is a step function over
residue positions. It was discarded for a stratified paired bootstrap of the one quantity the
MDE rests on, the SE of the difference, computed on the real scores:
**bootstrap/DeLong SE ratio median 0.984, range 0.958–1.021 across all 36 comparisons**, largest
MDE discrepancy 0.017 AUC. The SE is validated; no power curve is claimed.

---

## 2. `nearest_pos`: one integer, no training, same performance

The sharpest form of the position finding. Score = −(distance in residues to the nearest
resistant position in that fold's training set). No wild-type residue, no substituted residue,
no sequence, no fitted model. `extras.py` → `summary_with_nearestpos.csv`.

| target | setting | MultiPASS | **nearest_pos** | pos_only | full_gbm | subst_lr |
|---|---|---|---|---|---|---|
| ABL1 | full | 0.860 | 0.854 | 0.842 | 0.875 | 0.665 |
| ABL1 | region | 0.589 | 0.562 | 0.529 | 0.610 | 0.576 |
| ALK | full | 0.930 | **0.961** | 0.950 | 0.919 | 0.436 |
| ALK | region | 0.676 | 0.688 | 0.697 | 0.705 | 0.508 |
| EGFR | full | 0.922 | **0.958** | 0.921 | 0.892 | 0.631 |
| EGFR | region | 0.698 | **0.799** | 0.722 | 0.756 | 0.677 |
| KIT | full | 0.903 | 0.893 | 0.869 | 0.831 | 0.700 |
| KIT | region | 0.753 | 0.726 | 0.745 | 0.691 | 0.658 |
| NT5C2 | full | 0.696 | 0.651 | 0.660 | 0.741 | 0.487 |
| NT5C2 | region | 0.645 | 0.608 | 0.606 | 0.694 | 0.505 |
| SMO | full | 0.824 | 0.836 | 0.785 | 0.844 | 0.454 |
| SMO | region | 0.637 | 0.622 | 0.586 | 0.660 | 0.470 |

Mean difference **−0.002 AUC** (i.e. `nearest_pos` marginally ahead on average); MultiPASS ahead
in 7 of 12 cells; **0 of 12 differences survive Holm** (48 tests in that file). On EGFR in the
`region` setting `nearest_pos` is ahead by +0.101.

This upgrades the first report's conclusion from an inference about what `pos_only`'s gradient
booster might be using to a direct reading: *proximity to a known resistance position* accounts
for the performance.

---

## 3. Pa is not a probability — and that is fixable

The tool presents Pa as a probability of activity, so it was tested as one.
`extras.py` → `calibration.csv`, `calibration_reliability.csv`; `recalibrate.py` →
`recalibration.csv`.

**Raw Pa is a worse probability than predicting the base rate, in 10 of 12 cells.** It
systematically overestimates, by up to **5.9×** (ALK `full`: prevalence 2.1%, mean Pa 12.6%;
EGFR `full`: 1.4% vs 6.9%). In ALK's 0.1–0.2 reliability bin sit 581 variants with mean Pa 0.130
against an observed rate of 0.005. ABL1 is the one cell where Pa happens to be calibrated, and
only because its prevalence (14.2%) nearly equals Pa's mean (14.2%).

This is expected for a PASS-style Pa — an invariant-accuracy quantity, not a posterior for this
prevalence — and a monotone recalibration fitted out-of-fold (position-grouped, never fitted on
the variants it scores) repairs it:

| | mean Brier | beats base rate |
|---|---|---|
| raw Pa | 0.0938 | 2 / 12 |
| Platt | **0.0803** | **12 / 12** |
| isotonic | 0.0787 | 10 / 12 |
| base-rate forecast | 0.0832 | — |

**Platt is the recommendation**, not isotonic. An invariance check — a monotone map cannot change
AUC — failed at first, up to 0.098 pooled. Both maps are monotone within every fold (12/12), so
it was not a coding error: the pooled loss comes from stitching five fold-specific maps whose
outputs are not on a common scale, and isotonic additionally flattens variants into ties
(within-fold AUC loss up to **−0.053**, mean **−0.013** across the 12 cells — the −0.053 is the
worst single cell, **NT5C2 `full`**, and an earlier version of this line quoted it without saying
so; the next worst are SMO `full` −0.037 and EGFR `full` −0.030).
Platt's within-fold AUC change is **+0.00008** on average, 0.00012 in absolute value, and exactly
0.00000 in all six `full` cells. So Platt buys
calibration at no ranking cost; isotonic buys slightly better Brier at a real one.

For the web platform this is a concrete change: display a recalibrated value, not raw Pa.

### Operating points

AUC and AUPRC are threshold-free; a tool is used at a threshold. `operating_points.csv`, means
over the six targets in the `full` setting:

| model | precision@top-10 | precision@top-20 | sensitivity @95% spec | @99% spec |
|---|---|---|---|---|
| MultiPASS | 0.467 | 0.400 | 0.465 | 0.253 |
| full_gbm | 0.317 | 0.400 | 0.467 | 0.194 |
| nearest_pos | 0.317 | 0.275 | 0.480 | 0.084 |
| pos_only | 0.283 | 0.308 | 0.413 | 0.128 |
| subst_lr | 0.133 | 0.108 | 0.1205 | 0.015 |

At 95% specificity MultiPASS recovers under half the resistance variants, and `nearest_pos`
slightly more. MultiPASS's one clear edge is in the extreme tail (@99% specificity, 0.253 vs
0.084) — worth stating, since it is the only operating region where it separates.

---

## 4. The region collapse is not an artifact of the ±10 aa margin

`REPORT.md` conceded the margin was a design choice with no sensitivity analysis.
`margin_sweep.py` → `margin_sweep.csv`. The positive set never changes; the margin only admits
more distant negatives (verified).

MultiPASS trained protein-wide, scored inside the margin:

| target | ±5 | ±10 | ±20 | ±30 | ±50 | whole protein |
|---|---|---|---|---|---|---|
| ABL1 | 0.589 | 0.599 | 0.614 | 0.634 | 0.654 | 0.860 |
| ALK | 0.711 | 0.725 | 0.741 | 0.758 | 0.778 | 0.930 |
| EGFR | 0.666 | 0.679 | 0.683 | 0.704 | 0.757 | 0.922 |
| KIT | 0.760 | 0.770 | 0.795 | 0.814 | 0.829 | 0.903 |
| NT5C2 | 0.661 | 0.664 | 0.657 | 0.657 | 0.667 | 0.696 |
| SMO | 0.670 | 0.668 | 0.669 | 0.677 | 0.703 | 0.824 |

Monotone in the margin: Spearman ρ = **+1.00** for ABL1, ALK, EGFR and KIT, +0.83 SMO, +0.52
NT5C2. The ±10 figure is a point on a smooth dose–response curve, not a cliff. ρ between the
number of admitted negatives and AUC is +1.00 in four of six targets.

**The control that makes this an argument rather than an observation:** `subst_lr`, which uses
substitution chemistry and never sees position, has **mean ρ = −0.004** across the margin. A
model blind to position is indifferent to the margin, so the entire margin dependence sits in
the positional channel.

NT5C2 is nearly flat because it has only 124 protein-wide negatives — its "region" is almost the
whole protein, so the margin has little left to admit.

**Side result: retraining inside the region buys nothing.** At ±10, a protein-wide model
restricted to the window beats a model retrained inside the window by **+0.018 AUC on average**
(ahead in 5 of 6 targets). So the `region` setting is not under-reporting the method because of
a smaller training set; there is simply no signal there.

---

## 5. Therapy-matched CIViC, and a caution about the per-drug mode

Closes the "Not done" caveat in `civic/REPORT_CIVIC.md`, whose scored evaluation collapsed drug
identity. Each CIViC call is now scored against the Pa of the drug CIViC actually names.
`civic_therapy_eval.py` → `civic/civic_therapy_eval.csv`.

| target | drug | n+ | Pa(drug) | Pa(binary) | Δ |
|---|---|---|---|---|---|
| ABL1 | Imatinib | 32 | 0.871 | 0.865 | +0.006 |
| ABL1 | Dasatinib | 14 | 0.846 | 0.870 | −0.024 |
| ALK | Crizotinib | 11 | 0.913 | 0.905 | +0.008 |
| KIT | Imatinib | 13 | 0.930 | 0.944 | −0.014 |
| NT5C2 | Purine | 5 | 0.642 | 0.642 | 0.000 |
| ABL1 | Nilotinib | 13 | 0.677 | 0.876 | −0.199 |
| ABL1 | Bosutinib | 5 | 0.720 | 0.884 | −0.164 |
| ALK | Lorlatinib | 5 | 0.818 | 0.972 | −0.154 |
| ALK | Alectinib | 6 | **0.409** | 0.880 | −0.472 |

The split is exactly by training density. For drugs that cleared the ≥10-positive bar of the
per-drug table (first five rows) the drug-specific and pooled scores are indistinguishable,
mean Δ **−0.005** (−0.024…+0.008). For drugs with 1–9 COSMIC positives (nilotinib 9, lorlatinib
7, bosutinib 5, alectinib 5) the drug-specific score **collapses**, mean Δ **−0.247**, and on
ALK/alectinib it falls below chance. One test survives Holm: ALK/lorlatinib, p = 0.0012.

**This is a deployment caution with a clear rule.** In per-drug mode on a sparse drug MultiPASS
does not refuse — it returns a normal-looking Pa that is 0.15–0.47 AUC worse than the pooled
score. Per-drug output should be gated on the number of training positives.

EGFR and SMO are absent for a verifiable reason: only 2 CIViC calls name osimertinib (T790M,
C797S, both our positives) and 2 name vismodegib (L412F our negative, D473H positive) — below
the 5-positive floor. `tyrosine kinase inhibitor - NS` cannot be therapy-matched at all (a class
with no CIViC member list) and was excluded rather than guessed.

---

## 6. Leave-one-target-out: the model does not transfer to an unseen protein

This is the analysis `08_loto.csv` was meant to be. **That file was defective** and is replaced:
written 3 Oct, before the 4 Oct multi-record FASTA fix (KIT at 943 negatives instead of 936, SMO
at 464 instead of 453); its EGFR row had silently fallen back to the `region` variant table via
a `cands` fallback in `loto.py`, giving 17+/251− instead of 17+/1197−; and it contained no
MultiPASS row at all. The fallback is now a hard failure instead of a silent substitution.
`loto_train.py` + `loto_multipass.py` + `loto_summary.py` → `loto.csv`, `loto_baselines.csv`.

One MNA base per held-out protein, trained on every variant of the other five (3409–4481
variants, 178–185 resistant), then every variant of the held-out protein scored. No folds and no
seeds: the split is by protein, so there is nothing random to average. Verified that each base
was trained on exactly the other five targets' variant counts, and that the held-out variant set
is identical to the main table's.

| held out | within-protein | **LOTO, MNA 9** | 95% CI | **LOTO, MNA 4** | pos_only | full_gbm | subst_lr |
|---|---|---|---|---|---|---|---|
| ABL1 | 0.860 | 0.635 | 0.585–0.691 | 0.628 | 0.461 | 0.542 | 0.532 |
| ALK | 0.930 | 0.599 | 0.499–0.696 | 0.577 | 0.572 | 0.530 | 0.602 |
| EGFR | 0.922 | 0.613 | 0.493–0.721 | 0.668 | 0.595 | 0.634 | 0.618 |
| KIT | 0.903 | 0.504 | 0.384–0.603 | 0.619 | 0.366 | 0.305 | 0.609 |
| NT5C2 | 0.696 | 0.573 | 0.459–0.686 | 0.558 | 0.631 | 0.678 | 0.455 |
| SMO | 0.824 | 0.700 | 0.626–0.779 | 0.807 | 0.673 | 0.529 | 0.557 |
| **mean** | **0.856** | **0.604** | | **0.643** | 0.550 | 0.536 | 0.562 |

**Mean drop 0.856 → 0.604, i.e. −0.252 AUC.** Paired DeLong on identical variants is
Holm-significant for ABL1, ALK, EGFR and KIT (p = 6.3e-15, 1.5e-6, 5.6e-9, 3.5e-10). **Four of
six LOTO confidence intervals include 0.5**; KIT's point estimate is 0.504.

Cross-protein, MultiPASS is also indistinguishable from the baselines: 0.604 against
0.550/0.536/0.562, with **1 of 18 tests Holm-significant** — ABL1 against `pos_only` at AUC
0.461. That is not a win: residue numbering is not comparable across proteins, so `pos_only`
cross-protein is a floor, not a competitor.

**MNA level.** A cross-target base needs a single MNA level, which is a free choice, so both 9
and 4 were run for all six targets. No difference is Holm-significant anywhere (largest raw
p = 0.068, KIT; 0.077, SMO; the rest p > 0.26), so the choice does not demonstrably matter and
the honest presentation is the range across both levels. The one apparent exception —
SMO at MNA 4, 0.807 against 0.824 within-protein — has n+ = 18, a CI of 0.722–0.884, and 0.700
at MNA 9, so it reads as spread rather than transfer.

Taken with §2 and §4 the three measurements agree: ~0.92 within a protein, a monotone slide to
~0.6 as the window narrows, ~0.60 on an unseen protein, and an untrained
distance-to-nearest-known-position rule level with the model throughout.

---

## 7. Selection optimism: how much of the published AUC is configuration search?

`REPORT.md` conceded that Table 1 reports the CV value of the configuration chosen by maximising
that same CV. `nested_config.py` re-does the (window length L, MNA level) search *inside* each
outer training fold, so the configuration never sees the variants it is scored on. Grid
L ∈ {17, 25, 33} × MNA ∈ {4, 9}; outer 5 folds, inner 3 folds, both grouped by residue position;
one outer seed, since the question is selection bias and fold-assignment variance is already
reported as sd in the main table. Run on the largest and the smallest target.

| | ABL1 (101+/606−) | NT5C2 (18+/124−) |
|---|---|---|
| flat max over the grid — *what Table 1 reports* | 0.844 (L17/MNA9) | 0.710 (L17/MNA4) |
| median over the grid | 0.817 | 0.624 |
| **nested (configuration chosen honestly)** | **0.824** | **0.576** |
| **selection optimism** | **0.019** | **0.134** |

**Selection optimism scales with how little data there is.** On ABL1 the published configuration
is barely inflated — 0.019 AUC, and the inner CV picks the project's own L17/MNA9 in 3 of 5
folds. On NT5C2 it is 0.134 AUC, and the inner CV picks a *different* configuration in 4 of 5
folds (L25/MNA9, L25/MNA4, L33/MNA9, L17/MNA4, L33/MNA9) — the signature of selection noise
rather than a real optimum.

The implication for the paper is specific: ABL1's number can stand roughly as reported, but the
five small-n targets (ALK 24, EGFR 17, KIT 24, NT5C2 18, SMO 18 positives) sit closer to NT5C2's
situation than ABL1's, so their `full`-setting figures should be read as upper bounds. Measuring
each of them properly is the obvious next run and was not done here.

Two incidental findings. (a) **MultiPASS sometimes refuses to train at all** on small folds:
`B.HST` reports every activity type `Failed`, 0 selected, average IAP 0.0000 — 7 of NT5C2's 90
inner blocks. Recorded as an outcome rather than a crash. (b) The variant set depends on L
(windows running off either protein terminus are dropped), so all comparisons are made on the
intersection across the grid: ABL1 707 of 712 variants, NT5C2 142 of 143. Hence flat L17/MNA9 =
0.844 here against 0.860 in the main table, which pools 5 seeds over 712 variants.

---

## A defect of my own, found by the verification pass and fixed

While training the LOTO bases I launched a duplicate training loop without noticing the first was
still alive, and two processes wrote into the same SMO directory. I killed both, cleaned up and
retrained — but the `rm` that would have removed `train.sdf` was itself interrupted, and
`loto_train.py` only writes the SDF when absent. So the SMO MNA 9 base was retrained from a
**truncated** SDF the killed writer had left: 4093 of 4152 molecules, the missing 59 being the
tail of NT5C2 and including **8 of the 184 resistant training molecules**.

This was caught by the check that reads the substance count out of MultiPASS's own `B.HST` and
compares it with the intended training-set size — the only check in `verify_new.py` that failed
on the first run (1 of 86). After a clean retrain SMO's LOTO AUC at MNA 9 moved **0.665 → 0.700**,
the LOTO mean 0.598 → 0.604, and the mean drop −0.258 → −0.252. No conclusion changed, but the
earlier SMO figure was wrong and is superseded. The verification now passes 86/86.

---

## Scope of verification

**Verified** (`verify_new.py`, 86 checks): the MDE and exclusion-bound arithmetic; the DeLong SE
against a stratified paired bootstrap on all 36 comparisons; `nearest_pos` recomputed from
scratch with a different implementation for all 12 cells; the base-rate Brier identity
p(1−p); monotonicity of both recalibration maps within every fold; reliability-bin counts
against cell totals; that the margin only ever admits negatives and never changes the positive
set; the margin monotonicity and the position-free control; that every LOTO variant set is
identical to the main table's and every LOTO AUC recomputes; that each LOTO base was trained on
exactly the other five targets' variants, read from MultiPASS's own training log; the nested
optimism arithmetic, that `flat_max` is really the grid maximum, and that the five held-out
folds are disjoint and cover the common variant set.

**Not verified.** MultiPASS's internals — a closed Win32 binary. That MNA 9 and 4 are the right
levels for a cross-protein base: both were run and no difference is Holm-significant, but two
points are not a sweep. The nested configuration search covers two targets on a 6-point grid at
one outer seed; it measures selection optimism at two sample sizes, not for all six targets. No
power curve is claimed anywhere — only the SE of the AUC difference was validated, after the
binormal simulation was found unreliable. The CIViC therapy matching depends on drug-name
normalisation, and `tyrosine kinase inhibitor - NS` was excluded rather than mapped.

## Files

| file | contents |
|---|---|
| `power.csv`, `power_secheck.csv`, `power.py` | detectable effect, exclusion bounds, DeLong-vs-bootstrap SE check |
| `summary_with_nearestpos.csv`, `extras.py` | all five models incl. `nearest_pos`, with DeLong and Holm |
| `operating_points.csv` | precision@top-k and sensitivity at fixed specificity |
| `calibration.csv`, `calibration_reliability.csv` | Brier vs base rate, 10-bin reliability curves |
| `recalibration.csv`, `recalibrate.py` | nested Platt and isotonic recalibration, with the monotonicity diagnostics |
| `margin_sweep.csv`, `margin_sweep.py` | AUC vs resistance-span margin, all models |
| `loto.csv`, `loto_baselines.csv` | leave-one-target-out, both MNA levels, CIs, paired DeLong |
| `loto_train.py`, `loto_multipass.py`, `loto_summary.py` | the LOTO pipeline |
| `nested_ABL1_ENST00000318560.json`, `nested_NT5C2.json`, `nested_config.py` | selection-optimism measurement |
| `civic/civic_therapy_eval.csv`, `civic_therapy_eval.py` | therapy-matched CIViC evaluation |
| `verify_new.py`, `verify_new.log` | the 86 verification checks |
