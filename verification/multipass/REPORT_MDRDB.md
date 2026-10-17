# MdrDB — the first external benchmark that is quantitative and not COSMIC-derived

Generated 2026-10-09. Not a method comparison: MdrDB is a database. It is the first chance in
this project to ask whether our resistance score and our labels track a **continuous, externally
measured quantity** rather than another curated binary annotation.

**Why that matters here.** Every external check so far shares our labels' ancestry. DRSP agreed
with our COSMIC calls on 439 of 440 variants (κ = 0.995), and CIViC overlaps COSMIC heavily. MdrDB
is different: its ΔΔG values come from GDSC and DepMap IC50 ratios, the Platinum affinity
database, AIMMS, KinaseMD and a tyrosine-kinase-inhibitor set. On the mutations we share it
agrees with us at **κ = 0.28** — it is genuinely independent, and it disagrees a great deal.

**Two findings, one good for the paper and one not.**

1. **Our labels are externally validated.** On 168 shared mutations, the ones COSMIC calls
   resistant have measured ΔΔG of mean 0.838 / median 0.680 kcal/mol against mean 0.169 / median
   0.147 for the rest (Mann–Whitney p = 1.1 × 10⁻⁴, permutation p = 0.00015 over 20 000 shuffles).
   This is the first independent, quantitative support for the label set — which has been the
   weakest part of the manuscript, since our negatives mean "no annotation" rather than "tested
   and found non-resistant".
2. **The model adds nothing over the labels on this benchmark.** As a predictor of MdrDB's own
   resistance call, MultiPASS gets AUC **0.716** (CI 0.572–0.844) and **our binary COSMIC label,
   used directly as the predictor, gets 0.712** (CI 0.598–0.823). On the five targets with
   within-protein CV it is 0.710 against 0.688. The score is not measurably better than the
   labels it was trained on. This is the same shape as the `pos_only` result and should be
   reported in the same breath.

## The data, and how it was obtained

`MdrDB_CoreSet_release_v1.0.2022.tsv`, 5 032 052 bytes, sha256 `efa06906…abf2`, taken from the
project's own GitHub repository (`raw.githubusercontent.com/tencent-quantum-lab/MdrDB`). The
website's download page is gated behind a form asking for name, e-mail, telephone and
institution, so the repository copy — the same core set, ungated — is what was used.

Contents as read from the file: **4292 samples, 240 proteins, 410 drugs, 2503 mutations**, every
row carrying a ΔΔG, from seven sources (GDSC 2914, Platinum 792, AIMMS 276, TKI 144, KinaseMD 68,
DepMap 65, RET 33). The 240 proteins and 2503 mutations match the paper's headline; the 440 drugs
it advertises are for the full set, the core set holds 410.

MdrDB's own convention: **ΔΔG > 1.36 kcal/mol is a ten-fold affinity drop and is called
resistant**. Positive ΔΔG means weaker binding, i.e. resistance-ward.

### Our slice

**647 single-substitution rows** on **nine** of our eleven targets (nothing for NT5C2 or IDH1),
over 46 drugs. After mapping, **411 rows / 168 distinct mutations** matched, of which **18 are
resistant** by MdrDB's threshold. The 236 dropped rows are overwhelmingly mutations absent from
our COSMIC v100 missense export for that transcript; 11 failed the numbering check below.

| gene | rows | mutations | drugs | resistant |
|---|---|---|---|---|
| ABL1 | 114 | 45 | 5 | 19 |
| EGFR | 102 | 27 | 13 | 2 |
| MET | 69 | 20 | 6 | 2 |
| ALK | 35 | 20 | 4 | 4 |
| BTK | 33 | 13 | 3 | 0 |
| KIT | 24 | 18 | 8 | 0 |
| SMO | 14 | 10 | 2 | 0 |
| ESR1 | 10 | 5 | 2 | 0 |
| FLT3 | 10 | 10 | 1 | 0 |

## MdrDB mixes two ABL1 numberings inside one file

The wild-type check that every external join in this project goes through — this time against the
**reference sequence** rather than our variant table, so that "numbering disagrees" and "we have
no such variant" stay separate diagnoses — found that MdrDB's ABL1 entries agree with the
1130-residue canonical sequence for 108 of 122 distinct substitutions, and **nine resolve only on
the 1149-residue isoform b** (E216K, G340W, P460L, R153C, R185K, R351W, R502W, V390I, Y276C). The
shift is exactly 19, the length difference between the two isoforms, and each of the nine carries
the right wild-type residue at its own number on isoform b.

So the offset was resolved **per row**, against a short justified candidate list — the isoforms
that actually exist for that gene — and accepted only when the wild-type residue matched.
Allowing an arbitrary per-row offset would let one match in twenty through by chance. After that,
**every one of the 411 matched rows has the correct wild-type residue, 0 mismatches**. Eleven rows
matched neither numbering and were dropped by name (ABL1 A156T, two BTK, eight EGFR including
P1071L and S703F) — probably other isoforms or data errors on their side.

This is worth a sentence in the paper if MdrDB is cited: a database that mixes isoform numbering
within a single protein will silently corrupt any join that assumes one convention.

## Result

Nothing is refitted. Two score panels as elsewhere: within-protein CV (`full` setting, per-target
best configuration, out-of-fold, rank-averaged over 5 seeds) for ABL1/ALK/EGFR/KIT/SMO, and LOTO
(L = 25, MNA 9, held-out protein) for all nine.

| panel | n | Spearman ρ vs ΔΔG | 95% CI | p | resistant | AUC vs ΔΔG > 1.36 | 95% CI |
|---|---|---|---|---|---|---|---|
| pooled, within-protein CV | 120 | **0.317** | 0.139–0.474 | 0.00043 | 17 | **0.710** | 0.591–0.817 |
| pooled, LOTO | 168 | **0.220** | 0.062–0.372 | 0.0041 | 18 | **0.716** | 0.572–0.844 |
| TKI only, LOTO | 30 | 0.360 | 0.011–0.651 | 0.051 | 11 | 0.675 | 0.462–0.880 |
| AIMMS only, LOTO | 23 | 0.455 | 0.069–0.716 | 0.029 | 8 | 0.637 | 0.391–0.859 |
| GDSC only, LOTO | 122 | 0.109 | −0.071–0.281 | 0.23 | 6 | 0.671 | 0.472–0.843 |
| TKI only, within-protein CV | 30 | 0.143 | −0.211–0.427 | 0.45 | 11 | 0.474 | 0.268–0.694 |
| AIMMS only, within-protein CV | 23 | 0.137 | −0.361–0.561 | 0.53 | 8 | 0.433 | 0.208–0.650 |
| GDSC only, within-protein CV | 74 | 0.213 | −0.017–0.416 | 0.069 | 5 | 0.722 | 0.522–0.930 |

The correlation is positive and significant pooled, in both panels. By source it behaves as it
should for LOTO: the two measurement-based subsets (TKI 0.360, AIMMS 0.455) correlate better than
GDSC (0.109), whose ΔΔG is derived from averaged IC50 values and has a median of −0.10 — i.e.
mostly noise around zero.

**One thing that does not behave and is not explained.** On TKI and AIMMS the within-protein CV
scores correlate *worse* than the held-out-protein LOTO scores (0.143 vs 0.360, 0.137 vs 0.455),
which is backwards. On n = 30 and n = 23 with CIs that span zero and overlap each other, this is
within noise, but it is stated rather than smoothed over, and it means neither per-source figure
should be quoted on its own.

### Per gene, LOTO

| gene | n | ρ | 95% CI | p |
|---|---|---|---|---|
| FLT3 | 10 | 0.517 | −0.148–0.962 | 0.13 |
| ABL1 | 45 | 0.297 | −0.024–0.557 | 0.048 |
| EGFR | 27 | 0.253 | −0.145–0.577 | 0.20 |
| BTK | 13 | 0.187 | −0.466–0.747 | 0.54 |
| MET | 20 | 0.152 | −0.410–0.634 | 0.52 |
| SMO | 10 | 0.079 | −0.494–0.589 | 0.83 |
| ALK | 20 | −0.051 | −0.483–0.386 | 0.83 |
| KIT | 18 | −0.130 | −0.590–0.313 | 0.61 |

Only ABL1 reaches nominal significance and only just; every interval except FLT3's includes zero.
Eight genes at n = 10–45 cannot carry per-gene conclusions, and none is offered.

### Agreement with our labels

| | our label 0 | our label 1 |
|---|---|---|
| **MdrDB ΔΔG ≤ 1.36** | 122 | 28 |
| **MdrDB ΔΔG > 1.36** | 7 | 11 |

Raw agreement 0.792, **Cohen κ = 0.280**. Of our 39 positives here, 28 do not clear MdrDB's
threshold; of MdrDB's 18 resistant calls, 7 are our negatives. Both directions are informative:
a curated resistance annotation and a ten-fold affinity drop are not the same event, and a
variant can confer clinical resistance without a large measured ΔΔG (and the reverse).

## What this means for the manuscript

- **Use it for the labels, not for the model.** The defensible new sentence is that COSMIC-derived
  resistance annotations are supported by independently measured binding-affinity changes
  (p ≈ 10⁻⁴ on 168 mutations from a source agreeing with us only at κ = 0.28). That answers a
  reviewer who asks whether the positives mean anything.
- **Do not claim the model beats the label.** 0.716 against 0.712 is the whole story, and it
  belongs next to the `pos_only` result rather than in a separate optimistic paragraph.
- **MdrDB is also the obvious route to a much larger external test** if the project wants one:
  the full set is 100 537 samples, and the core set already covers nine of our eleven targets. The
  limit here is our own variant tables, which lack 225 of its 647 single substitutions on those
  genes.

## Scope of verification

`verify_mdrdb.py`, **58 checks, all passing**.

**Verified.** The file's byte size against the server's, with sha256 recorded; its 4292 rows, 167
columns, 240 proteins, 410 drugs, 2503 mutations, seven sources, and that every row carries a
ΔΔG. That the slice is 647 single-substitution rows over exactly nine of our targets, and 411
matched. That every matched row's wild-type residue equals the reference sequence's — 0
mismatches — and that the only offsets applied are 0 and +19, the +19 ones ABL1 only, nine
distinct substitutions, each independently confirmed to carry the right residue at its own number
on isoform b, with the isoform lengths (1130, 1149) and their difference checked. That the
collapse to one row per mutation takes the maximum ΔΔG and that the binary call is exactly
ΔΔG > 1.36. Both Spearman coefficients recomputed as Pearson correlations on ranks rather than
through `scipy.spearmanr`, and both AUCs from the Mann–Whitney U statistic rather than from
`sklearn`. The cross-tabulation, raw agreement and κ. The label-versus-score AUC comparison, and
that the score's margin over the label is at most 0.01. The ΔΔG group sizes, means, medians and
Mann–Whitney statistic, plus an assumption-free permutation test over 20 000 shuffles. That the
TKI and AIMMS subsets correlate better than GDSC, and that GDSC's ΔΔG is centred near zero.

**Not verified.** How MdrDB computed each ΔΔG — the GDSC and DepMap values are IC50-derived by
their pipeline and were taken as given; no attempt was made to re-derive one. Whether the 225
mutations absent from our tables are missing because they are not in COSMIC v100 or because of a
transcript choice on our side — only the numbering was checked, not the ingest. Whether MdrDB's
mixed ABL1 numbering extends to other proteins beyond the eleven rows that failed outright; nine
proteins were checked, not 240. The full set (19 GB of structures, 100 537 samples) was not
downloaded, so these figures are from the core set only. Whether the per-source differences are
real or sampling noise — the CIs say they cannot be separated at these n. And nothing here tests
any other method: MdrDB has no predictor to compare against.

## Files

| file | contents |
|---|---|
| `mdrdb_overlap.csv` | the 411 matched rows: MdrDB's mutation and ours, the offset applied, source, drug, ΔΔG, both labels, both scores |
| `mdrdb_compare.csv` | the score panels with Spearman, CIs, AUC and AUPRC |
| `mdrdb_pergene.csv` | the per-gene LOTO panel |
| `mdrdb_compare.py` | the analysis, including per-row offset resolution |
| `verify_mdrdb.py` | the 58 checks |
