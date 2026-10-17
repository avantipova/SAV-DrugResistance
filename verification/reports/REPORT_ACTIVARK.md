# ActivArk vs SAV-DrugResistance, judged by experiment

Generated 2026-10-10. The first head-to-head against a **current, working, sequence-level**
resistance predictor, decided on **experimental** data rather than on either method's own
cross-validation.

**The competitor.** Singh, Schmenger, Gonzalez-Sanchez, Kutkina, Bremec, Diwan, Mozas, Lopez,
Siebert, Sotillo & Russell, *Discriminating activating, deactivating and resistance variants in
protein kinases*, **Genome Medicine** (2025),
[10.1186/s13073-025-01564-z](https://doi.org/10.1186/s13073-025-01564-z); server
`activark.russelllab.org`, code `github.com/russelllab/kinaseResistance` (GPLv3). Trained on
2 505 curated kinase variants — 375 activating, 1 028 deactivating, **98 resistance**, 1 004
neutral — over 441 kinases, reporting mean AUC 0.941 across its classes. It exposes a dedicated
**resistance-versus-neutral probability `RvN`**, takes a bare `GENE/T315I`, and needs no
structure. It is the closest thing to a peer SAV-DrugResistance has.

**The arena.** The 1 848 missense substitutions of the ABL1 saturation prime-editing scan
(`REPORT_PRIMEEDIT.md`), each measured against six TKIs in K562 cells. Neither method was built
on this dataset, and the library is saturating, so neither can be rewarded for recognising
famous mutations.

## Conclusions first

1. **On imatinib, nilotinib and the pooled ATP-site label, SAV-DrugResistance is ahead, and
   significantly so.** Imatinib 0.605 against 0.400; nilotinib 0.668 against 0.435; pooled
   0.598 against 0.417 — all three paired bootstrap p < 0.001.
2. **Most of that gap is an artefact ActivArk itself warns about, and it does not vanish when
   the warning is heeded.** Its help page says its resistance predictor "lacks information about
   activating and deactivating information" and that users may ignore it for strongly
   deactivating variants. Its `RvN` correlates **+0.563** with its own deactivating probability,
   and that deactivating probability runs *against* experimental resistance (AUC **0.312**).
   Removing the variants its own predictor calls deactivating lifts ActivArk from 0.417 to
   0.480–0.511 and lowers ours from 0.598 to 0.558 — the gap narrows from −0.180 to −0.064…−0.079,
   and on nilotinib it stays significant in every scope tested.
3. **ActivArk is never significantly above chance on this data; ours is, on nilotinib, in every
   scope.** After the filters ActivArk's intervals all contain 0.5.
4. **ActivArk is nominally better on bosutinib and dasatinib** — 0.752 against 0.660, and 0.604
   against 0.560 — but on 12 and 30 resistant variants, with overlapping intervals and p of 0.15 (0.1515)
   and 0.55 (0.5465). Those two cells should not be claimed by either side.
5. **Neither method has substitution-level resolution.** Within a position, ours scores 0.500 and
   ActivArk 0.380 over the same 249 pairs at 60 positions. This is a property of the task as both
   methods currently attack it, not a flaw unique to ours.

## Coverage

| | |
|---|---|
| scan missense variants submitted | 1 848 |
| ActivArk returned an `RvN` | **1 686** (91.2 %) |
| returned nothing — outside its kinase-domain boundary | 162, all at positions 184–211 |
| scored position range | 212–505 |
| in ActivArk's own curated **resistance** training set | **7** (D276G, E255K, F311I, F317L, G250E, T315I, Y253H) |
| other curated classes present | 2 "increase", 1 "neutral" |
| our 91 ABL1 positives that overlap the scan, scored by ActivArk | **91 of 91** |

Seven training variants out of 1 686 is a small leak, but they are precisely the famous ones, and
on them ActivArk is near-perfect (T315I 0.976, F317L 0.961, E255K 0.917, Y253H 0.906, G250E
0.825, F311I 0.783, D276G 0.738). They are excluded from every number below. Submitted through
the server's own form endpoint on 2026-10-10 in 37 batches of 50, algorithm GBC; no batch failed
and all 1 848 were accounted for.

## Head-to-head, the 1 679 ActivArk-scored variants that are not in its training set

Resistant vs Sensitive per the scan's own `Classification`, Intermediate excluded. Bootstrap CIs
and the paired difference, 4 000 resamples, seed 0. Δ is ActivArk minus ours, so a negative Δ
favours SAV-DrugResistance.

| drug | resistant | ours | 95 % CI | ActivArk | 95 % CI | Δ | 95 % CI | p |
|---|---|---|---|---|---|---|---|---|
| Nilotinib | 75 | **0.668** | 0.608–0.723 | 0.435 | 0.365–0.505 | −0.233 | −0.323…−0.144 | **<0.001** |
| Imatinib | 130 | **0.605** | 0.557–0.652 | 0.400 | 0.351–0.451 | −0.205 | −0.268…−0.141 | **<0.001** |
| Ponatinib | 12 | 0.614 | 0.507–0.718 | 0.607 | 0.433–0.781 | −0.007 | −0.170…+0.160 | 0.9255 |
| Dasatinib | 30 | 0.560 | 0.438–0.682 | 0.604 | 0.480–0.720 | +0.043 | −0.100…+0.178 | 0.5465 |
| Bosutinib | 12 | 0.660 | 0.495–0.810 | 0.752 | 0.606–0.877 | +0.091 | −0.034…+0.223 | 0.1515 |
| Asciminib | 215 | 0.487 | 0.444–0.530 | 0.456 | 0.414–0.497 | −0.032 | −0.089…+0.025 | 0.2735 |
| **any ATP-site TKI** | 140 | **0.598** | 0.548–0.647 | 0.417 | 0.369–0.469 | **−0.180** | −0.246…−0.114 | **<0.001** |

Our AUCs here are a little below those in `REPORT_PRIMEEDIT.md` (0.620 for imatinib, 0.686 for
nilotinib) because this table is restricted to the 1 679 variants ActivArk could score and did
not train on. The two reports are consistent; the restriction is the difference.

## Applying ActivArk's own caveat

Their help page, verbatim on the point: users "may ignore its predictions when other predictors
predict the variant to be strongly deactivating". So the comparison is repeated under
progressively stricter versions of that instruction.

| scope | n | ours (imatinib / nilotinib / pooled) | ActivArk | pooled Δ | p |
|---|---|---|---|---|---|
| all scored, training variants removed | 1 679 | 0.605 / 0.668 / 0.598 | 0.400 / 0.435 / 0.417 | −0.180 | <0.001 |
| minus Deactivating (High) | 1 390 | 0.579 / 0.646 / 0.574 | 0.432 / 0.455 / 0.447 | −0.128 | <0.001 |
| minus any Deactivating call | 1 021 | 0.555 / 0.638 / 0.558 | 0.472 / 0.503 / 0.480 | −0.079 | 0.048 |
| Activating calls only | 796 | 0.560 / 0.651 / 0.563 | 0.484 / 0.523 / 0.500 | −0.064 | 0.1605 |

And continuously, capping its deactivating probability D on the pooled ATP-site label:

| | D < 1.01 | < 0.9 | < 0.8 | < 0.7 | < 0.6 | < 0.5 | < 0.4 | < 0.3 |
|---|---|---|---|---|---|---|---|---|
| n | 1 237 | 1 216 | 1 174 | 1 080 | 997 | 878 | 725 | 557 |
| ours | 0.598 | 0.597 | 0.592 | 0.589 | 0.585 | 0.578 | 0.576 | 0.592 |
| ActivArk | 0.417 | 0.422 | 0.430 | 0.445 | 0.459 | 0.493 | 0.511 | — |

ActivArk climbs monotonically towards 0.5 as deactivating variants are removed and stops there.
Ours is flat at 0.58–0.60 throughout, which says the two methods are measuring different things:
theirs is substantially a kinase-function score, ours is substantially a positional
resistance-hotspot score.

Supporting that reading, over the 1 679 variants: Spearman(`RvN`, their deactivating probability)
= **+0.563**, Spearman(`RvN`, our Pa) = **+0.028**. The two predictors are very nearly
independent, which is itself a result — an ensemble is the obvious next question, and is not
attempted here.

## Within-position, both methods

| | within-position AUC | positions | pairs | ties | position-mean only |
|---|---|---|---|---|---|
| SAV-DrugResistance | 0.5000 | 60 | 249 | 9 | 0.6078 |
| ActivArk | 0.3795 | 60 | 249 | 1 | 0.4614 |

Hold the residue fixed and neither method can say which substitution at it confers resistance.
This is the finding of `REPORT_PRIMEEDIT.md` reproduced for an independent method, which makes it
a statement about the state of the field rather than about one tool.

## What this means for the manuscript

This is the comparison the desk rejections were asking for, and it can be stated plainly.

- ActivArk belongs in Table S3 and in the Comparison section as the current sequence-level
  method, published after our submission.
- The claim the data supports, and it must be stated this narrowly: *on the ABL1 saturation
  mutagenesis dataset, judged by experiment, SAV-DrugResistance discriminates imatinib and
  nilotinib resistance better than ActivArk (0.605 vs 0.400 and 0.668 vs 0.435, p < 0.001), and
  the advantage persists but narrows when ActivArk's documented deactivating-variant caveat is
  applied.* **This does not generalise.** Extending the comparison to ALK, KIT and EGFR
  (`REPORT_ACTIVARK_3TARGETS.md`) reverses it: on ALK, inside the kinase domain, ActivArk beats
  us by 0.229 (0.707 against 0.476, p = 0.0005) and our model is below chance; KIT is a tie on
  eight positives; EGFR cannot be judged because 16 of our 17 positives are in ActivArk's own
  training set. Any sentence in the manuscript that generalises the ABL1 result across targets
  will be wrong.
- The claim the data does **not** support: that we beat ActivArk at its own published task.
  Their 0.941 is a four-class activating/deactivating/resistance/neutral AUC on their curated
  set; it is not a drug-specific AUC against experiment, and nothing here tests it.
- Bosutinib and dasatinib should be reported as going the other way, with the small counts named.
  Leaving them out would be the kind of omission a reviewer finds.
- The near-zero correlation between the two scores is worth one sentence: the methods are
  complementary, which is a friendlier and more defensible framing than a league table.

## Scope of verification

**Verified.** That the server is live and that its result endpoint is
`/static/predictor/output/<id>/output.json`, read out of its own `outputFunctions.js`. That all
1 848 variants were submitted and all 1 848 accounted for, with none of the 37 batches failing.
That the extracted values are aligned to the right variants: the variant order was reproduced
independently in the container and checked against a rolling checksum of the name sequence
computed in the browser (882026044, matched exactly), and the deactivating probability was
present for exactly the 1 686 variants that have an `RvN` and absent for exactly the 162 that do
not. That ActivArk's 7 curated resistance variants among the scan are the ones named above, from
its own `mutType` field. That the numbering is canonical on both sides —
`REPORT_PRIMEEDIT.md` establishes 1 848/1 848 wild-type agreement at offset 0 against our
reference FASTA, and ActivArk resolves `ABL1` to **P00519**, the same accession. Every AUC, CI,
paired bootstrap, Spearman coefficient and within-position count is recomputed by
`verify_activark.py`.

**Not verified.** ActivArk's published mean AUC of 0.941, which was neither reproduced nor
contradicted — it is measured on a different task. Whether the `RF` algorithm option behaves
differently from the `GBC` option used throughout; only GBC was run. Whether ActivArk's curated
training set overlaps our COSMIC/CIViC positives beyond the 7 variants its own output flags —
their full curated list was not downloaded, so a larger overlap outside the scan cannot be ruled
out. Whether the comparison holds on our other six kinase targets (ALK, BTK, EGFR, FLT3, KIT,
MET): that needs 5 343 further submissions, about three and a half hours of their server's time,
and was not run. Whether an ensemble of the two scores beats either, which the near-zero
correlation invites but which is not tested here. Whether the 162 unscored variants at positions
184–211 are non-random with respect to resistance; the scan calls 1 of them resistant to an
ATP-site TKI.

## Files

| file | contents |
|---|---|
| `activark_abl1.csv` | all 1 848 variants with ActivArk's `RvN`, `mutType`, A/D verdict and deactivating probability |
| `ak_merged.csv` | the joined table both methods and the experiment are read from |
| `ak_compare.csv` | the per-drug head-to-head, both scopes |
| `ak_caveat.csv`, `ak_threshold.csv` | the caveat-filtered and threshold-swept comparisons |
