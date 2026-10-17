# CIViC as an independent label source

CIViC release **01-Oct-2026**, `AcceptedClinicalEvidenceSummaries.tsv`, CC0.
Generated 2026-10-04, **recomputed after the FASTA fix described in `../verification/VERIFICATION.md`**. Scripts: `civic_overlap.py`, `civic_eval.py`. The 90 substitutions and the 70/14/6 split are unchanged by that fix; the AUCs below moved slightly.

> **Numbers corrected 2026-10-08.** The evaluation tables in this file were written on 04-Oct
> at 22:09, but `civic_external_eval.csv` was regenerated on 05-Oct at 16:53 — after the KIT and
> SMO FASTA re-run and after EGFR's fifth fold-assignment seed finished. The report therefore
> carried 12 stale per-target cells, both stale pooled tables and a stale negative-pool size
> (4425/4301 instead of 4407/4283). Every table below is now read straight from that CSV, and
> the conclusions that moved are flagged in place. The overlap table, the discordance list and
> the therapy-matched analysis were *not* affected — those were re-verified line by line against
> `civic_overlap.log` and `therapy_match.log` and match exactly.

## How CIViC was read

A CIViC resistance call = `evidence_type = Predictive` **and** `significance = Resistance`
**and** `evidence_direction = Supports`. The direction filter matters: of the 1050 Resistance
rows, **220 are "Does Not Support"**, i.e. evidence *against* resistance. Counting those would
have inflated the set by a fifth.

Fusion partners are normalised to the kinase (`BCR::ABL1` → `ABL1`, `EML4::ALK` → `ALK`, …)
and compound molecular profiles are mined for every `<gene> <X###Y>` token, so
`Fusion AND ABL1 T315I` is not lost — that pattern alone holds most of the BCR::ABL1 evidence.

**Numbering check:** all 90 CIViC substitutions were tested against the reference FASTA the
pipeline itself reads. **0 wild-type mismatches.** CIViC uses the same numbering we do, for all
six targets, including the EGFR extracellular-domain calls (V441F, G465R, K467T, S492R) which
are sometimes published on mature-protein numbering. No offset correction is needed anywhere.

## Agreement with the COSMIC-derived labels

| gene | CIViC SAVs | our label = positive | our label = negative | absent from our set | also carries "Does Not Support" |
|---|---|---|---|---|---|
| ABL1 | 45 | 39 | 1 | 5 | 33 |
| ALK | 13 | 11 | 2 | 0 | 0 |
| EGFR | 11 | 3 | 7 | 1 | 0 |
| KIT | 14 | 11 | 3 | 0 | 5 |
| NT5C2 | 5 | 5 | 0 | 0 | 5 |
| SMO | 2 | 1 | 1 | 0 | 0 |
| **total** | **90** | **70** | **14** | **6** | **43** |

The 14 discordant calls, with the therapy CIViC attaches:

    ABL1   p.D363G   C      Imatinib
    ALK    p.F1245C  C,D    Crizotinib
    ALK    p.R1275Q  D      ALK Inhibitor TAE684
    EGFR   p.G465R   C      Panitumumab, Cetuximab
    EGFR   p.K467T   D      Cetuximab
    EGFR   p.S492R   C,D    Cetuximab
    EGFR   p.S768I   C      Erlotinib
    EGFR   p.N826S   D      Erlotinib
    EGFR   p.L838V   D      Gefitinib, Erlotinib
    EGFR   p.L858R   C      Crizotinib
    KIT    p.V560D   C      Sunitinib
    KIT    p.L576P   C,D    Sorafenib, Imatinib, Nilotinib; Sunitinib
    KIT    p.D816V   B,D    Dasatinib, Imatinib; Imatinib
    SMO    p.L412F   C,D    Vismodegib

### Therapy-matched: the disagreement is mostly real, not a coverage gap

Our labels are per (variant, drug) — the `Drug` column of each `*_mutations.csv` — so each
CIViC call can be tested against the drugs we actually label for that target
(`therapy_match.csv`, script `therapy_match.py`):

    ABL1   bosutinib, dasatinib, imatinib, nilotinib, tyrosine kinase inhibitor - NS
    ALK    alectinib, brigatinib, ceritinib, crizotinib, lorlatinib
    EGFR   afatinib, erlotinib, gefitinib, osimertinib, tyrosine kinase inhibitor - NS
    KIT    imatinib, nilotinib, sunitinib
    NT5C2  purine   (CIViC names its members, mercaptopurine and thioguanine)
    SMO    vismodegib

Of the 14 discordant calls, **9 name a drug we do label** and only 5 fall outside our coverage:

| | variant | CIViC drug | verdict |
|---|---|---|---|
| **genuine disagreement** | ABL1 D363G | imatinib | we label imatinib |
| | ALK F1245C | crizotinib | we label crizotinib |
| | EGFR S768I | erlotinib | we label erlotinib |
| | EGFR N826S | erlotinib | we label erlotinib |
| | EGFR L838V | erlotinib, gefitinib | we label both |
| | KIT V560D | sunitinib | we label sunitinib |
| | KIT L576P | imatinib, nilotinib, sunitinib | we label all three |
| | KIT D816V | imatinib | we label imatinib |
| | SMO L412F | vismodegib | we label vismodegib |
| outside our coverage | EGFR G465R, K467T, S492R | cetuximab, panitumumab | **antibodies**; our EGFR labels are TKI-only |
| | ALK R1275Q | TAE684 | experimental compound, not in our set |
| | EGFR L858R | crizotinib | an ALK inhibitor against EGFR — cross-context, likely CIViC noise |

This is a sharper finding than a coverage artifact. **Nine well-characterised resistance
mutations sit in our data as negatives for drugs we do label**, including two textbook ones:
`KIT D816V` (imatinib) and `SMO L412F` (vismodegib, the canonical smoothened-inhibitor
resistance mutation). For EGFR specifically the picture splits: 3 of its 7 discordances are
antibody resistance we never set out to cover, 1 is cross-context noise, and **3 are real
misses in our own drug class** (erlotinib/gefitinib).

On the agreeing side the drug sets line up too: of the 70 CIViC calls our labels also call
positive, **67 share at least one drug**. The three that do not are ABL1 F311I, F317R and
M343T, where CIViC cites axitinib — a VEGFR inhibitor outside our panel.

**CIViC is not clean ground truth.** 43 of the 90 substitutions carry both supporting and
contradicting evidence — 33 of ABL1's 45, and all 5 of NT5C2's. Any use of CIViC as a gold
standard has to say so.

## External-label evaluation

No model is refitted here. The scores are the out-of-fold scores already computed for the
`full` setting, rank-averaged over the 5 fold-assignment seeds exactly as in the main table,
so nothing from a test fold entered its own score.

One shared negative pool: variants in our set with COSMIC label 0 **and** no CIViC resistance
evidence. Two positive sets:

- **A — all CIViC:** every CIViC resistance SAV present in our set (84).
- **B — discordant only:** just the 14 our COSMIC labels call negative. These were *negatives
  during training*, so ranking them highly is generalisation, not label recall.

### A — all CIViC positives (pooled: 84 pos / 4407 neg)

| model | AUC | 95% CI | AUPRC | DeLong p | Holm |
|---|---|---|---|---|---|
| full_gbm | **0.828** | 0.780–0.874 | 0.096 | 0.710 | 1.00 |
| **MultiPASS** | **0.821** | 0.774–0.863 | 0.104 | — | — |
| pos_only | 0.799 | 0.746–0.846 | 0.072 | 0.101 | 1.00 |
| subst_lr | 0.619 | 0.548–0.687 | 0.042 | <0.0001 | <0.001 |

Per target (AUC):

| target | pos / neg | MultiPASS | pos_only | full_gbm | subst_lr |
|---|---|---|---|---|---|
| ABL1 | 40 / 610 | 0.856 | 0.826 | **0.896** | 0.722 |
| ALK | 13 / 1098 | 0.906 | **0.940** | 0.903 | 0.550 |
| EGFR | 10 / 1190 | **0.641** | 0.580 | 0.582 | 0.433 |
| KIT | 14 / 933 | **0.886** | 0.843 | 0.788 | 0.614 |
| NT5C2 | 5 / 124 | 0.642 | 0.575 | **0.751** | 0.663 |
| SMO † | 2 / 452 | **0.838** | 0.832 | 0.836 | 0.038 |

† 2 positives — not interpretable, shown for completeness only.

**This reproduces the main finding on an independent label source**, and after the correction it
reproduces it slightly more bluntly than before. Pooled, MultiPASS (0.821) and the position-only
model (0.799) are statistically indistinguishable (p = 0.101 raw, 1.00 after Holm across the 39
tests in this file), and `full_gbm` is now **nominally ahead of MultiPASS** (0.828 vs 0.821,
p = 0.710 — a difference far inside the noise, but the sign has flipped relative to the earlier
numbers in this file). The only model MultiPASS clearly beats is chemistry-only (`subst_lr`,
Holm < 0.001). On ABL1 `full_gbm` is ahead at raw p = 0.0044 (Holm 0.150), and on ALK the
position-only model is ahead outright.

Per target, MultiPASS still has the highest point estimate on EGFR, KIT and SMO, but the KIT and
EGFR margins shrank in the correction (KIT 0.920 → 0.886, EGFR 0.686 → 0.641), and on SMO the
three non-chemistry models are now within 0.006 of one another (0.838 / 0.832 / 0.836) on two
positives — i.e. indistinguishable. Only the KIT comparisons against `full_gbm` (raw 0.0045) and
`subst_lr` (raw 0.0028) remain anywhere near significance, and neither survives Holm (0.150,
0.098).

### B — discordant positives only (pooled: 14 pos / 4283 neg)

| model | AUC | 95% CI | AUPRC | DeLong p | Holm |
|---|---|---|---|---|---|
| **MultiPASS** | **0.626** | 0.484–0.759 | 0.006 | — | — |
| pos_only | 0.596 | 0.412–0.765 | 0.006 | 0.568 | 1.00 |
| full_gbm | 0.579 | 0.410–0.740 | 0.007 | 0.543 | 1.00 |
| subst_lr | 0.566 | 0.379–0.751 | 0.008 | 0.645 | 1.00 |

**This is the one conclusion the correction actually overturns.** The earlier numbers in this
file supported a thin claim that MultiPASS ranks variants which were negatives in its own
training data but are CIViC-confirmed resistant above chance (0.674, CI 0.528–0.798). The
current figure is **0.626 with CI 0.484–0.759 — the interval now includes 0.5**, so the claim is
no longer available: on set B, MultiPASS is not demonstrably better than a coin flip, and no
model separates from any other (all Holm = 1.00). The ordering of the four models is unchanged
and the whole table still sits between 0.566 and 0.626.

The per-target B rows are in the CSV but with 1–7 positives each only the pooled row carries
information. With a single positive (ABL1, SMO) the "AUC" is that one variant's percentile rank,
and the bootstrap CI resamples the same positive every time, which makes it look far tighter
than it is — ABL1's `full_gbm` 0.970 [0.957–0.982] is one variant, not a result. One per-target
comparison does reach Holm significance, **KIT: MultiPASS 0.657 vs pos_only 0.583, raw p =
0.0011, Holm 0.0396** — but it rests on 3 positives and is not offered as a finding.

## Caveats

- **Set A is only partly independent.** 70 of its 84 positives are also COSMIC positives, so A
  is close to the main evaluation with a slightly different positive set. The genuinely
  independent test is B, and B is small.
- **n is small throughout.** 84 positives in A, 14 in B, spread over six proteins. Per-target
  CIs are wide and several per-target comparisons are uninterpretable on their own.
- **Pooling across targets** concatenates within-target rank-averaged scores. Ranks are
  comparable across proteins by construction, but a target with many negatives contributes
  more negatives to the pool.
- **CIViC evidence levels are mostly C and D** (485 D, 316 C, 241 B, 5 A, 3 E across all
  resistance rows) — case-report and preclinical tiers, not practice-guideline tiers.
- **The negative pool is "no evidence", not "evidence of no resistance".** Some of those 4407
  variants are presumably resistant and simply unannotated in both sources, which depresses
  every model's apparent performance equally.
- **Not done:** no therapy-matched analysis. Our labels are drug-specific and CIViC's are too,
  but the evaluation above collapses drug identity. Restricting to matched drug classes would
  remove the cetuximab/panitumumab EGFR mismatch — worth doing if EGFR is a focus.

## Files

| file | contents |
|---|---|
| `civic_resistance_vs_our_labels.csv` | all 90 CIViC SAVs: therapies, evidence levels, whether contradicted elsewhere in CIViC, our label |
| `civic_external_eval.csv` | both label sets × 6 targets + pooled × 4 models: `n_pos`, `n_neg`, AUC, CI, AUPRC, DeLong raw and Holm. There is **no** `interpretable` column — an earlier version of this file claimed one; interpretability has to be judged from `n_pos` (treat rows under 3 positives as uninterpretable) |
| `therapy_match.csv` | every CIViC call against the drugs our labels cover, per variant |
| `civic_overlap.py` | CIViC parsing, numbering check, overlap table |
| `therapy_match.py` | the therapy-matched analysis |
| `civic_eval.py` | the external-label evaluation |
| `verify_civic.py` | 159 checks of this report against the CSV and the two logs; run it after any re-run |
| `civic_overlap.log` | run log |
| `data/civic/01-Oct-2026-*.tsv` | the four downloaded release files |
