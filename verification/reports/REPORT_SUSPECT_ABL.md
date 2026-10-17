# MultiPASS vs SUSPECT-ABL — comparison against a resistance-specific method

Generated 2026-10-08. This is the comparison the manuscript's own plan asks for
(`SAV-DR_handover.md` §5.3, tier 3) and the first one in this project against a method that
predicts **drug resistance** rather than pathogenicity.

Everything previously in `level2/REPORT_LEVEL2.md` — AlphaMissense, REVEL, CADD, PolyPhen-2,
SIFT, ESM-2 650M and 150M — predicts **pathogenicity**. Those seven answer a different question
and should be presented only as the handover's tier 2 ("asymmetric in our favour, zero-shot,
lower bound"), never as the comparison that answers the supervisor.

## The comparator

SUSPECT-ABL (Zhou, Portelli, Pat, Nguyen, Rodrigues, Pires & Ascher;
`biosig.lab.uq.edu.au/suspect_abl`) predicts the resistance profile and binding free-energy
change ΔΔG of ABL1 mutations against eight inhibitors, using graph-based structural signatures.
Its own paper reports MCC up to 0.73 and Pearson up to 0.77.

Predictions were obtained through the server's mutation-list endpoint in six submissions of ≤50
mutations (its documented cap), selecting the four inhibitors it shares with our labels:
**imatinib, dasatinib, nilotinib, bosutinib**. Raw output in `suspect_abl_predictions.csv` —
one binary call and one ΔΔG per (variant, inhibitor).

## Comparability — checked, not assumed

- **Sequence.** SUSPECT-ABL's reference is UniProt P00519, 1130 aa, and is **byte-identical** to
  our `ABL1_ENST00000318560`: CRC-85FE6C1C0E483EA2 for both. Numbering needs no offset. (For
  contrast, our other ABL1 FASTA — isoform b, 1149 aa — agrees with our own variant table's
  wild-type residues only 7.23% of the time, so this check is not a formality.)
- **Window.** The server accepts positions **242–493** only. The comparison set is our ABL1
  variants in that window: 280 of 709, holding 97 of our 101 positives.
- **Task.** Resistant to ≥1 of the four shared inhibitors (87) vs no resistance annotation at
  all (183). The 10 variants annotated resistant only to the class "tyrosine kinase inhibitor
  - NS" are **excluded**, since their status for the four is unknown. Final set: 270 variants,
  160 positions.
- **Sanity.** T315I, the canonical pan-TKI gatekeeper, is called Resistant by SUSPECT-ABL for
  all four inhibitors (ΔΔG 2.32 / 5.08 / 3.75 / 2.45 kcal/mol). The server behaves as published.

**The comparison is biased in SUSPECT-ABL's favour.** It was trained on COSMIC-derived
resistance data that very likely overlaps COSMIC v100, the source of our labels, and it cannot
be retrained. Our figures are out-of-fold under position-grouped CV. A tie therefore favours
us, and a SUSPECT-ABL win is not evidence of a better method. This is the handover's tier-3
caveat and it must appear next to the table.

## Result

| model | AUC | 95% CI | AUPRC | DeLong vs MultiPASS | Holm |
|---|---|---|---|---|---|
| full_gbm | **0.618** | 0.538–0.691 | 0.451 | 0.082 | 0.49 |
| **SUSPECT-ABL (max ΔΔG)** | **0.595** | 0.520–0.662 | 0.407 | 0.452 | 1.00 |
| nearest_pos | 0.580 | 0.509–0.657 | 0.380 | 0.486 | 1.00 |
| SUSPECT-ABL (mean ΔΔG) | 0.562 | 0.493–0.633 | 0.372 | 0.921 | 1.00 |
| **MultiPASS** | **0.557** | 0.481–0.630 | 0.375 | — | — |
| pos_only | 0.552 | 0.480–0.623 | 0.338 | 0.889 | 1.00 |
| subst_lr | 0.546 | 0.469–0.625 | 0.432 | 0.845 | 1.00 |

87 positives / 183 negatives. **Nothing separates from anything**: every CI spans 0.48–0.70,
no comparison survives Holm, and the whole table sits between 0.546 and 0.618.

Per drug, imatinib (82 positives, the only shared drug with enough data). The base is the
270-variant set, i.e. after the 10 TKI-NS-only positives are already removed; of its 87
positives, 82 are imatinib-labelled and the other **5 are excluded** as resistant to one of the
remaining three drugs but not imatinib. Negatives stay at 183. (On the 280-variant base the
same exclusion reads 15, being those 5 plus the 10 TKI-NS-only — stated because this ambiguity
produced the one failure in the verification pass.)

| score | AUC | 95% CI | AUPRC |
|---|---|---|---|
| MultiPASS Pa(Imatinib) | 0.571 | 0.499–0.644 | 0.428 |
| SUSPECT-ABL ΔΔG(imatinib) | 0.560 | 0.489–0.639 | 0.350 |
| MultiPASS Pa(Resistant) | 0.556 | 0.485–0.625 | 0.366 |

Paired DeLong, MultiPASS Pa(Imatinib) vs SUSPECT-ABL: **p = 0.824**. Dead even.

### Why everything is near chance here, and why that is the point

The 242–493 window is essentially the span of our own positives (244–486), so this comparison
sits in the same regime as the `region` setting of the main analysis — the position shortcut is
gone. Our own margin sweep gives MultiPASS 0.589 at ±5 and 0.599 at ±10 on ABL1, against 0.860
on the whole protein; the 0.557 measured here, on a slightly different variant set (270 vs 311)
and label (87 four-drug positives vs 101 any-drug), is consistent with that.

So the finding is not "MultiPASS is weak and SUSPECT-ABL is strong". It is that **a dedicated,
structure-based, resistance-specific predictor also lands near chance on these labels in this
window** — while enjoying a training-data advantage. That relocates the difficulty from our
method to the task and the labels, which is a considerably stronger position for the manuscript
than the pathogenicity comparison ever was.

### What SUSPECT-ABL's own classifier says about our labels

Using its binary call rather than ΔΔG, on the 270-variant set:

| | value |
|---|---|
| TP / FP / FN / TN | 14 / 17 / 73 / 166 |
| sensitivity | 0.161 |
| specificity | 0.907 |
| precision | 0.452 |
| MCC | **0.100** |

Against its published MCC of up to 0.73 on its own data, 0.100 on ours. And drug-matched — for
each of our positives, did SUSPECT-ABL call Resistant for the *same* drug it is labelled for:

| drug | matched recall |
|---|---|
| imatinib | **5 / 82 = 0.061** |
| dasatinib | 5 / 10 = 0.500 |
| nilotinib | 4 / 9 = 0.444 |
| bosutinib | 2 / 5 = 0.400 |
| **total** | **16 / 106 = 0.151** |

SUSPECT-ABL calls Resistant for 11.4% of the 280 variants; our label prevalence in the same
window is 34.6%.

**The imatinib row is the one to look at.** An independent resistance-specific method, trained
on overlapping data, agrees with 6% of our imatinib-resistance labels. Two readings, and the
paper should state both rather than pick one: either SUSPECT-ABL is far too conservative for
imatinib, or a large share of our 82 imatinib positives are passengers — variants observed in a
resistant patient and curated as resistance mutations without being causal. The second reading
is the mirror image of the caveat already in `REPORT.md` about negatives being "no annotation"
rather than "tested and found non-resistant": the positives may be inflated the same way. It is
also testable, by the saturation prime-editing dataset discussed earlier.

## Scope of verification

**Verified.** Byte-identity of SUSPECT-ABL's reference sequence with ours (CRC64 both sides).
That the 50-mutation cap is the server's own documented limit, by reading its error
("Invalid number of mutations :70. Please provide less than or equal to 50"). That all 280
variants returned a ΔΔG and a class for all four inhibitors. The one-to-one merge of 280
variants between the two sources (asserted). T315I called Resistant for all four. Consistency
of MultiPASS's 0.557 with its own margin sweep at comparable window width. Drug-matched recall
recomputed per drug.

**Not verified.** Whether SUSPECT-ABL's training set actually overlaps COSMIC v100 — this is
inferred from its paper describing COSMIC-derived training data, not established, and it cannot
be checked without their training list. The contamination direction is therefore stated as
likely, not certain. Dr. Kinase, the other tier-3 comparator the handover names, has not been
run. The four inhibitors were matched by name; no check was made that SUSPECT-ABL's structural
model for each uses the same co-crystal our drug labels implicitly refer to. Positions outside
242–493 (429 of our 709 ABL1 variants, 4 positives) are untested by this comparison.

## Files

| file | contents |
|---|---|
| `suspect_abl_predictions.csv` | raw SUSPECT-ABL output: ΔΔG and class per variant × 4 inhibitors |
| `suspect_abl_compare.csv` | the AUC table with CIs, AUPRC, DeLong and Holm |
| `suspect_compare.py` | the comparison script |
