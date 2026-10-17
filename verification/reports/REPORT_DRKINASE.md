# Dr. Kinase — the second resistance-specific comparator, and why it cannot be benchmarked as a ranker

Generated 2026-10-08. Companion to `REPORT_SUSPECT_ABL.md`. Dr. Kinase is the other tier-3
comparator named in `SAV-DR_handover.md` §5.3.

**Conclusion first: Dr. Kinase cannot be evaluated as a per-variant ranker on our data through
its web server.** It returned a score at the queried residue for **9 of our 280** ABL1
kinase-domain variants (3.2%), covering **6 of 97** positives. Scoring the other 271 as zero
would leave 96.8% of the set tied at one value, so any AUC computed from it would be a
statement about the tie-handling convention and not about the method. No AUC is reported here,
deliberately.

## The tool

Dr. Kinase (Nucleic Acids Research 2025, 53:W258; `modinfor.com/drkinase`;
`github.com/BioDataStudy/DrKinase`) predicts **the loci of four drug-resistance hotspot types**
in protein kinases — gatekeeper, G-loop, αC-helix and A-loop. Its own framing is locus
prediction, not per-(variant, drug) resistance, and it never names a drug in its output. All
seven of our kinases are in its catalogue, under canonical UniProt accessions: ABL1 P00519,
ALK Q9UM73, BTK Q06187, EGFR P00533, FLT3 P36888, KIT P10721, MET P08581.

How the predictions were obtained: its submission form (`submit_job.php`, mutation mode, input
format `GENE MUTATION` per line) plus the `kinaseid[]` field, which the backend needs even in
mutation mode — without it the job fails server-side with
`Error: local variable 'outdir_pred' referenced before assignment`. Results come from the
page's own AJAX endpoint (`get_result.php`, `type=result`, `job_id=<id>`). Our 280 ABL1
variants in positions 242–493 — the same set used for SUSPECT-ABL — went in as 7 jobs of 40.
All 280 came back; per-job payloads were 2.4–3.4 MB.

## What it returns

Per submitted variant: a top-loci list for each of the four hotspot types, each entry
`[position, position, score, sequence context]`, plus one full per-residue score track
(`gk_scores`, length 1130 for ABL1 — matching our sequence length) for the gatekeeper type
only. The other three types have no full track.

The queried position appears in those lists only when it ranks among the method's top loci. For
most of our variants it does not: for example G251C, a classic P-loop position, gets no score at
251 — the nearest nominated G-loop locus is **249 at 0.9928** — and its gatekeeper track is all
zeros. Absence is therefore "not among this kinase's nominated hotspots", not "score 0".

**A correction to my own first reading.** I initially concluded from EGFR T790M and T790A
returning identical scores (0.8983 both) that the method ignores the substituted residue and is
purely position-level. That inference was wrong. At ABL1 position 315 the three substitutions
score **0.5557 (T315N), 0.5495 (T315I), 0.8997 (T315A)**, and at 351 they score **0.1506
(M351V), 0.1534 (M351T), 0.9016 (M351K)**. Dr. Kinase is substitution-aware; the EGFR tie was a
coincidence, not evidence. The earlier claim is withdrawn.

## The nine variants it scored

| variant | pos | Dr. Kinase (gatekeeper) | our label | drugs we annotate |
|---|---|---|---|---|
| M351K | 351 | **0.9016** | positive | Imatinib |
| T315A | 315 | 0.8997 | positive | Dasatinib; TKI-NS |
| T315N | 315 | 0.5557 | positive | Imatinib |
| **T315I** | 315 | **0.5495** | positive | Dasatinib; Imatinib; Nilotinib; TKI-NS |
| M351T | 351 | 0.1534 | positive | Bosutinib; Imatinib; Nilotinib; TKI-NS |
| M351V | 351 | 0.1506 | negative | — |
| T319I | 319 | 0.0715 | negative | — |
| L324P | 324 | 0.0132 | negative | — |
| L324Q | 324 | 0.0131 | positive | TKI-NS |

Six of the nine are our positives, so the loci it nominates are broadly the right ones. The
ordering within them is another matter: **T315I — the single most consequential
imatinib-resistance mutation in ABL1 — scores 0.5495**, below the 0.90 reference line that the
server's own result plot draws, while M351V, which we do not label, scores 0.1506 against
M351T's 0.1534 for the same position. Only M351K clears 0.90.

That the method nominates position 315 and 351 at all is a real signal. That it ranks T315I
mid-pack and leaves 97% of the kinase domain unscored is why it is not a ranker.

## What this means for the manuscript

The two resistance-specific comparators the manuscript names behave very differently, and
neither supports the comparison the Comparison section currently makes:

- **SUSPECT-ABL** is a genuine per-variant predictor and was benchmarked properly
  (`REPORT_SUSPECT_ABL.md`): AUC 0.595 against MultiPASS's 0.557 on the same 270 variants,
  p = 0.452, nothing significant, everything near chance.
- **Dr. Kinase** answers a different question — where the hotspots are — and gives no score for
  97% of our variants. It can be cited as evidence that the field treats kinase resistance
  largely as a locus-level problem, but it cannot be put in an AUC table next to MultiPASS.

The sentence flagged in the handover — "our ABL1-specific SSPR model achieved a significantly
higher AUC of 0.923" against 0.84 for SUSPECT-ABL — remains unsupportable, and now for a second
reason: on our own data, measured head to head, MultiPASS does **not** come out above
SUSPECT-ABL.

## Scope of verification

**Verified.** That `kinaseid[]` is required in mutation mode (the failing job's `error_pred.txt`
read directly). That all 280 variants returned results across 7 jobs. That `gk_scores` has
length 1130, matching our ABL1 sequence. That the scores differ between substitutions at the
same position (315 and 351, six values). That the 0.90 line is drawn by the server's own result
plot (read from its JavaScript) — this is a plotted reference line, and I did not confirm it is
the classification threshold used in the paper. The 9 scored variants were cross-checked against
our label table.

**Not verified.** Whether the full per-residue tracks for G-loop, αC-helix and A-loop can be
obtained at all — only the gatekeeper track is returned by the web endpoint, so three of the
four channels are only visible as top-loci lists. Whether running the GitHub implementation
locally would yield complete four-channel per-residue scores, which is the only route to a real
ranking comparison; that was not attempted. Whether Dr. Kinase's training data overlaps COSMIC
v100 — the handover assumes it does, and I did not confirm it. Nothing was run for ALK, EGFR,
KIT, BTK, FLT3 or MET, although all six are in its catalogue.

## Files

| file | contents |
|---|---|
| `drkinase_abl1.csv` | our 280 ABL1 variants in 242–493 with the Dr. Kinase gatekeeper score where one was returned |
