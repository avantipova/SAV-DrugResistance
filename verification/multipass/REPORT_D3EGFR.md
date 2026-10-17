# D3EGFR — the server is up, the prediction engine is not

Generated 2026-10-09. Fourth resistance-specific comparator from the manuscript's own Table S3,
after SUSPECT-ABL, Dr. Kinase and DRSP.

D3EGFR (Shi, Y. *et al.*, *Brief Bioinform* **25**, bbae121, 2024; `www.d3pharma.com/D3EGFR/`)
is the one comparator in Table S3 built exactly like SUSPECT-ABL but for a target of ours: deep
learning, **single mutation in, drug response out**, EGFR, eight approved drugs. On paper it was
the best remaining candidate for a second genuine head-to-head.

**Conclusion first: it cannot be benchmarked, because its prediction module does not return
results for mutations outside its database, and all seventeen of our EGFR resistance positives
are outside it.** Of 50 EGFR substitutions submitted, **1 returned a prediction**. Of our 17
positives, **0** did. No AUC is reported here, deliberately — the same decision as for Dr. Kinase
and for the same reason.

## What the server does, and what it does not

The prediction module is a plain GET form: `resultAi.php?task_mutation=<SAV>`, one substitution
per query, documented at about 10 seconds per new mutation. For a substitution the server already
holds, the page returns a clean table — this is **L858R**, the canonical sensitising mutation:

| drug | predicted response | P(CR/PR) | P(SD) | P(PD) | score |
|---|---|---|---|---|---|
| Dacomitinib | CR/PR | 92.0 % | 8.0 % | 0.0 % | −2.12 |
| Furmonertinib | CR/PR | 90.4 % | 9.5 % | 0.0 % | −1.99 |
| Almonertinib | CR/PR | 90.2 % | 9.8 % | 0.0 % | −1.97 |
| Osimertinib | CR/PR | 90.2 % | 9.8 % | 0.0 % | −1.97 |
| Afatinib | CR/PR | 87.7 % | 12.2 % | 0.1 % | −1.80 |
| Icotinib | CR/PR | 82.8 % | 17.1 % | 0.2 % | −1.53 |
| Gefitinib | CR/PR | 78.1 % | 21.6 % | 0.3 % | −1.33 |
| Erlotinib | CR/PR | 77.8 % | 21.8 % | 0.4 % | −1.32 |

Four of those eight — erlotinib, gefitinib, afatinib, osimertinib — are drugs our EGFR labels
cover, so the task would have lined up well. The resistance score would be P(PD), or −score.

For any other substitution the page instead says *"Sorry, the mutation is not in our database, we
will predict it now. Please wait seconds to refresh this page"*, with a ten-second countdown that
reloads the same URL. **That state never resolves.** T790M was re-requested at 0, 15, 30 and
76 seconds and again several minutes later; it stayed pending every time. The behaviour is
identical for G719A, G719C, G719S, S768I, C797S, D761Y, L747P, L792F, E746K, A750P, T854A,
L833V and V769L — in other words for every canonical EGFR resistance mutation tried.

To rule out a session or cookie effect, the T790M result page was then opened as a **real
navigation** rather than a background fetch, and reloaded four times at nine-second intervals,
exactly as the server's own countdown script does. All four reloads returned the pending page.

## The hit rate

Submitted: our **17 EGFR resistance positives** in the kinase domain (positions 718–797) plus a
seeded random sample of **33 negatives** from the same window, 50 in all.

| | submitted | returned a prediction |
|---|---|---|
| our EGFR positives | 17 | **0** |
| sampled negatives | 33 | **1** (L858R) |
| **first probe, total** | **50** | **1 (2 %)** |

A second probe of 15 well-known EGFR substitutions (6 of them already in the first set) returned
**L861Q** and nothing else — also a sensitising mutation, not a resistance one. Across both
probes: **59 distinct substitutions submitted, 2 answered**, and both answers are sensitising
mutations.

The pattern is consistent with D3EGFRdb being the patient-derived part of the resource: the two
substitutions that work, L858R and L861Q, are the two classic *sensitising* mutations that
dominate their 1158-patient cohort. The *predictive* half, D3EGFRAI, is what would have to run
for everything else, and it is the half that does not answer. The site is stamped
"Version 1.00, Last update 2022-11-1", which is consistent with a compute backend that has since
been switched off while the static pages stay served.

Scoring the 49 non-returning variants as zero would leave 98 % of the set tied at one value, so
any AUC computed from it would describe the tie-handling convention rather than the method. That
is the same situation as Dr. Kinase, and it gets the same treatment: no number.

## What this means for the manuscript

The Comparison section cites D3EGFR as a current method. That citation is fine, but the
Supplementary Table S3 row claiming "Accuracy ≈ 0.81–0.85" should carry the observation that the
figure cannot be reproduced through the server, because the server will not score a new mutation.
That is a statement about reproducibility, not about the method's quality — their published
accuracy may well be right, and nothing here tests it.

It also sharpens the manuscript's own argument. Of the five published web servers Table S3 lists
besides ours, **two are now unusable**: KRDS returns an empty response (`bcbl.kaist.ac.kr`,
published 2018) and D3EGFR's prediction module never completes. A third, Dr. Kinase, answers a
different question (hotspot loci, no per-variant score for 97 % of a kinase domain). That is a
concrete, checkable claim for a web-server paper to make, and it is a stronger argument for
SAV-DrugResistance existing than any accuracy comparison: the field's tools stop working.

## Scope of verification

**Verified.** That the form is a GET to `resultAi.php` with the single field `task_mutation`,
read out of the page's own DOM. That L858R returns a nine-row table with eight drugs, three
probabilities and a score, and that the three probabilities sum to 100 % in every row to within
display rounding (two rows read 99.9 % and 100.1 %). That L861Q returns the same shape. That 48 of the 50 submitted substitutions and 13 of the 15 probe
substitutions return the "not in our database, we will predict it now" page, identified by its
own text. That the pending state persists across re-requests at 0, 15, 30 and 76 seconds and in
a later call. That the 17 positives are exactly our EGFR `y = 1` variants in positions 688–875
and that all 17 of our EGFR positives fall in 718–797, so the kinase-domain window loses none of
them. That the negative sample was drawn with a fixed seed from the same window.

**Not verified.** Whether the prediction backend is permanently down or was down only during
this session — it was probed over roughly twenty-five minutes on one day, from one client. The
result page was opened as a real navigation and reloaded on the server's own schedule, so a
missing cookie or session is ruled out, but the form itself was never submitted by clicking its
button, so a hidden client-side step cannot be excluded absolutely. Whether their database holds
more of our variants under a different spelling — only the documented point-mutation format was
used. D3EGFR's published accuracy of 0.81–0.85, which was neither reproduced nor contradicted.
Whether the other seven drugs it reports would have added anything, since four was already enough
overlap for a comparison.

## Files

| file | contents |
|---|---|
| `d3egfr_probe.csv` | every substitution submitted, with the server's response class |
