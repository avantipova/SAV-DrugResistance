# Resistance substitutions the literature documents and our labels miss

Generated 2026-10-14. **This report replaces `REPORT_PROSPECTIVE.md`, whose framing was wrong.**
That version called this a prospective test of substitutions "reported after the training freeze".
Checking the first publication year of each variant shows that is false for six of the seven: the
date filter in the harvest constrained the *search window*, not when the variant was first
reported. What the experiment actually measures is **label coverage** — whether the models find
resistance substitutions that the literature documents and that COSMIC lists as missense variants
in our targets, but that our COSMIC resistance annotation did not label. That is a narrower claim
and a different one.

> **Corrected 2026-10-18.** The first-reported years for EGFR V802F and ESR1 L536R were wrong
> in the first version of this report and are corrected above. V802F was given as 2009; no V802F
> literature exists before 2017, where it first appears as a minor ctDNA variant in a patient
> progressing on osimertinib (Ou et al., *Lung Cancer* 108, 228-231), and the causal
> demonstration is later still (Tanaka et al., *Cancer Sci.* 117, 2426-2442, 2026). L536R was
> given as 2019; it was already shown to drive estrogen-independent growth in work published
> online in December 2016 (Toy et al., *Cancer Discov.* 7, 277-287), which itself attributes the
> variant to reports from 2013-2014. Which paper first named L536R is unresolved: it is in none
> of the three obvious 2013-2014 abstracts, and the supplementary tables could not be read.
> Neither correction changes the label-coverage conclusion, which is that our resistance
> annotation does not label either variant; it changes only how long each has been known.

## What was wrong, precisely

| variant | first reported | the earlier report implied |
|---|---|---|
| EGFR V802F | **2017** | post-freeze |
| EGFR L844V | **2014** | post-freeze |
| ESR1 L536R | **2016** | post-freeze |
| BTK T474I | 2024 | post-freeze |
| BTK L528W | 2024 | post-freeze |
| FLT3 Y842C | 2024 | post-freeze |
| **ESR1 L370F** | **2026** | correct |

Only **ESR1 L370F** was first reported after COSMIC v100. So there is exactly one prospective data
point, not seven.

A second error, raised as the right question before I caught it myself: **our models are
drug-specific, and three of the seven are resistant to drugs that are not in our panel.** A
variant resistant to pirtobrutinib says nothing about an ibrutinib model unless cross-resistance
is established.

## The audit, per variant

| target | variant | drug in the literature | in our panel? | first reported | evidence | Pa | Pi | call |
|---|---|---|---|---|---|---|---|---|
| ESR1 | **L370F** | endocrine therapy (fulvestrant) | yes | **2026** | **sound** | 0.012 | 0.900 | not resistant |
| ESR1 | **L536R** | endocrine therapy | yes | 2016 | **sound** | 0.523 | 0.007 | **resistant** |
| EGFR | **V802F** | osimertinib | yes | 2017 | **sound** | 0.095 | 0.141 | not resistant |
| BTK | L528W | pirtobrutinib; covalent class by cross-resistance | class only | 2024 | indirect | 0.103 | 0.576 | not resistant |
| BTK | T474I | pirtobrutinib; ibrutinib only as a **compound** mutant | **no** | 2024 | **excluded** | 0.180 | 0.132 | resistant |
| FLT3 | Y842C | investigational agents only | **no** | 2024 | **excluded** | 0.086 | 0.524 | not resistant |
| EGFR | L844V | osimertinib, review list only | yes | 2014 | weak | 0.003 | 0.892 | not resistant |

The three exclusions, with the reason:

- **BTK T474I.** *Sci Rep* 2025 reports ibrutinib resistance for the **compound** T474I–C481S, not
  for T474I alone; the single-substitution evidence is for pirtobrutinib, which is not in our
  panel. Our BTK model is ibrutinib-specific. **This was the only variant at a position the model
  had never seen that it called resistant** — so removing it removes the one apparent generalisation
  success.
- **FLT3 Y842C.** It appears only in panels of resistant cell lines used to test *new* agents —
  SGI-7079, LT-171-861, HDZ-740. No quizartinib or sorafenib claim, and those are the two FLT3
  drugs we model.
- **EGFR L844V.** A review's list of on-target mechanisms, no primary report found, earliest
  mention 2014 in an afatinib review.

**BTK L528W** is kept but marked indirect: *Ann Hematol* 2026 states cross-resistance to covalent
agents as a class, which covers ibrutinib, but not ibrutinib specifically.

## The result, corrected

**On the three with sound evidence and the drug or class in our panel: 1 of 3.** And the one hit,
ESR1 L536R, sits at a position that is **already a training positive** (L536H), so the model had
only to carry a known position to a different mutant residue. At positions it had not seen —
EGFR V802F at distance 5, ESR1 L370F at distance 10 — it called both not resistant, L370F
emphatically (Pa 0.012 against Pi 0.900).

**On the single genuinely new variant, ESR1 L370F: 0 of 1.**

**The control is unchanged and still matters: 33 of 33** training positives recalled by the same
four models, mean Pa 0.52 (BTK), 0.81 (EGFR), 0.92 (ESR1), 0.59 (FLT3). The models are not broken.
They recall what they were trained on and do not find what they were not.

Counting all seven as before gives 2 of 7, but that number should not be used: it rests on
T474I, whose drug is outside our panel.

## What this is still good for

The reframing costs the prospective claim but leaves something the manuscript needs:

- **A label-coverage finding.** EGFR V802F has been in the literature since 2017, in
  osimertinib-progression sequencing, and ESR1 L536R since 2016 as conferring estrogen-independent
  growth. Both are present in our ingested COSMIC missense set — COSMIC lists the substitution —
  but our resistance annotation does not label either. So **the COSMIC resistance annotation is
  incomplete relative to the published literature**, and the models inherit that incompleteness.
  This is a real, checkable limitation and it bounds every performance figure in the paper: the
  negative class contains known resistance variants.
- **One honest prospective data point.** ESR1 L370F, first reported in 2026, with fulvestrant
  resistance shown by ERα binding, transcriptional activity, receptor degradation and
  proliferation — missed, with Pa 0.012 against Pi 0.900. One point is not a validation, and it
  should be reported as one point.
- **A concrete next step rather than a caveat.** A proper prospective test needs a set assembled
  *before* seeing the scores and drawn only from reports postdating the label freeze, with the
  drug matched to our panel per variant. The BTK pirtobrutinib series (*Blood* 2026,
  10.1182/blood.2024027009 — T474I/F/S/Y/L 26 %, L528W 16 %, V416L, A428D, D539A/G/H, Y545N across
  39 patients) would be the natural basis, but it requires adding pirtobrutinib to the panel
  first, which means new labels, not a new evaluation.

## Addendum, 2026-10-14: is there anything first reported after May 2025?

Asked directly, and the answer is **no — not from the literature, not yet.**

The search was rebuilt so that both earlier mistakes were impossible. It is **keyed on our own
panel drugs**, so every candidate is tied from the start to a drug one of our models is trained on
(imatinib/nilotinib/dasatinib/bosutinib/ponatinib for ABL1, the five ALK inhibitors, ibrutinib,
osimertinib and the earlier EGFR TKIs, the endocrine agents, quizartinib and sorafenib,
ivosidenib, crizotinib, thiopurines, vismodegib). And novelty is tested by **counting** papers
that mention the variant with its gene **before 2025-05-01**, not by sorting on publication date —
the sort is unreliable and had reported EGFR L861Q, which has **2 068** earlier papers, as first
seen in 2026.

**44 candidates** came out of the harvest. Papers mentioning each before May 2025:

| | |
|---|---|
| MET M546L | **0** |
| ESR1 L370F | **1** |
| MET G401A | 2 |
| FLT3 F692L | 4 |
| EGFR M766T | 6 |
| EGFR V802F | 9 |
| KIT G680R | 10 |
| … the remaining 37 | 12 to 12 030 |

**Zero candidates qualify.** The single variant with no earlier papers, MET M546L, is an artefact
of my own harvest and is rejected by an existence guard added for exactly this: it has **no papers
at all**, because the context window that was supposed to require the gene name matched the
English word *"met"* rather than the gene. Without that guard the novelty test passes anything
that does not exist.

The closest thing to a genuinely new variant remains **ESR1 L370F**: one incidental earlier
mention (a 2023 endometrial-carcinoma review, where the string does not appear in the abstract at
all and so is almost certainly a table or full-text hit) against three papers in 2026 that
actually characterise it, including the fulvestrant-resistance work. The model scores it
**Pa 0.012 against Pi 0.900** — missed.

**So a prospective test set does not exist yet from this route**, and the manuscript should not
claim one. Two honest options follow. Report the single L370F data point as a single data point,
which is what it is. Or build the prospective set from the clinical sequencing series instead of
from novel-variant reports — the *Blood* 2026 pirtobrutinib cohort (T474I/F/S/Y/L 26 %, L528W
16 %, V416L, A428D, D539A/G/H, Y545N over 39 patients) is a genuine post-freeze clinical dataset,
but its drug is not in our panel, so using it means adding pirtobrutinib and relabelling, not
re-evaluating.

**Scope of this addendum.** Verified: the 44 candidates and their pre-May-2025 counts, the
existence guard and what it rejected, and that no candidate passes both gates. Not verified: a
variant reported only in a supplementary table of a closed-access paper would be missed, since the
search reaches abstracts and open-access full text; and the harvest's own recall was not measured,
so "nothing qualifies" is a statement about what this search found, not a proof that the
literature contains nothing.

## Scope of verification

**Verified.** The first publication year of each of the seven, from the earliest Europe PMC record
mentioning the variant with its gene, sorted ascending — this is what overturned the earlier
framing. The drug attribution for each, from the quoted sentence of the cited paper, and whether
that drug is among the ones our model for that target is trained on. That all seven are single
substitutions whose wild-type residue matches our reference FASTA, and that all seven are present
in our ingested COSMIC missense set but absent from our resistance positives. The Pa and Pi values
and the Pa > Pi calls, read from the executable's own two-row output at each target's Table 1
configuration. The 33-of-33 control. Recomputed by `verify_label_coverage.py`.

**Not verified.** Whether COSMIC v100 itself annotates V802F, L536R or L844V as resistant anywhere
in a field our extraction did not read — the claim here is that *our* labels miss them, not that
COSMIC does; separating those two needs the raw COSMIC export rather than our ingested tables.
Whether the harvest is complete: abstracts only, so a resistance substitution named only in a
figure or supplement is missed. Whether L528W's class-level cross-resistance to covalent agents
extends to ibrutinib specifically. Whether a different configuration or the charged-peptide option
would change any call — one configuration per target, the one Table 1 names. And the central
limitation of the whole exercise: with three sound variants and one prospective one, no rate
computed here has a usable confidence interval.
