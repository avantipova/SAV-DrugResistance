# PremPLI and mCSM-lig — input prepared and verified, neither server completes a job

Generated 2026-10-09. Fifth and sixth resistance-relevant comparators, after SUSPECT-ABL,
Dr. Kinase, DRSP and D3EGFR.

Both predict the **change in protein–ligand binding affinity caused by a substitution**, which is
the quantity SUSPECT-ABL reports as ΔΔG, so a result from either would have slotted straight into
the existing ABL1/imatinib table and given three structure-based methods on one task.

- **PremPLI** (Sun *et al.*, *Commun Biol* **4**, 1311, 2021) — in the manuscript's Table S3.
- **mCSM-lig** (Pires, Blundell & Ascher, *Sci Rep* **6**, 29575, 2016) — **not** in Table S3,
  and it should be: same laboratory as SUSPECT-ABL, trained on the Platinum database, and its
  own abstract is about "the emergence of drug resistance".

**Conclusion first: the input is prepared and triple-checked, but neither server returned a
single prediction.** PremPLI's structure-preparation step answers **504 Gateway Time-out** after
60 seconds — for our structure *and for the server's own example*. mCSM-lig accepts a job and
issues a results URL, but the job was **still "Processing" 7.6 hours later**, and three further
submissions never even returned their redirect. Nothing is claimed about either method's accuracy
here; this is a report on reachability, and the prepared input is left in place so the run can be
started the moment either server answers.

## The structure, and a third independent numbering check

**PDB 1IEP** — ABL1 kinase domain in complex with imatinib (ligand code **STI**), the standard
co-crystal for this pair. Read directly from `files.rcsb.org`:

| | |
|---|---|
| chains | A and B, identical |
| chain A resolved residues | **225–498, 274 residues, no gaps** |
| hetero groups | STI (imatinib), CL, HOH |
| numbering | matches ours exactly |

The numbering check is the point. Our ABL1 variant set is numbered on ENST00000318560 (1130
residues, = UniProt P00519); DRSP turned out to use the 1149-residue isoform, +19. In 1IEP,
residue 244 is M, 253 is Y, 255 is E, **315 is T**, 351 is M, 359 is F and 396 is H — all seven
canonical positions carry the residue our table expects, so 1IEP is on our numbering and needs no
offset. Then, exhaustively: **all 280** of our ABL1 variants in the 242–493 window were checked
against chain A, giving **0 wild-type mismatches and 0 positions missing from the structure**.
The window our SUSPECT-ABL comparison already uses is therefore fully covered by this structure,
which means a PremPLI or mCSM-lig run would land on exactly the same 270-variant task — the
cleanest possible three-way comparison.

`abl1_1iep_input_set.csv` holds that prepared set, and
`abl1_1iep_chainA_residues.txt` the structure's residue list used for the check.

## PremPLI

Django app, `POST /research/PremPLI/upload_pdb` with a CSRF token, `pdb_id`, `pdb_mol`
(bioassembly or asymmetric unit) and `bioassembly`. Three submissions:

| submission | result | after |
|---|---|---|
| 1IEP, bioassembly 1 | **504 Gateway Time-out** (nginx) | 72 s |
| **1OSS, the server's own example**, bioassembly 1 | **504 Gateway Time-out** (nginx) | 60 s |
| 1IEP, asymmetric unit | connection dropped | 103 s |

That the server's own documented example also times out is what makes this a server problem
rather than a problem with our structure. The front page carries a notice about monthly reboots
on the first Monday at 12:00 Beijing time; 2026-10-09 is a Friday, so that is not the explanation.

## mCSM-lig

Single mutation per job: `POST /mcsm_lig/prediction` with `pdb_code`, `mutation`, `chain`,
`lig_id`, `affin_wt` (wild-type affinity in nM) and `run=single`. Submitted through the form so
that the browser follows the redirect the way a user would.

- **T315I** (1IEP, chain A, STI, wild-type affinity 37 nM) was accepted and redirected to
  `https://biosig.lab.uq.edu.au/mcsm_lig/results_prediction/1791535819.8`. That page said
  "Processing Submission … will automatically refresh after 30 seconds". It still said so
  21 minutes later, and still said so when re-checked at **7.6 hours** (submitted 08:50 UTC,
  re-checked 16:27 UTC). No error, no result at any point. The URL is recorded here because the
  page itself invites bookmarking, so it can be checked again later.
- Three further jobs (M351T, E255K, Y253H), submitted through hidden iframes, had not returned
  their redirect after 435 seconds, and were still hanging at **7.5 hours** — the POST itself
  never returns.
- `GET /mcsm_lig/run_example` returns instantly with a complete result (2Z4O, D30N chain A,
  ligand 065: predicted affinity change **−2.056 log fold, destabilising**; DUET stability change
  −0.087 kcal/mol; distance to ligand 2.814 Å). That told us the output format but it is a cached
  page, not a fresh computation, so it is no evidence that the queue works.

The output format matters for later: the resistance-relevant number is the predicted affinity
change in log fold units, negative meaning weaker binding, which is directionally the same
quantity as SUSPECT-ABL's ΔΔG and would be compared the same way.

## What this means for the manuscript

This is now a pattern rather than an anecdote. Of the six external resistance-prediction servers
examined in this project, **four will not produce a prediction today**:

| server | state on 2026-10-08/09 |
|---|---|
| SUSPECT-ABL | works — 280 variants × 4 inhibitors returned |
| Dr. Kinase | works, but answers a different question (hotspot loci) |
| KRDS | **ERR_EMPTY_RESPONSE** |
| D3EGFR | page up, prediction engine never completes (2 of 59 substitutions) |
| PremPLI | **504 on structure preparation, including its own example** |
| mCSM-lig | job accepted, still processing 7.6 h later |

For a web-server paper this is the strongest honest argument available, and it is checkable by a
reviewer: the field's tools stop working, and a maintained public server is itself the
contribution. It should be stated as reachability on a dated attempt, not as a verdict on any
method's accuracy — none of these observations tests whether PremPLI or mCSM-lig predicts well.

## Scope of verification

**Verified.** 1IEP's chain composition, residue range, absence of gaps and ligand codes, read
from the RCSB file itself. That seven canonical ABL1 positions carry the expected wild-type
residue in the structure, and that all 280 of our 242–493 variants match the structure's
wild-type residue with no missing positions. PremPLI's form fields and endpoint, read from its
own DOM, and the 504 responses with their timings, including on the server's own example
structure. mCSM-lig's form fields and endpoint, read from its own DOM; that a submission is
accepted and issues a results URL; that the results page was still in the processing state 7.6
hours after submission, with the elapsed time taken from the job id's own epoch timestamp; that
three further submissions had not returned after 7.5 hours; the example's output fields and
values.

**Not verified.** Whether either server completes jobs on a timescale longer than the 7.6 hours
observed — one day, one client. Whether the mCSM-lig job 1791535819.8 eventually finished; the URL is
recorded so it can be checked. Whether mCSM-lig's `affin_wt` input affects the predicted affinity
*change* or only the reported absolute affinity — two jobs with 37 nM and 370 nM were submitted
to test exactly this and neither returned, so the wild-type affinity to use for imatinib remains
an open choice that would have to be justified from the literature. Whether PremPLI would accept
an uploaded structure where the PDB-code path times out. Neither method's accuracy, on our data
or any other.

## Files

| file | contents |
|---|---|
| `abl1_1iep_input_set.csv` | the 280 ABL1 variants prepared for 1IEP chain A, with the structure's wild-type residue and the match flag |
| `abl1_1iep_chainA_residues.txt` | chain A's resolved residues, used for the numbering check |
