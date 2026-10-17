# SPLDExtraTrees and the TKI benchmark — a fourth head-to-head, and we lose it

Generated 2026-10-10. Found by a second literature sweep after `REPORT_MISSED_METHODS.md`, and
runnable without the authors' code because their **test set** is public.

**SPLDExtraTrees** (Zhu *et al.*, *Briefings in Bioinformatics* 2022, bbac050; preprint
arXiv:2111.08008) predicts the change in inhibitor binding affinity caused by a kinase mutation
from 146 structure-derived features — Rosetta REF15 energies, FoldX, PLIP contacts, AutoDock
Vina. It is **not** in the manuscript's Table S3 and should be, and it publishes no code. But its
benchmark, the **TKI dataset** — 144 measured ΔΔG values over ABL kinase point mutations and
eight tyrosine-kinase inhibitors — is redistributed inside MdrDB, which we already have. So our
scores can be placed beside their published correlations without running their model.

**Conclusion first: on this benchmark our score carries essentially no signal.** Drug-matched,
our per-drug Pa correlates with the measured ΔΔG at Spearman **0.163 (95% CI −0.072 to 0.372,
p = 0.14)**, against their published Pearson of **0.54** in the comparison that is fair to us and
0.74 in their best case, with Rosetta REF15 at 0.67. Our binary resistance Pa over all 143 rows
correlates **−0.129**. This does not say our method is broken; it says our score and the
biophysics of binding are close to unrelated at the substitution level, which is consistent with
everything else in this project.

## Published performance on the TKI dataset

Quoted from the paper, not recomputed. Mean over 20 repetitions.

| model | RMSE (kcal/mol) | Pearson | AUPRC |
|---|---|---|---|
| SPLDExtraTrees, scenario 1 — trained on Platinum only, **no kinase data** | 0.73 | **0.54** | 0.50 |
| SPLDExtraTrees, scenario 2 — Platinum + 2 TKI samples per inhibitor | 0.70 | 0.60 | 0.55 |
| SPLDExtraTrees, scenario 3 — trained and tested on the TKI dataset | 0.58 | **0.74** | 0.56 |
| Rosetta REF15, for reference | 0.72 | 0.67 | 0.53 |

Scenario 1 is the comparison to use: neither model has seen kinase ΔΔG. Scenario 3 is
cross-validated on the very set it is scored on, so it is their upper bound, not a fair target.

## A data-quality finding about MdrDB that had to be fixed first

MdrDB spreads the TKI set across **four different UniProt accessions**, apparently by which
template structure was used:

| accession | inhibitors | rows |
|---|---|---|
| P00519 — human ABL1 | axitinib, bosutinib, dasatinib | 68 |
| **P00520 — mouse Abl1** | imatinib, ponatinib | 42 |
| **Q835X8** | nilotinib | 21 |
| **UNDEFINED** | erlotinib, gefitinib | 13 |

All 144 rows are in fact human ABL1 substitutions: each was checked against the human sequence
and **all 144 wild-type residues match, at offset 0**, so the set is numbered on the canonical
1130-residue sequence and the accession labels are simply wrong for 76 of the rows. Every one of
the 31 distinct substitutions lies in the kinase domain, positions 244–486, and 143 of 144 are in
our own COSMIC-derived variant table (only E292L is absent).

This matters beyond this report: filtering MdrDB by human accession — the conservative choice,
and what `mdrdb_compare.py` does — silently keeps only 68 of these 144 rows and 3 of the 8
inhibitors. `REPORT_MDRDB.md`'s TKI subset is that conservative slice. The numbers below use all
144, justified by the wild-type check rather than by the accession label.

The inhibitors themselves were resolved by **PubChem CID**, not by ligand code, because the
`DRUG` column holds PDB chemical-component ids: STI→imatinib (CID 5291), NIL→nilotinib,
1N1→dasatinib, DB8→bosutinib, AXI→axitinib, 0LI→ponatinib, AQ4→erlotinib, IRE→gefitinib. That
reproduces exactly the eight inhibitors the paper names.

## Result

Nothing is refitted. Per-drug Pa is the ABL1 multiclass run, rank-averaged over the five
fold-assignment seeds; the binary Pa is the main `full`-setting run.

| panel | n | Spearman | 95% CI | p | Pearson |
|---|---|---|---|---|---|
| **drug-matched Pa, the 4 inhibitors we label** | 84 | **0.163** | −0.072 .. 0.372 | 0.14 | 0.122 |
| — dasatinib | 21 | **0.507** | 0.066 .. 0.816 | 0.019 | 0.480 |
| — imatinib | 21 | 0.171 | −0.303 .. 0.605 | 0.46 | 0.156 |
| — bosutinib | 21 | −0.038 | −0.567 .. 0.515 | 0.87 | 0.061 |
| — nilotinib | 21 | **−0.406** | −0.751 .. 0.049 | 0.068 | −0.433 |
| binary resistance Pa, all 8 inhibitors | 143 | −0.129 | −0.286 .. 0.032 | 0.13 | −0.027 |
| binary Pa, the same 4 inhibitors | 84 | −0.117 | −0.320 .. 0.091 | 0.29 | −0.014 |
| binary Pa, the 4 inhibitors we do not label | 59 | −0.031 | −0.290 .. 0.211 | 0.81 | 0.091 |

Pooled, the interval includes zero. Per drug the four estimates range from +0.507 to −0.406: one
nominally positive at p = 0.019 (dasatinib), one nominally *negative* at p = 0.068 (nilotinib).
Across four tests on 21 points each that pattern is what noise looks like, and the honest reading
is that there is no consistent relationship, not that dasatinib works and nilotinib is inverted.

## Why this is expected rather than alarming, and what to write

Three things make this result the one the rest of the project predicts:

1. **The TKI set sits entirely inside the kinase domain** (positions 244–486). That is exactly
   the `region` regime where our own margin sweep already put ABL1 at 0.589–0.599 and where
   SUSPECT-ABL also lands near chance. The position shortcut is gone, and with it most of our
   signal.
2. **ΔΔG of binding is not the quantity our labels encode.** MdrDB's own threshold calls only
   8.2% of its samples resistant at a ten-fold affinity drop, and many clinical resistance
   mutations act by stabilising the active conformation or through allostery rather than by
   weakening inhibitor binding. A method trained on curated resistance annotations is not
   obliged to track ΔΔG, and ours does not.
3. **Structure-based methods should win on a structural quantity.** Rosetta REF15 at 0.67 is a
   physics calculation on the complex. That SPLDExtraTrees matches it with learned features is
   their contribution; that a sequence-only classifier does not is not a surprise.

For the manuscript:

- **Add SPLDExtraTrees to Table S3** and state plainly that on its TKI benchmark MultiPASS does
  not predict ΔΔG (Spearman 0.163, CI spanning zero) while it reports 0.54–0.74. Claiming
  anything else here is indefensible, and the comparison is easy for a reviewer to repeat because
  the data is public.
- **Use it to sharpen the scope claim, which is the real point.** SAV-DrugResistance predicts a
  curated resistance *annotation* from sequence alone, across 18 targets and 18 drugs, through a
  public server. SPLDExtraTrees predicts a *binding-affinity change* for one kinase and eight
  inhibitors, needs a co-crystal or docked complex plus Rosetta/FoldX/Vina features, and ships no
  code. These are different tasks with different inputs, and the paper is stronger saying so than
  pretending they are comparable on accuracy.
- **Correct `REPORT_MDRDB.md`'s TKI row in passing:** its ρ = 0.360 on 30 mutations was computed
  on the 68 human-accession rows after collapsing to one score per mutation. The drug-matched
  figure on all four labelled inhibitors is 0.163. Both are reported; neither is wrong, they
  answer slightly different questions, and the difference is itself a lesson about how much the
  accession filter and the per-mutation collapse move the answer.

## Scope of verification

`verify_tki.py`, **36 checks, all passing**.

**Verified.** That MdrDB's TKI slice holds 144 single-substitution rows, 31 distinct
substitutions and eight inhibitors, and that all 144 carried through. The four-accession split,
with imatinib and ponatinib on the mouse accession, nilotinib on Q835X8 and erlotinib/gefitinib
on none, and that a human-accession filter would keep only 68 rows. Every ligand code resolved to
a drug name by PubChem CID, with imatinib's full canonical SMILES checked against CID 5291 and
that CID appearing only as ligand STI. That all 144 wild-type residues match the human ABL1
sequence at offset 0, that 143 of 144 are in our variant table and the one absent is E292L, and
that the 31 substitutions span positions 244–486. That per-drug scores exist for exactly 84 rows
over exactly the four inhibitors we label, 21 rows each, and that every per-drug value equals the
five-seed rank average the main pipeline produces. The pooled Spearman recomputed as a Pearson
correlation on ranks rather than through `scipy.spearmanr`; all four per-drug coefficients; that
dasatinib is the only one nominally positive at p < 0.05; that the pooled interval includes zero;
and that our largest absolute correlation is below their smallest published Pearson.

**Not verified.** SPLDExtraTrees' own numbers — quoted from the paper, with no code to run and no
attempt to reproduce them. Whether their Pearson on the 144 rows would change if restricted to
the 84 we can match, which cannot be computed without their per-row predictions. Their feature
pipeline (Rosetta, FoldX, PLIP, Vina), none of it reproduced. How MdrDB derived each TKI ΔΔG, and
whether its redistribution of the TKI set is faithful to the original — only the mutations and
drug identities were checked, not the ΔΔG values themselves. Whether the four inhibitors we
label are the four for which our multiclass model is best trained; imatinib has 82 positives and
bosutinib 5, so the per-drug estimates rest on very different amounts of training signal, and
that was not factored in.

## Files

| file | contents |
|---|---|
| `tki_overlap.csv` | all 144 TKI rows: MdrDB's accession and mutation, ours, the resolved drug, ΔΔG, both of our scores |
| `tki_compare.csv` | the correlation panels with CIs |
| `tki_compare.py` | the comparison |
| `verify_tki.py` | the 36 checks |
