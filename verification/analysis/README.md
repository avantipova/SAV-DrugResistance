# SAV-DR re-analysis (no MultiPASS)

Data: `../data/{csv,fasta,mutations}` (from Drive folder AAS, 2026-09-30).

| script | output | what |
|---|---|---|
| 01_position_check.py | 01_position_check.csv | wt residue vs first FASTA record (as aas.ipynb uses) |
| 02_baselines.py | 02_baselines.csv | pos_only / subst_lr / full_gbm, random vs position-grouped CV, pooled + per-drug |
| 03_region_restricted.py | 03_region_restricted.csv | same, negatives restricted to resistance span +-10 aa |
| 04_sanity.py | stdout | no position leak in grouped folds; shuffled labels ~0.5; AUC cross-check |
| 05_sanity_all.py | 05_sanity_all.csv | the same checks on all 14 tasks + 6 region-restricted |

Run `python3 02_baselines.py <task idx...>` in chunks (each task ~40 s).

| 06_morgan.py | 06_morgan.csv | RDKit peptides (mutant residue substituted) -> Morgan count FP r=2/2048; BernoulliNB + L2 logistic; lengths 11/17/25 all reported; pos_only on the same variant subset |

06_morgan.py env vars: LENS=11,17,25 (which lengths), WITH_POS=0 (skip the pos_only reference).
Task idx 0-5 = pooled targets with full negatives, 6-11 = same targets, negatives restricted to
the resistance span +-10 aa. Morgan fingerprints are an OPEN ANALOGUE of atom-neighbourhood
descriptors, not MNA; no MultiPASS is involved anywhere in this folder.
