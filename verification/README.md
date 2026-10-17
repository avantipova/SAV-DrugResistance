# Verification

Everything in this directory exists to answer one question a reviewer is entitled to ask: *were
these checks actually run, and can I re-run them?*

Each number in the manuscript is re-derived here from stored data by a script that asserts the
exact value. A script exits non-zero if any assertion fails, so a clean run is a machine-checked
statement, not a claim. `notebook/verification.ipynb` carries the outputs as committed, so the
evidence can be read without installing or downloading anything.

```
verification/
├── README.md               this file
├── check_inputs.py         run this first: says what is present and what is missing
├── notebook/
│   └── verification.ipynb  runs every script and plots the four headline results
├── reports/                one REPORT_*.md per analysis, with its reasoning and caveats
├── results/data/           derived tables and run logs the scripts read
├── inputs/                 model scores, fold assignments and reference sequences
├── external/               empty: the datasets that are not ours to redistribute go here
└── <multipass|civic|level2|results|analysis|verification>/   the scripts themselves
```

## Quick start

```bash
pip install numpy pandas scipy scikit-learn            # that is the whole dependency list
python3 check_inputs.py                                # what is ready, what is missing
cd multipass && python3 verify_new.py                  # 179 assertions
```

`check_inputs.py` prints every verifier, whether it can run in this checkout, and for the ones
that cannot, exactly which file is missing and where to get it.

Each script must be run **from its own directory**, because a few read sibling files by relative
name. The scripts find the repository by walking up to the `.savdr-root` marker, so the checkout
can live anywhere; set `SAVDR_ROOT` to override.

## What runs with no downloads — 810 assertions

| from | script | what it checks |
| --- | --- | --- |
| `multipass/` | `verify_new.py` | position-grouped cross-validation against three baselines, calibration, leave-one-target-out, and that no held-out protein contributed to its own MNA base (179) |
| `multipass/` | `verify_external.py` | the SUSPECT-ABL and Dr. Kinase comparisons (49) |
| `civic/` | `verify_civic.py` | CIViC as an independent label source, including the 14 discordant calls (159) |
| `results/` | `verify_egfr_anomaly.py` | the EGFR activation-loop transfer result (112) |
| `results/code/` | `verify_optimism.py` | nested configuration-selection optimism on six targets, the three where it is not computable, and both correlations with exact permutation p-values (140) |
| `results/code/` | `verify_label_coverage.py` | documented resistance variants absent from the label set (83) |
| `results/code/` | `audit_cross.py` | the reports **against each other**, so a figure quoted twice cannot drift (88) |

## What needs a dataset first

None of these are ours to redistribute. Download the file, put it at the path given, and the
script runs unchanged.

| script | needs | where to get it |
| --- | --- | --- |
| `results/code/verify_primeedit.py`, `verify_activark.py`, `verify_ak3.py` | `external/primeedit/RS_<drug>.csv` for Imatinib, Nilotinib, Dasatinib, Bosutinib, Ponatinib, Asciminib | supplementary data of Jung *et al.*, *Nature Biomedical Engineering* (2025), [doi:10.1038/s41551-025-01531-4](https://doi.org/10.1038/s41551-025-01531-4) |
| `results/code/verify_activark.py`, `verify_ak3.py` | `external/activark/` — predictions for our ABL1, ALK, EGFR and KIT variants | submit the variant lists to <https://activark.russelllab.org> at the GBC setting (Singh *et al.*, *Genome Medicine* **17**, 133, 2025) |
| `multipass/verify_mdrdb.py`, `verify_tki.py` | `external/mdrdb/core.tsv`, `external/mdrdb/tki_verified.csv` | MdrDB release v1.0 (Yang *et al.*, *Commun. Chem.* **6**, 123, 2023), <https://quantum.tencent.com/mdrdb/> |
| `multipass/verify_drsp.py` | `external/drsp/s001.zip` and the `Dataset.xlsx` inside it | supplementary data of Lin *et al.*, *Pharmaceuticals* **15**, 136 (2022) |
| `level2/verify_level2.py` | `external/source-data/dbnsfp/dbnsfp_myvariant_6genes.tsv` and `external/source-data/alphamissense/AF-<ACC>-F1-aa-substitutions.csv` for the six accessions | dbNSFP via <https://myvariant.info> (`dbnsfp.genename:<gene>`); AlphaMissense per-protein files from AlphaFold DB (`amAnnotationsUrl`), CC BY-NC-SA 4.0 |

## Source data

**COSMIC exports are not redistributed in this directory.** COSMIC requires registration and its
licence does not permit redistribution — obtain the v100 exports from
<https://cancer.sanger.ac.uk/cosmic> as described in the repository root `README.md`.

What *is* here are the derived products of our own analysis: the per-variant model scores, the
cross-validation fold assignments, the curated variant tables the models were trained on, the
summary tables, and the run logs. Those are also archived, under CC BY 4.0, in the Zenodo record
for this work: <https://zenodo.org/records/21166107>.

Reference protein sequences in `inputs/fasta/` are from UniProt and NCBI Protein (UniProt data is
CC BY 4.0).

`inputs/suspect/` holds SUSPECT-ABL's predictions for our 280 ABL1 variants — the result of
submitting our own variant list to its public server — so that `verify_external.py` runs without
a re-submission. ActivArk's predictions are not included: its output flags which variants belong
to its own curated training set, which is its data rather than a prediction on ours, so the
scripts that use it ask you to query the server yourself.

## Why the fold assignment is a stored file

`inputs/cv_folds/fold_assignments.csv` holds the position-grouped fold assignment that produced
the published numbers, rather than letting each script rebuild it.
`sklearn.model_selection.StratifiedGroupKFold` does not guarantee the same split across
scikit-learn versions: rebuilding it made 12 of `verify_new.py`'s assertions disagree with the
published values on scikit-learn 1.7.2 while passing on 1.8.0. Reading the stored assignment,
the script returns 179/0 on both. If you regenerate the model scores, regenerate this file with
them — do not mix a new split with the stored scores.

The suite has been run on scikit-learn 1.7.2 / pandas 2.3.3 / numpy 2.2.6 / scipy 1.15.3 and on
scikit-learn 1.8.0 / pandas 3.0.2 / numpy 2.4.4 / scipy 1.17.1, with identical counts on both.

## What cannot be re-run here

The scripts that *generate* the model scores need MultiPASS, which is proprietary and not
distributed (available from A. A. Lagunin in a demo version on reasonable request) running under
Wine, plus ESM-2 weights for the sequence-model comparison. They are included for inspection and
still carry the working-directory paths of the machine they ran on:
`multipass/loto11.py`, `multipass/build11.py`, `multipass/loto11_summary.py`,
`multipass/drsp_compare.py`, `multipass/mdrdb_compare.py`, `results/code/pe_score.py`,
`results/code/prospective_score.py`, `results/code/within_target_cv.py`, `level2/esm2_score.py`.
`verification/verify_all.py` is also outside this path: it rebuilds each variant table from the
COSMIC source export and reads all 280 training SDFs, so it needs the COSMIC data.

Their outputs are what `inputs/` and `results/data/` contain, so every verifier above runs
without them.
