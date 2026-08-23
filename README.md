# SAV-DrugResistance

Data preparation and analysis code for **SAV-DrugResistance** — a web
application that predicts whether a somatic amino acid substitution (AAS) in an
anticancer drug target confers resistance to a given drug.

Somatic missense mutations in drug targets are turned into short peptide
fragments centred on the substituted residue, converted to structural formulas,
and used to train sequence–structure–property relationship (SSPR)
classification models with Multilevel Neighborhoods of Atoms (MNA) descriptors.

The service is freely available at
<https://way2drug.com/SAV-DrugResistance/>.

> This repository accompanies a manuscript under submission. Please cite the
> paper once published.

## Pipeline

| Stage | What happens |
| --- | --- |
| Mutation prep | Missense substitutions from COSMIC exports are matched against known resistance mutations, annotated with the causing drug, and cut into peptide windows centred on the mutated residue |
| Filtering | Targets with enough annotated mutations are kept; repeated records of the same mutation are collapsed to one set of peptide lengths |
| Structure generation | Peptides → RDKit molecules → SDF, optionally with formal charges assigned (`chargePeptide`) |
| Model training | **MultiPASS** trains one classification model per target–drug pair, sweeping peptide length and MNA descriptor level |
| Model evaluation | MultiPASS run logs are parsed for AUC (IAP) per target × drug × peptide length × MNA level |
| Figures | Mutation counts per target, AUC vs MNA scatter grids, performance heatmaps |

## Layout

```
notebooks/aas_analysis.ipynb   the analysis
scripts/multipass/*.ps1        drive MultiPASS over the SDF files
data/                          see data/README.md
results/                       figures
requirements.txt
```

The notebook does not run the whole pipeline on its own. It prepares peptides
and writes `data/sdf/`, then `scripts/multipass/run_mna_models.ps1` feeds those
to MultiPASS on Windows and collects `ConsolidatedOutput.txt`, and the notebook
picks the results back up from `data/multipass results/`.

## Data

Amino acid substitutions and their drug resistance annotations come from
**COSMIC** (Catalogue Of Somatic Mutations In Cancer), release **v100** for
model training, with **v103/v104** used for independent validation of variants
not seen during training. Full-length protein sequences are taken from UniProt
and NCBI Protein according to the corresponding transcripts.

COSMIC requires registration and its licence does not permit redistribution, so
the exports are not included here — obtain them from
<https://cancer.sanger.ac.uk/cosmic>. The datasets generated in this study are
deposited in Zenodo: <https://zenodo.org/records/21166107>.

Place the data at `data/`, or point `AAS_DATA_DIR` at wherever it lives. See
[data/README.md](data/README.md) for the expected layout.

## MultiPASS

Model training itself is **not** performed by this repository. MultiPASS is a
specialised extension of the PASS (Prediction of Activity Spectra for
Substances) methodology implementing a modified naïve Bayes approach; it runs
on Windows and is not publicly distributed.

For access to MultiPASS, please contact **Alexey A. Lagunin**
(<alexey.lagunin@ibmc.msk.ru>), Institute of Biomedical Chemistry / Pirogov
Russian National Research Medical University, Moscow.

The `scripts/multipass/*.ps1` scripts included here only automate calling
`OLMPASSdoMNABase.exe` over the SDF files and collecting its logs.

## Setup

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name sav-drugresistance
```

The SDF-generation section additionally needs pyPept:

```bash
pip install "pypept @ git+https://github.com/Boehringer-Ingelheim/pyPept.git"
```

Its command-line interface needs a source checkout — clone it to `external/pyPept`
(or set `PYPEPT_DIR`):

```bash
git clone https://github.com/Boehringer-Ingelheim/pyPept.git external/pyPept
```

## Running

```bash
jupyter lab notebooks/aas_analysis.ipynb
```

The first cell resolves paths and fails fast with a clear message if `data/` is
empty. Override locations with environment variables:

| Variable | Default | Meaning |
| --- | --- | --- |
| `AAS_DATA_DIR` | `<repo>/data` | the dataset |
| `AAS_RESULTS_DIR` | `<repo>/results` | where figures and derived files go |
| `PYPEPT_DIR` | `<repo>/external/pyPept` | pyPept source checkout |
