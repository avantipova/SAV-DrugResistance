# External source data

The files in this directory are **not** redistributed with the repository, because they belong to
other projects under their own terms. Each one is a single download, and `check_inputs.py` at the
repository root tells you which are missing and which verifiers they unlock.

Everything the repository *can* ship is already here or under `inputs/` — see `REPRODUCING.md`.

## Already included (licence permits redistribution)

| path | source | licence |
|---|---|---|
| `source-data/civic/01-Oct-2026-*.tsv` | CIViC release 01-Oct-2026 | CC0 |
| `../inputs/fasta/*.faa` | UniProt / NCBI Protein reference sequences | CC BY 4.0 (UniProt) |

## You need to download these

### `primeedit/RS_<drug>.csv` — six files
ABL1 kinase-domain prime-editing saturation scan. One file per inhibitor:
`RS_Imatinib.csv`, `RS_Nilotinib.csv`, `RS_Dasatinib.csv`, `RS_Bosutinib.csv`,
`RS_Ponatinib.csv`, `RS_Asciminib.csv`.
Source: the supplementary data of the published prime-editing scan (Jung et al.).
Unlocks: `results/code/verify_primeedit.py`, `results/code/verify_activark.py`.

### `mdrdb/core.tsv` and `mdrdb/tki_verified.csv`
MdrDB measured binding free-energy changes for mutant protein–ligand complexes.
Source: the MdrDB project's own distribution.
Unlocks: `multipass/verify_mdrdb.py`, `multipass/verify_tki.py`.

### `drsp/s001.zip` and `drsp/x/pharmaceuticals-1557619-supplementary/Dataset.xlsx`
The DRSP published dataset. `s001.zip` is the supplementary archive as downloaded;
`Dataset.xlsx` is the workbook inside it, at the path the archive itself creates.
Unlocks: `multipass/verify_drsp.py`.

### `source-data/alphamissense/AF-<ACC>-F1-aa-substitutions.csv` — six files
Per-protein AlphaMissense predictions from AlphaFold DB (`amAnnotationsUrl`), for the six
accessions of the position-grouped comparison. **CC BY-NC-SA 4.0** — a non-commercial
share-alike licence, which is why these are not carried in an MIT-licensed repository.
Unlocks, together with the next entry: `level2/verify_level2.py`.

### `source-data/dbnsfp/dbnsfp_myvariant_6genes.tsv`
REVEL, CADD, PolyPhen-2 and SIFT scores for the six genes, as retrieved from MyVariant.info
(`dbnsfp.genename:<gene>`, scrolled in full). Subject to dbNSFP's own terms.
Note that the version is whatever MyVariant served at retrieval time; it is not pinned, so a
fresh retrieval may not reproduce the published figures exactly.
Unlocks, together with the previous entry: `level2/verify_level2.py`.
