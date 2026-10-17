# Table S7 — verified and corrected

Rebuilt 2026-10-04. Replaces the earlier length-based audit: every row is now settled by
**sequence checksum**, so no assignment rests on two sequences happening to be the same length.

## Method

For each FASTA the pipeline actually reads, the CRC64 of the first record is computed
(`Bio.SeqUtils.CheckSum.crc64`) and compared with the CRC64 UniProt publishes for the
accession Table S7 cites (`rest.uniprot.org/uniprotkb/<acc>?format=json&fields=sequence`).
Where they differ, the sequence is identified through **UniParc** by checksum
(`rest.uniprot.org/uniparc/search?query=checksum:<CRC64>`), which returns the UniProtKB
isoform holding exactly that sequence. RefSeq protein↔mRNA pairings and isoform labels are
quoted from each UniProt entry's own cross-references, so no NCBI lookup is involved.

Script: `s7_fix.py`. Full table: `Table_S7_corrected.csv`.

## The result that unblocks the external comparison

**All six analysed targets sit on the canonical UniProt sequence, checksum-exact:**

| target | aa | CRC64 | UniProt canonical |
|---|---|---|---|
| ABL1_ENST00000318560 | 1130 | 85FE6C1C0E483EA2 | P00519 |
| ALK | 1620 | 0733D6C4FD212F41 | Q9UM73 |
| EGFR | 1210 | D8A2A50B4EFB6ED2 | P00533 |
| KIT | 976 | 81B0CD76817F3454 | P10721 |
| NT5C2 | 561 | 4C27D762575E0EA2 | P49902 |
| SMO | 787 | 8B4C459B34D13F83 | Q99835 |

And the variant numbering agrees with those sequences: re-checking every `p.X<pos>Y` call
against the FASTA the pipeline reads gives **0 mismatches** out of 2676 parseable records
(ABL1 1272, EGFR 653, KIT 493, ALK 150, NT5C2 65, SMO 43; 19 records unparseable, the
`p.E355X`-type calls).

So AlphaMissense / dbNSFP / ESM-2 scores keyed on UniProt canonical join to our variants with
**zero offset** for all six targets. The feared +19 ABL1 shift does not touch the analysis.

## Five rows to fix

The errors are all in rows that appear in S7 but carry no variant data. Three are the same
mistake: the **canonical** accession printed against the **non-canonical isoform** row.

1. **ABL1** (1149 aa) — S7 prints `P00519`; correct is **`P00519-2`** (isoform IB,
   UniParc UPI000013E4DE, RefSeq NP_009297.2 ↔ NM_007313.3).
   `P00519` is the 1130 aa isoform IA and belongs to the `ABL1_ENST00000318560` row.
2. **ABL1_ENST00000318560** (1130 aa) — S7 prints only `NM_005157.4`; should carry
   **`P00519`**. The RefSeq protein is right (NP_005148.2) but UniProt now pairs it with
   NM_005157**.6**, not .4.
3. **BTK** (693 aa) — S7 prints `Q06187`; correct is **`Q06187-2`**
   (UPI000021216F, NP_001274273.1 ↔ NM_001287344.2).
   **BTK_ENST00000308731** (659 aa) is the one that is `Q06187`; S7 writes `Q06187.3`, which
   is an entry-version string, not an isoform identifier, and should not appear.
4. **MET** (1408 aa) — S7 prints `P08581`; correct is **`P08581-2`**
   (UPI000014033E, NP_001120972.1 ↔ NM_001127500.3).
   **MET_ENST00000397752** (1390 aa) is `P08581`; its `NM_000245.4` citation is correct, it
   just needs the UniProt accession added. *(This row was not flagged in the first audit —
   the checksum check found it.)*
5. **EGFR_ENST00000454757** (1165 aa) — S7 cites `NM_005228.3`, but UniProt pairs
   NM_005228 with NP_005219.2, the **1210 aa** protein, i.e. the other EGFR row. The 1165 aa
   sequence is UniParc **UPI000387BFB2**, whose only UniProtKB mapping is **E9PFD7**, now an
   **inactive entry, deleted from its Ensembl source**. This row therefore has no current
   UniProt or RefSeq accession at all and cannot be joined to any UniProt-keyed resource.

## Two pairs of duplicate rows

Identical sequences, same CRC64, listed twice:

    5455C57AB0CCCAA7   ESR1  = ESR1_ENST00000206249    595 aa
    4C27D762575E0EA2   NT5C2 = NT5C2_ENST00000404739   561 aa

Any count or average over S7 rows double-counts these two proteins. (NT5C2 is one of the six
analysed targets, so this also means `NT5C2` and `NT5C2_ENST00000404739` are the same
experiment under two names.)

## Rows that are correct as printed

ALK (Q9UM73), ALK_ENST00000618119 (A0A087WZL3), EGFR (P00533), ESR1 / ESR1_ENST00000206249
(P03372), FLT3 (P36888), IDH1_ENST00000415913 (O75874), KIT (P10721), NT5C2 /
NT5C2_ENST00000404739 (P49902), SMO (Q99835) — all checksum-exact against the cited accession.

## Sequence relationships between the paired rows (offline, exact)

    BTK   / BTK_ENST     693 / 659   second = first minus 34 N-terminal residues
    ALK   / ALK_ENST    1620 / 1243  second = first minus 377 N-terminal residues
    ESR1  / ESR1_ENST     595 / 595  identical
    NT5C2 / NT5C2_ENST    561 / 561  identical
    ABL1  / ABL1_ENST   1149 / 1130  1 aa common prefix, 1104 aa common suffix
    EGFR  / EGFR_ENST   1210 / 1165  141 aa prefix, 1023 aa suffix
    MET   / MET_ENST    1408 / 1390  755 aa prefix, 636 aa suffix

The last three are alternative splice forms, so no constant shift maps a position in one onto
the other — cross-transcript pooling of those rows is wrong by construction.

## What was not verified here

- **Ensembl transcript ids were not checked.** `rest.ensembl.org` and `www.ensembl.org` are
  both unreachable from this environment, so whether `ENST00000454757` is still current (and
  what it now encodes) is open. The E9PFD7 deletion notice suggests it changed.
- **RefSeq versions** are reported as UniProt currently pairs them. S7's `NM_005157.4` and
  `NM_005228.3` are older versions; whether the protein changed between versions was not
  checked against NCBI (eutils needed an approval that was not given).
- **Only the first FASTA record** of each file is used, matching what the pipeline does.
  FLT3 (11 records), KIT (8) and SMO (3) have further records that were not examined.

## Also noticed while reading the Supplementary

Supplementary Table 2 lists a different best configuration for most targets than Table 1
(ABL1 15/6 vs 33/5, ALK 13/7 vs 33/4, FLT3 37/5, …) because it optimises the "Not resistant"
class instead. Two different "best" configurations for the same data is itself evidence that
the reported numbers are selection-driven; whichever is kept, the Methods must say which class
the configuration was selected on.
