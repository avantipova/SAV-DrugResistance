# Running MultiPASS without Windows (Wine)

Verified 2026-10-02 in the Cowork cloud container.

Setup
    apt-get install -y wine
    dpkg --add-architecture i386 && apt-get update && apt-get install -y wine32:i386
    export WINEARCH=win32 WINEPREFIX=<dir>/.wine WINEDEBUG=-all
    wineboot -i
Binaries: Downloads/2024_MultiPASS/*.exe (32-bit console, i386).

Works   OLMPASSdoMNABase.exe <task.txt>   full cycle: create/add/train/save/20-fold validation/saveAsSAS
        OLMPASSofMOL.exe <in.sdf> <out.txt> /B<base>.MSAS   -> Pa/Pi, FIRST record only, ~5 s per call
Known   doMNABase raises EAccessViolation on exit AFTER every file is written (same on Windows per the
        project's own logs) -- cosmetic.
Broken  OLMPASSofSDF.exe runs, exits 0, writes no output file. Unresolved.

build_sdf.py  rebuilds the pipeline's SDF (port of extract_substring + the cell-3 labelling).
              Multi-drug variants: the Drug field holds one drug per line; MultiPASS then counts the
              variant under each single drug AND creates a combined class, so per-drug N matches Table 1.

## Reproduction of Table 1, ABL1_ENST00000318560, peptide length 17, MNA level 9

drug                             N (mine/Table 1)   20-fold IAP mine   Table 1 20-fold
Imatinib                         83 / 87            0.8865             0.8887
Tyrosine kinase inhibitor - NS   48 / 48            0.8420             0.8435
Nilotinib                         9 / 9             0.8694             0.8710
Dasatinib                        10 / 10            0.8796             0.9144  (Table 1 config 11/6)
Bosutinib                         5 / 5             0.7813             0.7819  (Table 1 config 17/7)

The N gap for Imatinib is the four records p.E279Z, p.E355X, p.E450X, p.F359X: X/Z are not amino
acids and RDKit returns None for them, so they cannot enter an SDF here -- yet the published N is 87.

## MultiPASS on the same position-grouped folds as the baselines

folds.py trains an MSAR/MSAS base on the training folds only and scores held-out variants with
OLMPASSofMOL /C, taking Pa('Resistant'). Binary task: resistant to any drug vs not.
compare.py then runs the baselines on exactly those variants, same split (StratifiedGroupKFold,
5 folds, seed 0, one repeat). ABL1_ENST00000318560, peptide length 17, MNA level 9.

Negatives from the whole protein (101 pos / 611 neg)
  pos_only   AUC 0.830   AUPRC 0.318
  subst_lr   AUC 0.681   AUPRC 0.272
  full_gbm   AUC 0.869   AUPRC 0.446
  MultiPASS  AUC 0.851   AUPRC 0.382

Negatives from the resistance span +-10 aa only (101 pos / 210 neg)
  pos_only   AUC 0.454   AUPRC 0.299
  subst_lr   AUC 0.574   AUPRC 0.405
  full_gbm   AUC 0.554   AUPRC 0.375
  MultiPASS  AUC 0.555   AUPRC 0.354

Table 1 reports 0.894 (LOO) / 0.889 (20-fold) for ABL1_ENST/Imatinib at this configuration, but that
is a random-split CV on the full-protein negatives with the configuration chosen by maximising the
same metric. Under a position-grouped split MultiPASS gets 0.851, and inside the resistance region
it gets 0.555.

## 5 seeds, bootstrap CIs, paired DeLong vs MultiPASS

folds.py <target> <length> <MNA> [region] [seed] -- seeds 0..4, StratifiedGroupKFold(5, seed).
summarize.py rank-averages each model's out-of-fold scores across seeds, bootstraps the AUC CI over
variants (600 resamples, stratified) and runs a paired DeLong test of MultiPASS against each baseline.
The DeLong code is self-tested: identical scores give p = 1.000 and AUCs matching sklearn.
Full table: summary.csv.

ABL1_ENST00000318560, L17/MNA9, negatives from the whole protein (101 pos / 611 neg)
  model      AUC mean+-sd over seeds   pooled [95% CI]        AUPRC   DeLong p vs MultiPASS
  pos_only   0.842 +- 0.011            0.842 [0.810, 0.870]   0.354   0.078
  subst_lr   0.657 +- 0.023            0.665 [0.600, 0.730]   0.281   <0.001
  full_gbm   0.866 +- 0.007            0.875 [0.847, 0.901]   0.504   0.209
  MultiPASS  0.857 +- 0.005            0.860 [0.833, 0.888]   0.400   --

ABL1_ENST00000318560, L17/MNA9, negatives from the resistance span +-10 aa (101 pos / 210 neg)
  pos_only   0.523 +- 0.038            0.529 [0.463, 0.597]   0.336   0.028
  subst_lr   0.570 +- 0.008            0.576 [0.503, 0.651]   0.414   0.796
  full_gbm   0.597 +- 0.025            0.610 [0.547, 0.675]   0.429   0.491
  MultiPASS  0.582 +- 0.024            0.589 [0.518, 0.658]   0.396   --

ALK, L33/MNA4, region (24 pos / 137 neg)
  pos_only   0.686 +- 0.024            0.697 [0.589, 0.799]   0.246   0.601
  subst_lr   0.510 +- 0.020            0.508 [0.409, 0.623]   0.148   0.059
  full_gbm   0.677 +- 0.029            0.705 [0.611, 0.816]   0.292   0.612
  MultiPASS  0.665 +- 0.016            0.676 [0.562, 0.789]   0.225   --

EGFR, L37/MNA30, region (17 pos / 251 neg)
  pos_only   0.717 +- 0.055            0.722 [0.542, 0.877]   0.304   0.455
  subst_lr   0.669 +- 0.043            0.677 [0.524, 0.827]   0.141   0.854
  full_gbm   0.740 +- 0.022            0.756 [0.603, 0.884]   0.202   0.373
  MultiPASS  0.701 +- 0.025            0.698 [0.515, 0.859]   0.520   --

Reading: with whole-protein negatives MultiPASS (0.857) is indistinguishable from the gradient
boosting baseline (p = 0.21) and not significantly above the position-only model (p = 0.08); it beats
only the substitution-feature model that has no access to position (p < 0.001). Inside the resistance
region every method falls to 0.58-0.70 with CIs that include or nearly include 0.5, and no difference
between MultiPASS and either tree baseline is significant. The one significant region-restricted
result is MultiPASS > position-only on ABL1_ENST (p = 0.028), where position-only is at chance (0.523).
