#!/usr/bin/env python3
"""Report what this checkout can run, and what each remaining script is waiting for.

Looks at the filesystem only: computes nothing, changes nothing. Run it first.
"""
import os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))

OFFLINE = [
    ('multipass/verify_new.py',             'multipass',    179, 'position-grouped CV, baselines, calibration, LOTO'),
    ('multipass/verify_external.py',        'multipass',     49, 'SUSPECT-ABL and Dr. Kinase'),
    ('civic/verify_civic.py',               'civic',        159, 'CIViC as an independent label source'),
    ('results/verify_egfr_anomaly.py',      'results',      112, 'the EGFR activation-loop transfer result'),
    ('results/code/verify_optimism.py',     'results/code', 140, 'configuration-selection optimism'),
    ('results/code/verify_label_coverage.py','results/code',  83, 'literature variants absent from the labels'),
    ('results/code/audit_cross.py',         'results/code',  90, 'the reports against each other'),
]

NEEDS = [
    ('results/code/verify_primeedit.py', 'results/code',
     ['external/primeedit/RS_%s.csv' % d for d in
      ('Imatinib','Nilotinib','Dasatinib','Bosutinib','Ponatinib','Asciminib')],
     'ABL1 prime-editing saturation scan, Jung et al. Nat. Biomed. Eng. 2025, doi:10.1038/s41551-025-01531-4'),
    ('results/code/verify_activark.py', 'results/code',
     ['external/primeedit', 'external/activark'],
     'the scan, plus ActivArk predictions from activark.russelllab.org (GBC setting)'),
    ('results/code/verify_ak3.py', 'results/code',
     ['external/activark'],
     'ActivArk predictions for our ALK, EGFR and KIT variants'),
    ('multipass/verify_mdrdb.py', 'multipass', ['external/mdrdb/core.tsv'],
     'MdrDB v1.0, Yang et al. Commun. Chem. 6, 123 (2023), quantum.tencent.com/mdrdb'),
    ('multipass/verify_tki.py', 'multipass',
     ['external/mdrdb/core.tsv', 'external/mdrdb/tki_verified.csv'], 'MdrDB v1.0'),
    ('multipass/verify_drsp.py', 'multipass',
     ['external/drsp/s001.zip',
      'external/drsp/x/pharmaceuticals-1557619-supplementary/Dataset.xlsx'],
     'DRSP supplementary data, Lin et al. Pharmaceuticals 15, 136 (2022)'),
    ('level2/verify_level2.py', 'level2',
     ['external/source-data/dbnsfp/dbnsfp_myvariant_6genes.tsv',
      'external/source-data/alphamissense', 'level2/esm2_650M_scores.csv'],
     'dbNSFP via myvariant.info, AlphaMissense from AlphaFold DB, and the ESM-2 scores'),
]

SHIPPED = [
    ('inputs/variants',           'curated per-target variant tables'),
    ('inputs/fasta',              'UniProt / NCBI reference sequences'),
    ('inputs/loto_scores',        'leave-one-target-out per-variant scores'),
    ('inputs/cv_folds',           'cross-validation fold scores and the stored fold assignment'),
    ('inputs/loto_bases',         'LOTO scores and MNA base histories'),
    ('inputs/nested_checkpoints', 'per-fold checkpoints of the picked configurations'),
    ('results/data',              'derived result tables and run logs'),
    ('reports',                   'one report per analysis'),
    ('notebook/verification.ipynb','the notebook, with its outputs as committed'),
]

def present(rel):
    p = os.path.join(ROOT, rel)
    if os.path.isdir(p):
        return bool(os.listdir(p))
    return os.path.exists(p)

bad = 0
print('Shipped with the repository')
for rel, what in SHIPPED:
    ok = present(rel); bad += not ok
    print(f'  [{"ok " if ok else "MISSING"}] {rel:30s} {what}')

print('\nRuns with no downloads')
total = 0
for rel, cwd, n, what in OFFLINE:
    ok = os.path.exists(os.path.join(ROOT, rel)); bad += not ok; total += n
    print(f'  [{"ok " if ok else "MISSING"}] cd {cwd:14s} python3 {os.path.basename(rel):26s} {n:>4} assertions  {what}')
print(f'  -> {total} assertions when all seven are run')

print('\nWaiting on a dataset (see README.md for the source of each)')
for rel, cwd, needs, src in NEEDS:
    miss = [r for r in needs if not present(r)]
    if not miss:
        print(f'  [ok ] cd {cwd:14s} python3 {os.path.basename(rel):26s} inputs present')
        continue
    print(f'  [--- ] cd {cwd:14s} python3 {os.path.basename(rel)}')
    print(f'         {src}')
    for m in miss:
        print(f'           missing: {m}')

print('\nCOSMIC exports are not redistributed here; see README.md and the repository root README.')
if bad:
    print(f'\n{bad} shipped item(s) missing - this checkout looks incomplete.')
    sys.exit(1)
print('\nCheckout looks complete.')
