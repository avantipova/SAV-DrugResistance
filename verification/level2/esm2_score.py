"""ESM-2 variant effect scores for the six targets, wt-marginal scoring.

score(variant) = log p(mut | full wild-type sequence, position i)
               - log p(wt  | full wild-type sequence, position i)

This is the "wt-marginals" estimator of Meier et al. (ESM-1v): one forward pass over the
unmasked wild-type sequence, then read the masked-LM head's distribution at each position.
It needs no label and no fine-tuning, so ESM-2 enters the comparison exactly like
AlphaMissense and REVEL do -- as an off-the-shelf predictor.

The model is run on the FULL canonical protein, not the peptide window the MNA base uses,
because that is how ESM-2 is meant to be applied; the window is a MultiPASS implementation
detail, not a property of the variant.
"""
import os as _os

def _savdr_root():
    """Locate the repository root by walking up to the .savdr-root marker."""
    d = _os.path.dirname(_os.path.abspath(__file__))
    while True:
        if _os.path.exists(_os.path.join(d, '.savdr-root')):
            return d
        parent = _os.path.dirname(d)
        if parent == d:
            raise SystemExit('cannot locate repository root: .savdr-root marker not found')
        d = parent

_SAVDR_ROOT = _os.environ.get('SAVDR_ROOT') or _savdr_root()

import os, sys, glob, numpy as np, pandas as pd, torch
from transformers import AutoTokenizer, EsmForMaskedLM

U = f'{_SAVDR_ROOT}/external/source-data'
MODEL = sys.argv[1] if len(sys.argv) > 1 else '/home/claude/esm2_650M'
OUT = sys.argv[2] if len(sys.argv) > 2 else f'{_SAVDR_ROOT}/level2/esm2_scores.csv'
TARGETS = ['ABL1_ENST00000318560', 'ALK', 'EGFR', 'KIT', 'NT5C2', 'SMO']

def first_record(path):
    s, started = [], False
    for line in open(path):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)

torch.set_num_threads(2)
tok = AutoTokenizer.from_pretrained(MODEL)
model = EsmForMaskedLM.from_pretrained(MODEL, dtype=torch.float32).eval()
print('model loaded:', sum(p.numel() for p in model.parameters()) / 1e6, 'M params', flush=True)

rows = []
for t in TARGETS:
    seq = first_record(f'{U}/fasta/{t}.faa')
    # ESM-2's positional limit is 1024 tokens; longer proteins are scored in overlapping
    # windows of 1022 residues with 256-residue overlap, each position taken from the window
    # where it sits furthest from an edge.
    W, OV = 1022, 256
    best = {}
    starts = list(range(0, max(1, len(seq) - W + 1), W - OV))
    if starts[-1] + W < len(seq): starts.append(len(seq) - W)
    if len(seq) <= W: starts = [0]
    for st in starts:
        chunk = seq[st:st + W]
        enc = tok(chunk, return_tensors='pt')
        with torch.no_grad():
            logits = model(**enc).logits[0]            # (1 + len + 1, vocab)
        lp = torch.log_softmax(logits, dim=-1)
        for j in range(len(chunk)):
            gpos = st + j                              # 0-based position in the protein
            edge = min(j, len(chunk) - 1 - j)
            if gpos in best and best[gpos][0] >= edge: continue
            best[gpos] = (edge, lp[j + 1])             # +1 for the <cls> token
    print(f'{t}: {len(seq)} aa, {len(starts)} window(s), {len(best)} positions scored', flush=True)
    for gpos, (_, dist) in best.items():
        wt = seq[gpos]
        wt_id = tok.convert_tokens_to_ids(wt)
        for mu in 'ACDEFGHIKLMNPQRSTVWY':
            if mu == wt: continue
            mu_id = tok.convert_tokens_to_ids(mu)
            rows.append((t, f'{wt}{gpos+1}{mu}', float(dist[mu_id] - dist[wt_id])))
df = pd.DataFrame(rows, columns=['target', 'protein_variant', 'esm2_llr'])
df.to_csv(OUT, index=False)
print('wrote', OUT, len(df), 'rows', flush=True)
