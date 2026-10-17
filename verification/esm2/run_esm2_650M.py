#!/usr/bin/env python3
"""ESM-2 650M wt-marginal variant scores for the six SAV-DrugResistance targets.

Run on a machine that has the weights locally. It needs only torch + transformers +
safetensors; the model directory is assembled here from the weights file you downloaded,
so nothing is fetched from the network.

    pip install torch transformers safetensors
    python3 run_esm2_650M.py

Writes esm2_650M_scores.csv next to itself: one row per (target, substitution) with
    esm2_llr = log p(mut | wild-type sequence, position i) - log p(wt | ..., i)
Lower = more deleterious. This is the "wt-marginals" estimator (Meier et al., ESM-1v):
one forward pass over the unmasked wild-type sequence, reading the masked-LM head at each
position, over the FULL canonical protein. Proteins longer than ESM-2's 1022-token context
are scored in overlapping windows, each position taken from the window where it sits
furthest from an edge.
"""
import os, sys, json, csv

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
FASTA = os.path.join(PROJECT, 'data', 'fasta')
WEIGHTS = os.path.expanduser('~/Downloads/esm2_650M.safetensors')
MODELDIR = os.path.join(HERE, 'esm2_650M')
OUT = os.path.join(HERE, 'esm2_650M_scores.csv')
TARGETS = ['ABL1_ENST00000318560', 'ALK', 'EGFR', 'KIT', 'NT5C2', 'SMO']

CONFIG = {"architectures": ["EsmForMaskedLM"], "attention_probs_dropout_prob": 0.0,
          "classifier_dropout": None, "emb_layer_norm_before": False, "esmfold_config": None,
          "hidden_act": "gelu", "hidden_dropout_prob": 0.0, "hidden_size": 1280,
          "initializer_range": 0.02, "intermediate_size": 5120, "is_folding_model": False,
          "layer_norm_eps": 1e-05, "mask_token_id": 32, "max_position_embeddings": 1026,
          "model_type": "esm", "num_attention_heads": 20, "num_hidden_layers": 33,
          "pad_token_id": 1, "position_embedding_type": "rotary", "token_dropout": True,
          "torch_dtype": "float32", "use_cache": True, "vocab_list": None, "vocab_size": 33}
VOCAB = "<cls>\n<pad>\n<eos>\n<unk>\nL\nA\nG\nV\nS\nE\nR\nT\nI\nD\nP\nK\nQ\nN\nF\nY\nM\nH\nW\nC\nX\nB\nU\nZ\nO\n.\n-\n<null_1>\n<mask>"
TOKCFG = {"model_max_length": 1000000000000000019884624838656, "tokenizer_class": "EsmTokenizer"}
SPECIAL = {"cls_token": "<cls>", "eos_token": "<eos>", "mask_token": "<mask>",
           "pad_token": "<pad>", "unk_token": "<unk>"}

if not os.path.exists(WEIGHTS):
    sys.exit(f'weights not found at {WEIGHTS}')
os.makedirs(MODELDIR, exist_ok=True)
json.dump(CONFIG, open(os.path.join(MODELDIR, 'config.json'), 'w'), indent=2)
open(os.path.join(MODELDIR, 'vocab.txt'), 'w').write(VOCAB)
json.dump(TOKCFG, open(os.path.join(MODELDIR, 'tokenizer_config.json'), 'w'))
json.dump(SPECIAL, open(os.path.join(MODELDIR, 'special_tokens_map.json'), 'w'))
link = os.path.join(MODELDIR, 'model.safetensors')
if not os.path.exists(link):
    try: os.symlink(WEIGHTS, link)
    except OSError:
        import shutil; shutil.copy(WEIGHTS, link)

import torch
from transformers import AutoTokenizer, EsmForMaskedLM

def first_record(path):
    s, started = [], False
    for line in open(path):
        if line.startswith('>'):
            if started: break
            started = True; continue
        if started: s.append(line.strip())
    return ''.join(s)

torch.set_num_threads(os.cpu_count() or 4)
tok = AutoTokenizer.from_pretrained(MODELDIR)
import transformers as _tf
_major = int(_tf.__version__.split('.')[0])
_kw = {'dtype': torch.float32} if _major >= 5 else {'torch_dtype': torch.float32}
model = EsmForMaskedLM.from_pretrained(MODELDIR, **_kw).eval()
print('transformers', _tf.__version__, '| torch', torch.__version__, flush=True)
print('model loaded:', round(sum(p.numel() for p in model.parameters()) / 1e6, 1), 'M params', flush=True)

rows = []
for t in TARGETS:
    seq = first_record(os.path.join(FASTA, t + '.faa'))
    W, OV = 1022, 256
    starts = list(range(0, max(1, len(seq) - W + 1), W - OV))
    if len(seq) <= W: starts = [0]
    elif starts[-1] + W < len(seq): starts.append(len(seq) - W)
    best = {}
    for st in starts:
        chunk = seq[st:st + W]
        enc = tok(chunk, return_tensors='pt')
        with torch.no_grad():
            lp = torch.log_softmax(model(**enc).logits[0], dim=-1)
        for j in range(len(chunk)):
            g = st + j
            edge = min(j, len(chunk) - 1 - j)
            if g in best and best[g][0] >= edge: continue
            best[g] = (edge, lp[j + 1])
    print(f'{t}: {len(seq)} aa, {len(starts)} window(s), {len(best)} positions', flush=True)
    for g, (_, dist) in best.items():
        wt = seq[g]; wid = tok.convert_tokens_to_ids(wt)
        for mu in 'ACDEFGHIKLMNPQRSTVWY':
            if mu == wt: continue
            rows.append((t, f'{wt}{g+1}{mu}', float(dist[tok.convert_tokens_to_ids(mu)] - dist[wid])))
with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['target', 'protein_variant', 'esm2_llr']); w.writerows(rows)
print('wrote', OUT, len(rows), 'rows', flush=True)
