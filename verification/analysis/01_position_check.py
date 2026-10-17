"""Check that the wild-type residue of every missense variant matches the reference FASTA."""
import pandas as pd
from common import *

rows = []
for name in targets():
    fa = read_fasta(name)
    if not fa:
        rows.append(dict(target=name, note='no fasta')); continue
    header, seq = fa[0]
    for src, df in (('csv', load_all_variants(name)), ('mutations', load_resistance(name))):
        if df is None: continue
        if src == 'csv':
            df = df[df['Type'] == 'Substitution - Missense']
        muts = df['AA Mutation'].dropna().unique()
        parsed = [parse_missense(m) for m in muts]
        ok = [p for p in parsed if p]
        match = sum(1 for w, pos, _ in ok if pos <= len(seq) and seq[pos-1] == w)
        beyond = sum(1 for w, pos, _ in ok if pos > len(seq))
        # best alternative offset, to detect numbering shifts
        best = (0, match)
        for off in range(-400, 401):
            if off == 0: continue
            c = sum(1 for w, pos, _ in ok if 0 < pos+off <= len(seq) and seq[pos+off-1] == w)
            if c > best[1]: best = (off, c)
        rows.append(dict(target=name, source=src, fasta_header=header[:60], fasta_len=len(seq),
                         n_fasta_records=len(fa), unique_aas=len(muts), clean_missense=len(ok),
                         unparsable=len(muts)-len(ok), wt_match=match,
                         match_frac=round(match/len(ok), 3) if ok else None,
                         pos_beyond_seq=beyond, best_offset=best[0], best_offset_match=best[1]))
out = pd.DataFrame(rows)
out.to_csv('01_position_check.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20)
print(out.drop(columns=['fasta_header']).to_string(index=False))
