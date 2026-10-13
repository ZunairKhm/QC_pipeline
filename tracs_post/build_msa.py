#!/usr/bin/env python3
"""Write one reference's MSA from per-sample TRACS posterior_counts fastas.

Same output as tracs combine's write_alignment (">sample" + one-line sequence, gzipped)
but for a single reference, so the background run only reads the fastas it uses.
Also writes each sample's frac_N and length. Empty or zero-length sequences are left
out: tracs distance segfaults on them.

usage: build_msa.py <out.fasta.gz> <out.fracn.tsv> <sample>=<fasta> ...
Kept Python 3.6-compatible (the system python3 on Setonix).
"""
import gzip
import sys

out_fa, out_fn = sys.argv[1], sys.argv[2]
n_written = 0
with gzip.open(out_fa, "wt", compresslevel=1) as fa, open(out_fn, "w") as fn:
    fn.write("sample\tlength\tfrac_N\n")
    for arg in sys.argv[3:]:
        sample, path = arg.split("=", 1)
        seq, n_records = [], 0
        try:
            with open(path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        n_records += 1
                    else:
                        seq.append(line.strip())
        except FileNotFoundError:
            sys.stderr.write("WARNING: missing {}\n".format(path))
            continue
        if n_records != 1:
            sys.exit("ERROR: {} has {} sequences, expected 1".format(path, n_records))
        seq = "".join(seq)
        if not seq:
            sys.stderr.write("WARNING: empty sequence in {}\n".format(path))
            continue
        fa.write(">{}\n{}\n".format(sample, seq))
        fn.write("{}\t{}\t{:.6f}\n".format(sample, len(seq), seq.upper().count("N") / len(seq)))
        n_written += 1
print(n_written)
