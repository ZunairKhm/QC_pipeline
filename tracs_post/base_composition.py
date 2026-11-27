#!/usr/bin/env python3
"""Letter counts (A C G T N X other) of every TRACS align consensus, per sample and reference.

TRACS combine only records frac_N. But in the consensus `N` is low depth, and `X` is a
covered site where no allele passed the posterior-frequency threshold (reads, nothing
called). A sample can be 60% N + 40% X with no A/C/G/T at all and still pass a frac_N cut,
giving pairs with 0 sites considered (see README, "X sites"). This table gives frac_X and
the called share so samples can be filtered on what actually gets compared.

The table is a property of the images, not of a run, so it is shared, like the listing
cache: results/tracs_distance/.base_composition/<image key>.tsv, one per .sqfs, rebuilt when
the image is newer. Columns:
  unit reference length A C G T N X other frac_N frac_X frac_called
(reference is versioned, as in the fasta name; other = IUPAC ambiguity codes, lower case, etc.)

From QC_pipeline/:
  python3 tracs_post/base_composition.py build          # chunks for images without a table; prints sbatch
  ... submit, wait ...
  python3 tracs_post/base_composition.py combine        # parts -> <key>.tsv, re-run until all done
  python3 tracs_post/base_composition.py table          # cat every image's table to stdout
`count` is what the sbatch task runs: unsquashfs -cat streams the fastas of one chunk through
a counter, nothing is written to disk but the part file.
"""
import argparse
import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import sources as S

COMP_DIR = "results/tracs_distance/.base_composition"
COLS = ["unit", "reference", "length", "A", "C", "G", "T", "N", "X", "other",
        "frac_N", "frac_X", "frac_called"]

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("cmd", choices=["build", "count", "combine", "table"])
p.add_argument("--sources", default="tracs_post/sources.tsv")
p.add_argument("--comp-dir", default=COMP_DIR)
p.add_argument("--per-chunk", type=int, default=3000, help="fastas per array task (~1 min each)")
p.add_argument("--chunk", help="count: the chunk file to process")
p.add_argument("--processors", type=int, default=4, help="count: unsquashfs -processors")
args = p.parse_args()

comp = os.path.realpath(args.comp_dir)
chunk_dir, part_dir, log_dir = (os.path.join(comp, d) for d in ("chunks", "parts", "logs"))


def table_path(img):
    return os.path.join(comp, S.cache_key(img) + ".tsv")


def table_fresh(img):
    t = table_path(img)
    return os.path.exists(t) and os.path.getmtime(t) >= os.path.getmtime(img)


def all_images():
    imgs = []
    for s in S.read_sources(args.sources):
        imgs += S.images(s)
    return sorted(set(imgs))


# --- build: chunk files, one image per chunk ---------------------------------------------
if args.cmd == "build":
    for d in (chunk_dir, part_dir, log_dir):
        os.makedirs(d, exist_ok=True)
    for f in glob.glob(os.path.join(chunk_dir, "*.tsv")):
        os.remove(f)
    n_chunks, n_fastas = 0, 0
    for img in all_images():
        if table_fresh(img):
            continue
        key = S.cache_key(img)
        fastas = S.list_fastas(img)
        parts = sorted(glob.glob(os.path.join(part_dir, key + ".*.tsv")))
        if parts:  # stale parts from an older image
            if os.path.getmtime(parts[0]) < os.path.getmtime(img):
                for f in parts:
                    os.remove(f)
                parts = []
        done = {os.path.basename(f)[len(key) + 1:-4] for f in parts}
        k = 0
        for i in range(0, len(fastas), args.per_chunk):
            tag = f"{k:04d}"
            k += 1
            if tag in done:
                continue
            with open(os.path.join(chunk_dir, f"chunk_{n_chunks:04d}.tsv"), "w") as fh:
                for unit, ref in fastas[i:i + args.per_chunk]:
                    fh.write(f"{img}\t{key}\t{tag}\t{unit}\t{ref}\n")
            n_chunks += 1
            n_fastas += len(fastas[i:i + args.per_chunk])
        print(f"{key}: {len(fastas)} fastas, {k} chunks ({len(done)} already counted)")
    print(f"to count: {n_fastas} fastas in {n_chunks} chunks (~{n_fastas / 80 / 60:.0f} task-minutes)")
    if n_chunks:
        print(f"\nsubmit with:\n  sbatch --array=0-{n_chunks - 1}%40 --export=ALL,COMP_DIR={comp} \\\n"
              f"       --output={log_dir}/%A_%a.out tracs_post/base_composition.sbatch")

# --- count: one chunk, streamed through unsquashfs -cat ----------------------------------
elif args.cmd == "count":
    rows = [l.rstrip("\n").split("\t") for l in open(args.chunk)]
    img, key, tag = rows[0][:3]
    assert all(r[0] == img and r[2] == tag for r in rows), "a chunk holds one image"
    want = {(r[3], r[4]) for r in rows}
    paths = [f"{u}/{u}_posterior_counts_ref_{r}.fasta" for u, r in sorted(want)]
    out = os.path.join(part_dir, f"{key}.{tag}.tsv")
    os.makedirs(part_dir, exist_ok=True)
    proc = subprocess.Popen([S.UNSQUASHFS, "-cat", "-processors", str(args.processors), img, *paths],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    got = set()
    hdr = None
    with open(out + ".tmp", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for line in proc.stdout:
            if line[:1] == b">":
                hdr = line[1:].strip().decode()
                continue
            if hdr is None:
                sys.exit("ERROR: sequence before header")
            unit, ref = hdr.rsplit("_GC", 1)
            ref = "GC" + ref
            seq = line.rstrip(b"\r\n").upper()
            n = len(seq)
            c = {b: seq.count(b.encode()) for b in "ACGTNX"}
            other = n - sum(c.values())
            called = c["A"] + c["C"] + c["G"] + c["T"]
            fh.write("\t".join(map(str, [unit, ref, n, c["A"], c["C"], c["G"], c["T"], c["N"], c["X"], other,
                                         f"{c['N'] / n:.6f}", f"{c['X'] / n:.6f}", f"{called / n:.6f}"])) + "\n")
            got.add((unit, ref))
            hdr = None
    err = proc.stderr.read().decode()
    if proc.wait() != 0:
        sys.exit(f"ERROR: unsquashfs failed: {err[-2000:]}")
    if got != want:
        sys.exit(f"ERROR: {len(want - got)} fastas missing from output, {len(got - want)} unexpected; "
                 f"e.g. {sorted(want - got)[:3]} {sorted(got - want)[:3]}")
    os.replace(out + ".tmp", out)
    print(f"{key} chunk {tag}: {len(got)} fastas")

# --- combine: parts -> one table per image -----------------------------------------------
elif args.cmd == "combine":
    for img in all_images():
        key = S.cache_key(img)
        if table_fresh(img):
            print(f"{key}: done")
            continue
        n_expected = -(-len(S.list_fastas(img)) // args.per_chunk)
        parts = sorted(glob.glob(os.path.join(part_dir, key + ".*.tsv")))
        if len(parts) < n_expected:
            print(f"{key}: {len(parts)}/{n_expected} parts, waiting")
            continue
        n = 0
        with open(table_path(img) + ".tmp", "w") as fh:
            fh.write("\t".join(COLS) + "\n")
            for f in parts:
                with open(f) as ph:
                    next(ph)
                    for line in ph:
                        fh.write(line)
                        n += 1
        os.replace(table_path(img) + ".tmp", table_path(img))
        for f in parts:
            os.remove(f)
        print(f"{key}: {n} rows -> {table_path(img)}")

# --- table: everything, to stdout --------------------------------------------------------
elif args.cmd == "table":
    print("\t".join(COLS))
    for img in all_images():
        if not table_fresh(img):
            sys.exit(f"ERROR: no table for {img}; run build/combine")
        with open(table_path(img)) as fh:
            next(fh)
            sys.stdout.write(fh.read())
