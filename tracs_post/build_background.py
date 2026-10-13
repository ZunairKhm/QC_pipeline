#!/usr/bin/env python3
"""Pick between-person background samples per reference and pack them into array chunks.

For each reference, at most K samples are chosen so that all K*(K-1)/2 pairs stay
within --max-pairs (10000 -> K=141). tracs distance compares every pair in an MSA, so
one MSA of K samples gives the most samples for the pair budget. Every pair is between
different people:
  - one sample per individual (units with no individual are left out, since they
    could be the same person)
  - individuals are taken round-robin across studies (the "Study_YEAR" prefix of the
    individual), so as many studies as possible are represented. Most pairs are then
    between studies; same-study, different-person pairs are kept but flagged at merge
    time, because some studies sampled families.
If the within-individual run has finished, its combined_metadata.csv files give each
sample's frac_N per reference, and --max-frac-n drops poorly covered samples before
they use up the pair budget. Samples with unknown frac_N go after those with known
frac_N.

Needs index.tsv from build_chunks.py (run that first).

Outputs (under <outdir>/background/):
  selection.tsv          ref, unit, individual, study, source
  chunks/chunk_NNNN.tsv  the same columns, grouped into array tasks by estimated cost
A reference is skipped if per_ref/<ref>.units already matches its selection, or with
--keep-done, if it has any finished output at all.
"""
import argparse
import csv
import glob
import math
import os
import random
import shutil
import subprocess
import sys
from collections import defaultdict

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--outdir", default="results/tracs_distance")
p.add_argument("--max-pairs", type=int, default=10000, help="per reference")
p.add_argument("--max-frac-n", type=float, default=0.75,
               help="drop samples whose frac_N for the reference is above this (needs the "
                    "within-individual run's combined_metadata.csv)")
p.add_argument("--seed", type=int, default=1)
p.add_argument("--seconds-per-chunk", type=float, default=4 * 3600,
               help="estimated single-stream seconds per array task; the task runs "
                    "PARALLEL (default 4) references at once, so wall time is ~1/4 of this")
p.add_argument("--keep-done", action="store_true",
               help="treat every reference with finished output as done, even if new samples "
                    "would change its selection (only new/unfinished references run)")
p.add_argument("--limit-refs", type=int, default=0, help="only the N largest refs (for a test run)")
args = p.parse_args()

# Per-reference cost in seconds, fitted to test runs with --filter -D 5000 at 16 threads:
# K=141 refs took 59-159s (staging 40-117s, distance 19-42s), K=40 16s, K=10 6s.
# Staging the fastas out of squashfs dominates; distance is small once -D prunes pairs.
SEC_PER_REF, SEC_PER_FASTA, SEC_PER_PAIR = 5.0, 0.6, 0.004

outdir = os.path.realpath(args.outdir)
bgdir = os.path.join(outdir, "background")
os.makedirs(os.path.join(bgdir, "logs"), exist_ok=True)
K = max(k for k in range(2, 10000) if k * (k - 1) // 2 <= args.max_pairs)

# --- unit -> individual, source (from build_chunks.py) --------------------------------
index = {}
with open(os.path.join(outdir, "index.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        index[r["unit"]] = (r["individual"], r["source"])
if not index:
    sys.exit("ERROR: index.tsv is empty - run build_chunks.py first")

# --- which references each unit has a posterior fasta for -----------------------------
unsquashfs = "/usr/bin/unsquashfs" if os.access("/usr/bin/unsquashfs", os.X_OK) else "unsquashfs"
listing_dir = os.path.join(outdir, ".sqfs_listing")
os.makedirs(listing_dir, exist_ok=True)


def fasta_listing(src):
    """(unit, ref) pairs held by one .sqfs image, cached until the image changes."""
    stem = os.path.basename(src)[: -len(".sqfs")]
    cache = os.path.join(listing_dir, stem + ".fastas.tsv")
    if not os.path.exists(cache) or os.path.getmtime(cache) < os.path.getmtime(src):
        print(f"listing files in {src} ...", file=sys.stderr)
        out = subprocess.run([unsquashfs, "-l", src], check=True, capture_output=True, text=True).stdout
        with open(cache, "w") as fh:
            for line in out.splitlines():
                parts = line.split("/")
                if len(parts) == 3 and "_posterior_counts_ref_" in parts[2] and parts[2].endswith(".fasta"):
                    ref = parts[2].split("_posterior_counts_ref_", 1)[1][: -len(".fasta")]
                    fh.write(f"{parts[1]}\t{ref}\n")
    return [tuple(l.split()) for l in open(cache)]


ref_units = defaultdict(list)
images = sorted({src for _, src in index.values() if src.endswith(".sqfs")})
for img in images:
    for unit, ref in fasta_listing(img):
        if unit in index and index[unit][1] == img:
            ref_units[ref].append(unit)
for unit, (_, src) in index.items():
    if not src.endswith(".sqfs"):  # unsquashed sample dir
        for f in os.listdir(src):
            if "_posterior_counts_ref_" in f and f.endswith(".fasta"):
                ref_units[f.split("_posterior_counts_ref_", 1)[1][: -len(".fasta")]].append(unit)

# --- frac_N from the within-individual run, if available ------------------------------
frac_n = {}
for f in glob.glob(os.path.join(outdir, "per_individual", "*", "combined_metadata.csv")):
    with open(f) as fh:
        next(fh)
        for line in fh:
            c = line.split(",")
            if c[9] not in ("NA", ""):
                frac_n[(c[0], c[1])] = float(c[9])
if args.max_frac_n is not None and not frac_n:
    print("WARNING: --max-frac-n given but no combined_metadata.csv found; not filtering", file=sys.stderr)

# --- selection ------------------------------------------------------------------------
study_of = lambda ind: ind.split("::")[0]
selection = {}
n_dropped_fracn = 0
for ref in sorted(ref_units):
    rng = random.Random(f"{args.seed}:{ref}")  # per-ref seed: stable if other refs change
    by_ind = defaultdict(list)
    for unit in ref_units[ref]:
        fn = frac_n.get((unit, ref))
        if args.max_frac_n is not None and fn is not None and fn > args.max_frac_n:
            n_dropped_fracn += 1
            continue
        by_ind[index[unit][0]].append(unit)
    if len(by_ind) < 2:
        continue
    # one unit per individual: best known coverage, else random
    pick = {}
    for ind, units in by_ind.items():
        known = [u for u in units if (u, ref) in frac_n]
        pick[ind] = min(known, key=lambda u: frac_n[(u, ref)]) if known else rng.choice(sorted(units))
    # round-robin across studies; within a study, known-coverage individuals first
    by_study = defaultdict(list)
    for ind in sorted(pick):
        by_study[study_of(ind)].append(ind)
    for inds in by_study.values():
        rng.shuffle(inds)
        inds.sort(key=lambda i: (pick[i], ref) not in frac_n)
    studies = sorted(by_study)
    rng.shuffle(studies)
    chosen = []
    while len(chosen) < K and any(by_study[s] for s in studies):
        for s in studies:
            if by_study[s] and len(chosen) < K:
                chosen.append(by_study[s].pop(0))
    selection[ref] = sorted(pick[i] for i in chosen)

with open(os.path.join(bgdir, "selection.tsv"), "w") as fh:
    fh.write("ref\tunit\tindividual\tstudy\tsource\n")
    for ref, units in selection.items():
        for u in units:
            ind, src = index[u]
            fh.write(f"{ref}\t{u}\t{ind}\t{study_of(ind)}\t{src}\n")


# --- chunks ---------------------------------------------------------------------------
def is_done(ref):
    d = os.path.join(bgdir, "per_ref")
    try:
        done = open(os.path.join(d, ref + ".units")).read().split()
    except FileNotFoundError:
        return False
    if not os.path.exists(os.path.join(d, ref + ".csv")):
        return False
    return args.keep_done or done == selection[ref]


def cost(ref):
    k = len(selection[ref])
    return SEC_PER_REF + SEC_PER_FASTA * k + SEC_PER_PAIR * k * (k - 1) / 2


todo = [r for r in selection if not is_done(r)]
if args.limit_refs:
    todo = sorted(todo, key=cost, reverse=True)[: args.limit_refs]
budget = args.seconds_per_chunk
total = sum(cost(r) for r in todo)
n_chunks = max(1, math.ceil(total / budget)) if todo else 0
chunks, loads = [[] for _ in range(n_chunks)], [0] * n_chunks
for ref in sorted(todo, key=cost, reverse=True):
    k = loads.index(min(loads))
    chunks[k].append(ref)
    loads[k] += cost(ref)

chunk_dir = os.path.join(bgdir, "chunks")
shutil.rmtree(chunk_dir, ignore_errors=True)
os.makedirs(chunk_dir)
for k, refs in enumerate(chunks):
    with open(os.path.join(chunk_dir, f"chunk_{k:04d}.tsv"), "w") as fh:
        for ref in refs:
            for u in selection[ref]:
                ind, src = index[u]
                fh.write(f"{ref}\t{u}\t{ind}\t{study_of(ind)}\t{src}\n")

sizes = [len(v) for v in selection.values()]
pairs = sum(k * (k - 1) // 2 for k in sizes)
cross = 0
for ref, units in selection.items():
    st = defaultdict(int)
    for u in units:
        st[study_of(index[u][0])] += 1
    n = len(units)
    cross += n * (n - 1) // 2 - sum(c * (c - 1) // 2 for c in st.values())
print(f"K = {K} samples per reference (<= {args.max_pairs} pairs)")
print(f"references with >= 2 individuals: {len(selection)}  ({sum(s == K for s in sizes)} capped at K)")
print(f"selected fastas: {sum(sizes)};  pairs: {pairs}  ({cross / max(pairs, 1):.0%} between studies)")
if args.max_frac_n is not None:
    print(f"samples dropped for frac_N > {args.max_frac_n}: {n_dropped_fracn}")
print(f"frac_N known for {len(frac_n)} sample/reference combinations")
print(f"to run: {len(todo)} references in {n_chunks} chunks; estimated {total / 3600:.1f} "
      f"stream-hours (~{total / 3600 / 4 * 16:.0f} core-hours at 4 refs x 4 threads per 16-core task)")
if n_chunks:
    print(f"\nsubmit with:\n  sbatch --array=0-{n_chunks - 1}%40 tracs_post/tracs_background.sbatch")
