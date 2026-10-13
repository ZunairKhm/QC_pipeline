#!/usr/bin/env python3
"""Group TRACS align outputs by individual and pack them into slurm array chunks.

Only within-individual comparisons are wanted, so combine+distance run per individual:
~9k sample pairs in total, versus millions if every sample sharing a reference were
compared. Each individual's samples are read once, so I/O is the same as a global run.

Inputs:
  sample_individual.tsv        run -> individual
  results/metadata/units.*.tsv unit_id (= tracs sample dir name) -> runs
  results/tracs_align_squash/*.sqfs  (and optionally results/tracs_align/ for
                                      samples not yet squashed)

Outputs (under --outdir):
  index.tsv              every aligned unit with its individual and where it lives
  chunks/chunk_NNNN.tsv  one line per unit: safe_id, individual, unit, source
                         (source = the .sqfs holding it, or the unsquashed sample dir)

Individuals already finished (per_individual/<safe_id>/units.txt matches their current
unit list) are skipped, so re-running this after more samples land only redoes the
individuals whose sample set changed.
"""
import argparse
import csv
import glob
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--individuals", default="sample_individual.tsv")
p.add_argument("--units", default="results/metadata/units.batch*.tsv", help="glob")
p.add_argument("--squash-dir", default="results/tracs_align_squash")
p.add_argument("--unsquashed-dir", default="results/tracs_align",
               help="used only with --include-unsquashed")
p.add_argument("--include-unsquashed", action="store_true",
               help="also use sample dirs not yet packed into a .sqfs (only do this once "
                    "the batch that wrote them has finished, or dirs may be incomplete)")
p.add_argument("--outdir", default="results/tracs_distance")
p.add_argument("--units-per-chunk", type=int, default=150)
p.add_argument("--min-units", type=int, default=2, help="skip individuals with fewer units")
p.add_argument("--limit-individuals", type=int, default=0,
               help="only chunk the N smallest todo individuals (for a test run)")
args = p.parse_args()

outdir = os.path.realpath(args.outdir)
os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)  # sbatch --output lands here


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s.replace("::", "__"))


# --- unit -> individual -------------------------------------------------------------
run_ind = {r["run"]: r["individual"]
           for r in csv.DictReader(open(args.individuals), delimiter="\t")}
unit_runs = defaultdict(set)
for f in sorted(glob.glob(args.units)):
    for r in csv.DictReader(open(f), delimiter="\t"):
        unit_runs[r["unit_id"]].add(r["run"])

unit_ind, no_ind = {}, []
for unit, runs in unit_runs.items():
    inds = {run_ind[r] for r in runs if r in run_ind}
    if len(inds) > 1:
        sys.exit(f"ERROR: unit {unit} maps to several individuals: {sorted(inds)}")
    if inds:
        unit_ind[unit] = inds.pop()
    else:
        no_ind.append(unit)

# --- where each unit's tracs output lives -------------------------------------------
unsquashfs = "/usr/bin/unsquashfs" if os.access("/usr/bin/unsquashfs", os.X_OK) else "unsquashfs"
listing_dir = os.path.join(outdir, ".sqfs_listing")
os.makedirs(listing_dir, exist_ok=True)

unit_path = {}  # unit -> .sqfs image, or the sample dir itself if not squashed
for img in sorted(glob.glob(os.path.join(args.squash_dir, "*.sqfs"))):
    img = os.path.realpath(img)
    stem = os.path.basename(img)[: -len(".sqfs")]
    # listing an image takes ~15s, so cache it until the image changes
    cache = os.path.join(listing_dir, stem + ".txt")
    if not os.path.exists(cache) or os.path.getmtime(cache) < os.path.getmtime(img):
        print(f"listing {img} ...", file=sys.stderr)
        out = subprocess.run([unsquashfs, "-l", "-max-depth", "1", img],
                             check=True, capture_output=True, text=True).stdout
        names = [l.split("/", 1)[1] for l in out.splitlines() if l.startswith("squashfs-root/")]
        with open(cache, "w") as fh:
            fh.write("\n".join(names) + "\n")
    for name in open(cache).read().split():
        unit_path.setdefault(name, img)

if args.include_unsquashed and os.path.isdir(args.unsquashed_dir):
    src = os.path.realpath(args.unsquashed_dir)
    for name in os.listdir(src):
        if os.path.isdir(os.path.join(src, name)):
            unit_path.setdefault(name, os.path.join(src, name))  # squashed copy wins

# --- per-individual groups ----------------------------------------------------------
groups = defaultdict(list)
for unit, path in unit_path.items():
    if unit in unit_ind:
        groups[unit_ind[unit]].append((unit, path))

safe_of, seen = {}, {}
for ind in groups:
    s = safe_name(ind)
    if s in seen:
        sys.exit(f"ERROR: individuals {seen[s]!r} and {ind!r} both sanitise to {s}")
    seen[s] = ind
    safe_of[ind] = s

with open(os.path.join(outdir, "index.tsv"), "w") as fh:
    fh.write("unit\tindividual\tsafe_id\tsource\n")
    for ind in sorted(groups):
        for unit, path in sorted(groups[ind]):
            fh.write(f"{unit}\t{ind}\t{safe_of[ind]}\t{path}\n")


def is_done(ind):
    d = os.path.join(outdir, "per_individual", safe_of[ind])
    try:
        done_units = open(os.path.join(d, "units.txt")).read().split()
    except FileNotFoundError:
        return False
    return (sorted(done_units) == sorted(u for u, _ in groups[ind])
            and os.path.exists(os.path.join(d, "distances.csv")))


eligible = [i for i in groups if len(groups[i]) >= args.min_units]
todo = [i for i in eligible if not is_done(i)]
if args.limit_individuals:
    todo = sorted(todo, key=lambda i: len(groups[i]))[: args.limit_individuals]

# --- pack into chunks (largest first onto the lightest chunk) -----------------------
chunk_dir = os.path.join(outdir, "chunks")
shutil.rmtree(chunk_dir, ignore_errors=True)
os.makedirs(chunk_dir)
n_units = sum(len(groups[i]) for i in todo)
n_chunks = max(1, math.ceil(n_units / args.units_per_chunk)) if todo else 0
chunks = [[] for _ in range(n_chunks)]
loads = [0] * n_chunks
for ind in sorted(todo, key=lambda i: -len(groups[i])):
    k = loads.index(min(loads))
    chunks[k].append(ind)
    loads[k] += len(groups[ind])
for k, inds in enumerate(chunks):
    with open(os.path.join(chunk_dir, f"chunk_{k:04d}.tsv"), "w") as fh:
        for ind in inds:
            for unit, path in sorted(groups[ind]):
                fh.write(f"{safe_of[ind]}\t{ind}\t{unit}\t{path}\n")

aligned_no_ind = sum(1 for u in unit_path if u not in unit_ind)
print(f"aligned units found:         {len(unit_path)}  ({aligned_no_ind} have no individual, skipped)")
print(f"individuals with >= {args.min_units} units: {len(eligible)}  "
      f"(+{len(groups) - len(eligible)} with fewer, skipped)")
print(f"already done:                {len(eligible) - len([i for i in eligible if not is_done(i)])}")
print(f"to run:                      {len(todo)} individuals, {n_units} units, {n_chunks} chunks")
if n_chunks:
    print(f"\nsubmit with:\n  sbatch --array=0-{n_chunks - 1}%50 tracs_post/tracs_combine_distance.sbatch")
