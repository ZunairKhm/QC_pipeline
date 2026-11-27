#!/usr/bin/env python3
"""Group TRACS align outputs by individual and pack them into slurm array chunks.

Only within-individual comparisons are wanted, so combine+distance run per individual:
~9k sample pairs in total, versus millions if every sample sharing a reference were
compared. Each individual's samples are read once, so I/O is the same as a global run.

Inputs (from the run folder made by prepare_run.py):
  <run>/sample_individual.tsv  run -> individual
  <run>/sources.tsv            units.batch*.tsv globs (unit_id -> runs) and .sqfs images per source
  (and optionally results/tracs_align/ for samples not yet squashed)

Outputs (under <run>/):
  index.tsv              every aligned unit with its individual and where it lives
  chunks/chunk_NNNN.tsv  one line per unit: safe_id, individual, unit, source
                         (source = the .sqfs holding it, or the unsquashed sample dir)

Individuals already finished (per_individual/<safe_id>/units.txt matches their current
unit list) are skipped. With --reuse-from <dir>, finished individuals in <dir>/per_individual
whose unit list matches are copied in first (the within-individual distances don't depend on
frac_N or on other people, so they carry over between runs).
"""
import argparse
import csv
import math
import os
import re
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import sources as S

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run-dir", required=True)
p.add_argument("--reuse-from", default=None,
               help="earlier run folder (or legacy results/tracs_distance) whose per_individual/ to reuse")
p.add_argument("--unsquashed-dir", default="results/tracs_align",
               help="used only with --include-unsquashed")
p.add_argument("--include-unsquashed", action="store_true",
               help="also use sample dirs not yet packed into a .sqfs (only do this once "
                    "the batch that wrote them has finished, or dirs may be incomplete)")
p.add_argument("--units-per-chunk", type=int, default=150)
p.add_argument("--min-units", type=int, default=2, help="skip individuals with fewer units")
p.add_argument("--limit-individuals", type=int, default=0,
               help="only chunk the N smallest todo individuals (for a test run)")
args = p.parse_args()

outdir = os.path.realpath(args.run_dir)
srcs = S.run_sources(outdir)
S.read_config(outdir)  # must exist


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s.replace("::", "__"))


# --- unit -> individual -------------------------------------------------------------
with open(os.path.join(outdir, "sample_individual.tsv")) as fh:
    run_ind = {r["run"]: r["individual"] for r in csv.DictReader(fh, delimiter="\t")}
unit_runs, unit_src = S.unit_runs(srcs)

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
unit_path = {}  # unit -> .sqfs image, or the sample dir itself if not squashed
for s in srcs:
    imgs = S.images(s)
    if not imgs:
        print(f"WARNING: source {s['name']}: no image matches {s['squash']}", file=sys.stderr)
    for img in imgs:
        for name in S.list_units(img):
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


def finished(base, ind):
    d = os.path.join(base, "per_individual", safe_of[ind])
    try:
        done_units = open(os.path.join(d, "units.txt")).read().split()
    except FileNotFoundError:
        return False
    return (sorted(done_units) == sorted(u for u, _ in groups[ind])
            and os.path.exists(os.path.join(d, "distances.csv"))
            and os.path.exists(os.path.join(d, "combined_metadata.csv")))


eligible = [i for i in groups if len(groups[i]) >= args.min_units]
n_reused = 0
if args.reuse_from:
    reuse = os.path.realpath(args.reuse_from)
    if reuse == outdir:
        sys.exit("ERROR: --reuse-from is the run folder itself")
    for ind in eligible:
        if not finished(outdir, ind) and finished(reuse, ind):
            dst = os.path.join(outdir, "per_individual", safe_of[ind])
            shutil.rmtree(dst, ignore_errors=True)
            os.makedirs(dst)
            src = os.path.join(reuse, "per_individual", safe_of[ind])
            for f in ("distances.csv", "combined_metadata.csv", "combined.tar"):
                if os.path.exists(os.path.join(src, f)):
                    shutil.copy2(os.path.join(src, f), dst)
            shutil.copy2(os.path.join(src, "units.txt"), dst)  # last = done
            n_reused += 1

done = [i for i in eligible if finished(outdir, i)]
todo = [i for i in eligible if not finished(outdir, i)]
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
not_aligned = sum(1 for u in unit_ind if u not in unit_path)
print(f"aligned units found:         {len(unit_path)}  ({aligned_no_ind} have no individual, skipped)")
print(f"units with an individual but no TRACS output in any image: {not_aligned}")
per_src = defaultdict(lambda: [0, set()])
for ind in eligible:
    for u, _ in groups[ind]:
        per_src[unit_src.get(u, "?")][0] += 1
        per_src[unit_src.get(u, "?")][1].add(ind)
for s, (nu, inds) in sorted(per_src.items()):
    print(f"  {s}: {nu} units in {len(inds)} individuals with >= {args.min_units} units")
print(f"individuals with >= {args.min_units} units: {len(eligible)}  "
      f"(+{len(groups) - len(eligible)} with fewer, skipped)")
print(f"already done:                {len(done)}" + (f"  ({n_reused} copied from {args.reuse_from})" if args.reuse_from else ""))
print(f"to run:                      {len(todo)} individuals, {n_units} units, {n_chunks} chunks")
if n_chunks:
    print("\nsubmit with:\n  " + S.sbatch_cmd("tracs_post/tracs_combine_distance.sbatch", n_chunks, outdir, "logs"))
