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
With min_called > 0 (run_config.tsv or --min-called) samples are filtered on their called share
from base_composition.py, and each person's best-called sample is used; that needs no
within-individual run. Otherwise, if the within-individual run has finished, its combined_metadata.csv files give each
sample's frac_N per reference, and --max-frac-n drops poorly covered samples before
they use up the pair budget. Samples with unknown frac_N go after those with known
frac_N.

Needs index.tsv from build_chunks.py (run that first). The frac_N cut comes from the run
folder's run_config.tsv (max_frac_n) unless --max-frac-n is given.

Outputs (under <run>/background/):
  selection.tsv          ref, unit, individual, study, source
  chunks/chunk_NNNN.tsv  the same columns, grouped into array tasks by estimated cost
A reference is skipped if per_ref/<ref>.units already matches its selection, or with
--keep-done, if it has any finished output at all. With --reuse-from <dir> [<dir> ...], a reference whose
<dir>/background/per_ref/<ref>.units equals the new selection is copied in instead of rerun
(same samples -> same MSA -> same distances).
"""
import argparse
import csv
import glob
import math
import os
import random
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import sources as S

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run-dir", required=True)
p.add_argument("--reuse-from", nargs="+", default=[],
               help="earlier run folder(s) (or legacy frac_n_* dirs) whose background/per_ref to reuse; "
                    "checked in the order given")
p.add_argument("--max-pairs", type=int, default=10000, help="per reference")
p.add_argument("--max-frac-n", type=float, default=None,
               help="drop samples whose frac_N for the reference is above this (default: the run's "
                    "max_frac_n; needs the within-individual run's combined_metadata.csv)")
p.add_argument("--min-called", type=float, default=None,
               help="drop samples whose called share (A/C/G/T, base_composition.py) for the reference is "
                    "below this or unknown (default: the run's min_called; 0 = off). Needs no "
                    "within-individual run, so the background can start alongside it")
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
# Without a -D cut every pair is filtered: 101 samples (5,050 pairs) took 950 s at 16 threads
# and 2.0 CPU-hours (1.5 CPU-s per pair), so ~0.5 s per pair at the 4 threads a reference gets.
SEC_PER_REF, SEC_PER_FASTA, SEC_PER_PAIR = 5.0, 0.6, 0.004
SEC_PER_PAIR_UNCAPPED = 1.5  # run-2026-11-12: K=141 refs took a median 3.9 h (max 11.4 h) at 4 threads

outdir = os.path.realpath(args.run_dir)
cfg = S.read_config(outdir)
if args.max_frac_n is None:
    args.max_frac_n = float(cfg["max_frac_n"])
if args.min_called is None:
    args.min_called = float(cfg.get("min_called", 0))
if int(cfg.get("snp_threshold", 5000)) > 5000:  # pairs up to the run's -D are filtered
    SEC_PER_PAIR = SEC_PER_PAIR_UNCAPPED
called = S.load_called() if args.min_called > 0 else {}
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
ref_units = defaultdict(list)
images = sorted({src for _, src in index.values() if src.endswith(".sqfs")})
for img in images:
    for unit, ref in S.list_fastas(img):
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
if args.max_frac_n is not None and not frac_n and not called:
    print("WARNING: --max-frac-n given but no combined_metadata.csv found; not filtering", file=sys.stderr)

# --- selection ------------------------------------------------------------------------
study_of = lambda ind: ind.split("::")[0]
selection = {}
n_dropped_fracn = n_dropped_called = 0
for ref in sorted(ref_units):
    rng = random.Random(f"{args.seed}:{ref}")  # per-ref seed: stable if other refs change
    by_ind = defaultdict(list)
    for unit in ref_units[ref]:
        fn = frac_n.get((unit, ref))
        if args.max_frac_n is not None and fn is not None and fn > args.max_frac_n:
            n_dropped_fracn += 1
            continue
        if called and called.get((unit, ref.split(".")[0]), -1) < args.min_called:
            n_dropped_called += 1
            continue
        by_ind[index[unit][0]].append(unit)
    if len(by_ind) < 2:
        continue
    # one unit per individual: best known coverage (largest called share if filtering on it,
    # else lowest frac_N), else random
    pick = {}
    for ind, units in by_ind.items():
        known = [u for u in units if (u, ref) in frac_n]
        if called:
            pick[ind] = max(sorted(units), key=lambda u: called[(u, ref.split(".")[0])])
        else:
            pick[ind] = min(known, key=lambda u: frac_n[(u, ref)]) if known else rng.choice(sorted(units))
    # round-robin across studies; within a study, known-coverage individuals first
    by_study = defaultdict(list)
    for ind in sorted(pick):
        by_study[study_of(ind)].append(ind)
    for inds in by_study.values():
        rng.shuffle(inds)
        if not called:  # with the called filter every pick is known; frac_N must not change the selection
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
def is_done(ref, base=bgdir, keep_done=None):
    d = os.path.join(base, "per_ref")
    try:
        done = open(os.path.join(d, ref + ".units")).read().split()
    except FileNotFoundError:
        return False
    if not all(os.path.exists(os.path.join(d, ref + x)) for x in (".csv", ".fracn.tsv")):
        return False
    return (args.keep_done if keep_done is None else keep_done) or done == selection[ref]


n_reused = {}
for reuse in args.reuse_from:
    old_bg = os.path.join(os.path.realpath(reuse), "background")
    if old_bg == bgdir:
        sys.exit("ERROR: --reuse-from is the run folder itself")
    os.makedirs(os.path.join(bgdir, "per_ref"), exist_ok=True)
    n_reused[reuse] = 0
    for ref in selection:
        if not is_done(ref) and is_done(ref, old_bg, keep_done=False):
            for x in (".csv", ".fracn.tsv", ".units"):  # .units last = done
                shutil.copy2(os.path.join(old_bg, "per_ref", ref + x), os.path.join(bgdir, "per_ref"))
            n_reused[reuse] += 1


def cost(ref):
    k = len(selection[ref])
    return SEC_PER_REF + SEC_PER_FASTA * k + SEC_PER_PAIR * k * (k - 1) / 2


todo = [r for r in selection if not is_done(r)]
if args.limit_refs:
    todo = sorted(todo, key=cost, reverse=True)[: args.limit_refs]
budget = args.seconds_per_chunk
total = sum(cost(r) for r in todo)
n_chunks = min(len(todo), max(1, math.ceil(total / budget))) if todo else 0  # never more chunks than refs
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
if called:
    print(f"samples dropped for frac_called < {args.min_called:g} (or no count): {n_dropped_called}")
print(f"frac_N known for {len(frac_n)} sample/reference combinations")
for reuse, n in n_reused.items():
    print(f"copied unchanged references from {reuse}: {n}")
print(f"to run: {len(todo)} references in {n_chunks} chunks; estimated {total / 3600:.1f} "
      f"stream-hours (~{total / 3600 / 4 * 16:.0f} core-hours at 4 refs x 4 threads per 16-core task)")
if n_chunks:
    print("\nsubmit with:\n  " + S.sbatch_cmd("tracs_post/tracs_background.sbatch", n_chunks, outdir,
                                             "background/logs"))
