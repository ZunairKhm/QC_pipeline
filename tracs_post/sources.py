"""Shared helpers for the tracs_post build scripts: data sources, run folders, image listings.

A run folder (results/tracs_distance/run-YYYY-MM-DD/) holds everything one run needs and writes;
prepare_run.py creates it with a copy of tracs_post/sources.tsv and run_config.tsv, and every
later step reads those two files from the run folder, so a run is reproducible on its own.

sources.tsv, one line per data source (tab-separated, '#' lines ignored):
  name        label, written to sample_individual.tsv
  units_glob  units.batch*.tsv files the pipeline wrote for it (unit_id -> run)
  squash      glob of its .sqfs image(s) holding the TRACS align sample dirs
  metadata    cmd     -> individuals/times from curatedMetagenomicData (export_cmd_times.R)
              metalog -> from metalog human_extended_wide (joined via ENA sample accession)
Paths are relative to QC_pipeline/ (the directory everything is run from).
"""
import csv
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

SQUASH_ROOT = "results/tracs_align_squash"
LISTING_DIR = "results/tracs_distance/.sqfs_listing"  # shared: describes images, not a run
UNSQUASHFS = "/usr/bin/unsquashfs" if os.access("/usr/bin/unsquashfs", os.X_OK) else "unsquashfs"


def read_sources(path):
    with open(path) as fh:
        rows = [l for l in fh if l.strip() and not l.startswith("#")]
    srcs = list(csv.DictReader(rows, delimiter="\t"))
    for s in srcs:
        if s.get("metadata") not in ("cmd", "metalog"):
            sys.exit(f"ERROR: {path}: source {s.get('name')!r} has metadata {s.get('metadata')!r} "
                     "(expected cmd or metalog)")
    names = [s["name"] for s in srcs]
    if len(set(names)) != len(names):
        sys.exit(f"ERROR: {path}: duplicate source names")
    return srcs


def run_sources(run_dir):
    p = os.path.join(run_dir, "sources.tsv")
    if not os.path.exists(p):
        sys.exit(f"ERROR: {p} missing - run tracs_post/prepare_run.py first")
    return read_sources(p)


def read_config(run_dir):
    p = os.path.join(run_dir, "run_config.tsv")
    if not os.path.exists(p):
        sys.exit(f"ERROR: {p} missing - run tracs_post/prepare_run.py first")
    with open(p) as fh:
        return {r["key"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}


def units_files(src):
    return sorted(glob.glob(src["units_glob"]))


def images(src):
    return sorted({os.path.realpath(p) for p in glob.glob(src["squash"])})


def unit_runs(srcs):
    """unit -> set of runs, and unit -> source name, over every source's units files."""
    runs, unit_src = defaultdict(set), {}
    for s in srcs:
        for f in units_files(s):
            with open(f) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    runs[r["unit_id"]].add(r["run"])
                    unit_src.setdefault(r["unit_id"], s["name"])
    return runs, unit_src


def cache_key(img):
    """Unique listing-cache name for an image: its path under SQUASH_ROOT with / -> __
    (tracs_align_batch00.sqfs -> tracs_align_batch00, yaffe/tracs_align.sqfs -> yaffe__tracs_align)."""
    img = os.path.realpath(img)
    root = os.path.realpath(SQUASH_ROOT)
    rel = os.path.relpath(img, root) if img.startswith(root + os.sep) else img.lstrip(os.sep)
    rel = rel[: -len(".sqfs")] if rel.endswith(".sqfs") else rel
    return re.sub(r"[^A-Za-z0-9_.-]", "_", rel.replace(os.sep, "__"))


def _cached(img, suffix, args):
    os.makedirs(LISTING_DIR, exist_ok=True)
    cache = os.path.join(LISTING_DIR, cache_key(img) + suffix)
    if not os.path.exists(cache) or os.path.getmtime(cache) < os.path.getmtime(img):
        print(f"listing {img} ...", file=sys.stderr)
        out = subprocess.run([UNSQUASHFS, "-l", *args, img], check=True, capture_output=True, text=True).stdout
        return cache, out
    return cache, None


def list_units(img):
    """Top-level sample dirs in an image (cached until the image changes)."""
    cache, out = _cached(img, ".txt", ["-max-depth", "1"])
    if out is not None:
        names = [l.split("/", 1)[1] for l in out.splitlines() if l.startswith("squashfs-root/")]
        with open(cache, "w") as fh:
            fh.write("\n".join(names) + "\n")
    return open(cache).read().split()


def list_fastas(img):
    """(unit, ref) for every posterior_counts fasta in an image (cached until the image changes)."""
    cache, out = _cached(img, ".fastas.tsv", [])
    if out is not None:
        with open(cache, "w") as fh:
            for line in out.splitlines():
                parts = line.split("/")
                if len(parts) == 3 and "_posterior_counts_ref_" in parts[2] and parts[2].endswith(".fasta"):
                    ref = parts[2].split("_posterior_counts_ref_", 1)[1][: -len(".fasta")]
                    fh.write(f"{parts[1]}\t{ref}\n")
    return [tuple(l.split()) for l in open(cache)]


def sbatch_cmd(script, n_tasks, run_dir, log_subdir, throttle=40, extra_export=""):
    """The sbatch line to print: run dir exported and logs inside the run dir."""
    run_dir = os.path.realpath(run_dir)
    log_dir = os.path.join(run_dir, log_subdir)
    os.makedirs(log_dir, exist_ok=True)
    arr = f"0-{n_tasks - 1}" + (f"%{throttle}" if throttle else "")
    return (f"sbatch --array={arr} --export=ALL,RUN_DIR={run_dir}{extra_export} \\\n"
            f"       --output={log_dir}/%A_%a.out {script}")


COMP_DIR = "results/tracs_distance/.base_composition"  # base_composition.py output, shared like LISTING_DIR


def load_called(comp_dir=COMP_DIR):
    """(unit, unversioned reference) -> frac_called, from every image's base_composition table."""
    import glob
    out = {}
    files = glob.glob(os.path.join(comp_dir, "*.tsv"))
    if not files:
        sys.exit(f"ERROR: no tables in {comp_dir}; run tracs_post/base_composition.py build/combine")
    for f in files:
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out[(r["unit"], r["reference"].split(".")[0])] = float(r["frac_called"])
    return out
