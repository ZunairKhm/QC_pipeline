#!/usr/bin/env python3
"""Regroup the within-individual distances by reference, with sampling times, for mxsure.

per_individual/<id>/distances.csv has every within-person pair on every reference; this
writes the same rows split by reference into per_reference/<ref>.csv, named with the
versioned accession like background/per_ref/<ref>.csv, so the two can be paired up per
reference. The columns are tracs distance's own, with "date difference" filled in
(days, from <run>/sample_times.tsv: cMD via export_cmd_times.R, metalog via prepare_run.py),
so a file can go straight into
mxsure_input_tracs(), plus:
  individual, study, body_site_A/B, time_source, frac_N_A/B
  longitudinal   TRUE if the pair counts as a longitudinal comparison for mxsure:
                 sampled on different days (time > 0, same time scale), same body site,
                 frac_N <= --max-frac-n for both samples, sites considered >= --min-sites,
                 raw SNP distance <= --snp-threshold (the runs' -D, as in merge_distances.py)
Time is days_from_first_collection; where a study only has infant_age (days of life),
that is used instead, but only if both samples have it.

per_reference/summary.tsv, one row per reference: pairs, longitudinal pairs and the
individuals they come from, the background's state, and `eligible` (>= --min-comparisons
longitudinal pairs and a finished background). mxsure_rates.R runs on the eligible ones.

Run from QC_pipeline/ after the within-individual run. --max-frac-n and --min-sites default to the
run's run_config.tsv. Always rebuilt from everything on disk, so re-run it once more of
the background has finished.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import sources as S

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run-dir", required=True)
p.add_argument("--times", default=None, help="default: <run>/sample_times.tsv")
p.add_argument("--max-frac-n", type=float, default=None,
               help="longitudinal pairs need frac_N <= this in both samples (default: the run's max_frac_n)")
p.add_argument("--min-sites", type=int, default=None,
               help="longitudinal pairs need sites considered >= this (default: the run's min_sites, else 0)")
p.add_argument("--snp-threshold", type=int, default=None,
               help="raw SNP distance cap (default: the run's snp_threshold, i.e. its -D)")
p.add_argument("--min-called", type=float, default=None,
               help="longitudinal pairs need frac_called >= this in both samples (default: the run's "
                    "min_called, else 0)")
p.add_argument("--min-comparisons", type=int, default=30,
               help="longitudinal pairs needed for a reference to be eligible")
args = p.parse_args()
args.outdir = args.run_dir
cfg = S.read_config(args.outdir)
if args.snp_threshold is None:
    args.snp_threshold = int(cfg.get("snp_threshold", 5000))
if args.max_frac_n is None:
    args.max_frac_n = float(cfg["max_frac_n"])
if args.min_sites is None:
    args.min_sites = int(cfg.get("min_sites", 0))
if args.min_called is None:
    args.min_called = float(cfg.get("min_called", 0))
called = S.load_called()

unver = lambda acc: acc.split(".")[0]
per_ref_dir = os.path.join(args.outdir, "per_reference")
bg_dir = os.path.join(args.outdir, "background", "per_ref")

# --- unit -> runs -> sampling time ------------------------------------------------------
unit_runs = {u: sorted(rs) for u, rs in S.unit_runs(S.run_sources(args.outdir))[0].items()}

run_time = {}
with open(args.times or os.path.join(args.outdir, "sample_times.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        run_time[r["run"]] = r

na = lambda v: v in ("", "NA", None)
unit_ind = {}
with open(os.path.join(args.outdir, "index.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        unit_ind[r["unit"]] = r["individual"]


def unit_info(unit):
    """metadata row for a unit (its first run with metadata; a unit is one biological sample)."""
    for run in unit_runs.get(unit, [unit]):
        if run in run_time:
            return run_time[run]
    return None


def pair_time(a, b):
    """(days apart, which time scale) or (None, reason)."""
    ia, ib = unit_info(a), unit_info(b)
    if ia is None or ib is None:
        return None, "no_metadata"
    for col in ("days_from_first_collection", "infant_age"):
        if not na(ia[col]) and not na(ib[col]):
            return abs(float(ia[col]) - float(ib[col])), col
    return None, "no_time"


# --- frac_N, versioned accession and species per reference ----------------------------
frac_n, versioned, species = {}, {}, {}
done = [d for d in sorted(glob.glob(os.path.join(args.outdir, "per_individual", "*")))
        if os.path.exists(os.path.join(d, "units.txt"))]
for d in done:
    with open(os.path.join(d, "combined_metadata.csv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split(",")
            acc = unver(f[1])
            frac_n[(f[0], acc)] = f[9]
            versioned.setdefault(acc, f[1])
            species.setdefault(acc, ",".join(f[10:]))

# --- regroup ----------------------------------------------------------------------------
TRACS_COLS = ["sampleA", "sampleB", "date difference", "SNP distance", "transmission distance",
              "expected K", "filtered SNP distance", "sites considered", "MSA file"]
EXTRA = ["individual", "study", "body_site_A", "body_site_B", "time_source", "frac_N_A", "frac_N_B",
         "frac_called_A", "frac_called_B", "longitudinal"]
rows = defaultdict(list)
untimed = defaultdict(int)
for d in done:
    with open(os.path.join(d, "distances.csv")) as fh:
        for r in csv.DictReader(fh):
            a, b, acc = r["sampleA"], r["sampleB"], r["MSA file"]
            t, src = pair_time(a, b)
            if t is None:
                untimed[src] += 1
            ia, ib = unit_info(a) or {}, unit_info(b) or {}
            fa, fb = frac_n.get((a, acc), "NA"), frac_n.get((b, acc), "NA")
            ca, cb = called.get((a, unver(acc))), called.get((b, unver(acc)))
            sites = int(r["sites considered"])
            longitudinal = (
                t is not None and t > 0
                and ia.get("body_site", "NA") == ib.get("body_site", "NA")
                and not na(fa) and not na(fb)
                and float(fa) <= args.max_frac_n and float(fb) <= args.max_frac_n
                and sites >= args.min_sites
                and ca is not None and cb is not None
                and ca >= args.min_called and cb >= args.min_called
                and not na(r["filtered SNP distance"])
                and int(r["SNP distance"]) <= args.snp_threshold
            )
            ind = unit_ind.get(a, "NA")
            r["date difference"] = "NA" if t is None else f"{t:g}"
            rows[acc].append([r[c] for c in TRACS_COLS] + [
                ind, ind.split("::")[0], ia.get("body_site", "NA"), ib.get("body_site", "NA"),
                src if t is not None else "NA", fa, fb,
                "NA" if ca is None else f"{ca:.4f}", "NA" if cb is None else f"{cb:.4f}",
                "TRUE" if longitudinal else "FALSE"])

os.makedirs(per_ref_dir, exist_ok=True)
for f in glob.glob(os.path.join(per_ref_dir, "*.csv")):
    os.remove(f)  # rebuilt from scratch
summary = []
for acc in sorted(rows):
    ref = versioned.get(acc, acc)
    with open(os.path.join(per_ref_dir, f"{ref}.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(TRACS_COLS + EXTRA)
        w.writerows(rows[acc])
    lon = [x for x in rows[acc] if x[-1] == "TRUE"]
    bg = os.path.join(bg_dir, ref)
    bg_done = os.path.exists(bg + ".units") and os.path.exists(bg + ".csv")
    bg_samples = sum(1 for _ in open(bg + ".fracn.tsv")) - 1 if bg_done else 0
    bg_pairs = sum(1 for _ in open(bg + ".csv")) - 1 if bg_done else 0
    summary.append([
        ref, species.get(acc, "NA"), len(rows[acc]), len({x[9] for x in rows[acc]}),
        len(lon), len({x[9] for x in lon}),
        f"{max(float(x[2]) for x in lon):g}" if lon else "NA",
        "TRUE" if bg_done else "FALSE", bg_samples, bg_pairs,
        "TRUE" if len(lon) >= args.min_comparisons and bg_done and bg_pairs > 0 else "FALSE",
    ])

summary.sort(key=lambda s: -s[4])
with open(os.path.join(per_ref_dir, "summary.tsv"), "w") as fh:
    fh.write("\t".join(["reference", "species", "n_pairs", "n_individuals", "n_longitudinal",
                        "n_longitudinal_individuals", "max_time_days", "background_done",
                        "background_samples", "background_pairs", "eligible"]) + "\n")
    for s in summary:
        fh.write("\t".join(map(str, s)) + "\n")

n_pairs = sum(len(v) for v in rows.values())
enough = [s for s in summary if s[4] >= args.min_comparisons]
print(f"max_frac_n = {args.max_frac_n:g}, min_sites = {args.min_sites}, min_called = {args.min_called:g}")
print(f"{len(done)} individuals, {n_pairs} pairs over {len(rows)} references -> {per_ref_dir}/")
print(f"  pairs without a sampling time: {sum(untimed.values())} {dict(untimed)}")
print(f"  references with >= {args.min_comparisons} longitudinal pairs: {len(enough)}; eligible "
      f"(background has pairs): {sum(s[-1] == 'TRUE' for s in enough)}; background with no pairs: "
      f"{sum(s[7] == 'TRUE' and s[9] == 0 for s in enough)}; background not finished: "
      f"{sum(s[7] == 'FALSE' for s in enough)}")
if any(s[7] == "FALSE" for s in enough):
    print("  (re-run once the background job finishes to pick up the rest)", file=sys.stderr)
