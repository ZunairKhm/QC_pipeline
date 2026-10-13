#!/usr/bin/env python3
"""Regroup the within-individual distances by reference, with sampling times, for mxsure.

per_individual/<id>/distances.csv has every within-person pair on every reference; this
writes the same rows split by reference into per_reference/<ref>.csv, named with the
versioned accession like background/per_ref/<ref>.csv, so the two can be paired up per
reference. The columns are tracs distance's own, with "date difference" filled in
(days, from cmd_sample_times.tsv, see export_cmd_times.R), so a file can go straight into
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

Run from QC_pipeline/ after the within-individual run (and export_cmd_times.R). Always
rebuilt from everything on disk, so re-run it once more of the background has finished.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--outdir", default="results/tracs_distance")
p.add_argument("--times", default=None, help="default: <outdir>/cmd_sample_times.tsv")
p.add_argument("--units-glob", default="results/metadata/units.batch*.tsv")
p.add_argument("--max-frac-n", type=float, default=0.75,
               help="longitudinal pairs need frac_N <= this in both samples (as build_background.py)")
p.add_argument("--min-sites", type=int, default=0, help="longitudinal pairs need sites considered >= this")
p.add_argument("--snp-threshold", type=int, default=5000, help="raw SNP distance cap (the runs' -D)")
p.add_argument("--min-comparisons", type=int, default=30,
               help="longitudinal pairs needed for a reference to be eligible")
args = p.parse_args()

unver = lambda acc: acc.split(".")[0]
per_ref_dir = os.path.join(args.outdir, "per_reference")
bg_dir = os.path.join(args.outdir, "background", "per_ref")

# --- unit -> runs -> sampling time ------------------------------------------------------
unit_runs = defaultdict(list)
for f in sorted(glob.glob(args.units_glob)):
    with open(f) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            unit_runs[r["unit_id"]].append(r["run"])

run_time = {}
with open(args.times or os.path.join(args.outdir, "cmd_sample_times.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        run_time[r["run"]] = r

na = lambda v: v in ("", "NA", None)
unit_ind = {}
with open(os.path.join(args.outdir, "index.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        unit_ind[r["unit"]] = r["individual"]


def unit_info(unit):
    """cMD row for a unit (its first run with metadata; a unit is one cMD sample)."""
    for run in unit_runs.get(unit, [unit]):
        if run in run_time:
            return run_time[run]
    return None


def pair_time(a, b):
    """(days apart, which time scale) or (None, reason)."""
    ia, ib = unit_info(a), unit_info(b)
    if ia is None or ib is None:
        return None, "no_cmd_metadata"
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
         "longitudinal"]
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
            sites = int(r["sites considered"])
            longitudinal = (
                t is not None and t > 0
                and ia.get("body_site", "NA") == ib.get("body_site", "NA")
                and not na(fa) and not na(fb)
                and float(fa) <= args.max_frac_n and float(fb) <= args.max_frac_n
                and sites >= args.min_sites
                and not na(r["filtered SNP distance"])
                and int(r["SNP distance"]) <= args.snp_threshold
            )
            ind = unit_ind.get(a, "NA")
            r["date difference"] = "NA" if t is None else f"{t:g}"
            rows[acc].append([r[c] for c in TRACS_COLS] + [
                ind, ind.split("::")[0], ia.get("body_site", "NA"), ib.get("body_site", "NA"),
                src if t is not None else "NA", fa, fb, "TRUE" if longitudinal else "FALSE"])

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
print(f"{len(done)} individuals, {n_pairs} pairs over {len(rows)} references -> {per_ref_dir}/")
print(f"  pairs without a sampling time: {sum(untimed.values())} {dict(untimed)}")
print(f"  references with >= {args.min_comparisons} longitudinal pairs: {len(enough)}; eligible "
      f"(background has pairs): {sum(s[-1] == 'TRUE' for s in enough)}; background with no pairs: "
      f"{sum(s[7] == 'TRUE' and s[9] == 0 for s in enough)}; background not finished: "
      f"{sum(s[7] == 'FALSE' for s in enough)}")
if any(s[7] == "FALSE" for s in enough):
    print("  (re-run once the background job finishes to pick up the rest)", file=sys.stderr)
