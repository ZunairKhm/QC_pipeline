#!/usr/bin/env python3
"""Merge within-individual and between-person background distances into one table.

One row per (sample pair, reference), with `comparison` set to
  within_individual    same person (per_individual/ run, every pair written)
  between_same_study   different people, same study (some studies sampled families)
  between_study        different people, different studies
and each sample's frac_N for that reference, the species name, and filtered SNPs per Mb
of sites compared. Low coverage makes SNP distances look small, so filter on
sites_considered and frac_N before calling anything the same strain.

Both runs use tracs distance -D 5000, which drops pairs whose RAW snp_distance is above
5000. Individuals finished before -D was added to the within run have uncapped rows,
so the same cap (--snp-threshold) is applied here to every row. background_summary.tsv
gives, per reference, how many background pairs were compared in total, so the dropped
pairs can be counted.

Note: tracs distance writes the reference as basename.split(".")[0], which drops the
accession version (GCF_000012825.1 -> GCF_000012825), so joins use the unversioned
accession.
"""
import argparse
import csv
import glob
import gzip
import os
import sys
from collections import defaultdict

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--outdir", default="results/tracs_distance")
p.add_argument("--snp-threshold", type=int, default=5000,
               help="drop pairs with raw SNP distance above this (match the runs' -D)")
p.add_argument("-o", "--output", default=None, help="default: <outdir>/all_distances.tsv.gz")
args = p.parse_args()

out_path = args.output or os.path.join(args.outdir, "all_distances.tsv.gz")
unver = lambda acc: acc.split(".")[0]
study_of = lambda ind: ind.split("::")[0]

unit_ind, safe_to_ind = {}, {}
with open(os.path.join(args.outdir, "index.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        unit_ind[r["unit"]] = r["individual"]
        safe_to_ind[r["safe_id"]] = r["individual"]

# frac_N per (sample, ref) and species per ref, from every combined_metadata.csv.
# combined_metadata.csv has 10 fixed columns then the species name, which may itself
# contain commas (tracs writes it unquoted).
frac_n, species = {}, {}
done_ind = []
for d in sorted(glob.glob(os.path.join(args.outdir, "per_individual", "*"))):
    if not os.path.exists(os.path.join(d, "units.txt")):
        continue  # unfinished
    done_ind.append(d)
    with open(os.path.join(d, "combined_metadata.csv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split(",")
            acc = unver(f[1])
            frac_n[(f[0], acc)] = f[9]
            species.setdefault(acc, ",".join(f[10:]))

bg_dir = os.path.join(args.outdir, "background", "per_ref")
done_ref = [f[: -len(".units")] for f in sorted(glob.glob(os.path.join(bg_dir, "*.units")))]
for ref in done_ref:  # the background measures frac_N itself
    with open(ref + ".fracn.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            frac_n.setdefault((r["sample"], unver(os.path.basename(ref))), r["frac_N"])

cols = ["comparison", "reference", "species", "sampleA", "sampleB", "individualA", "individualB",
        "studyA", "studyB", "snp_distance", "filtered_snp_distance", "sites_considered",
        "filtered_snps_per_Mb", "frac_N_A", "frac_N_B"]
counts = defaultdict(int)


def write_rows(out, csv_path, comparison=None):
    with open(csv_path) as fh:
        for r in csv.DictReader(fh):
            if int(r["SNP distance"]) > args.snp_threshold:
                counts["dropped_over_threshold"] += 1
                continue
            a, b, ref = r["sampleA"], r["sampleB"], r["MSA file"]
            ia, ib = unit_ind.get(a, "NA"), unit_ind.get(b, "NA")
            comp = comparison or ("between_same_study" if study_of(ia) == study_of(ib) else "between_study")
            sites = int(r["sites considered"])
            filt = r["filtered SNP distance"]
            per_mb = (f"{int(filt) / sites * 1e6:.3f}"
                      if sites > 0 and filt not in ("", "NA") else "NA")
            out.write("\t".join([
                comp, ref, species.get(ref, "NA"), a, b, ia, ib, study_of(ia), study_of(ib),
                r["SNP distance"], filt, str(sites), per_mb,
                frac_n.get((a, ref), "NA"), frac_n.get((b, ref), "NA"),
            ]) + "\n")
            counts[comp] += 1


with gzip.open(out_path, "wt") as out:
    out.write("\t".join(cols) + "\n")
    for d in done_ind:
        write_rows(out, os.path.join(d, "distances.csv"), "within_individual")
    for ref in done_ref:
        write_rows(out, ref + ".csv")

if done_ref:
    summ = os.path.join(args.outdir, "background_summary.tsv")
    with open(summ, "w") as fh:
        fh.write("reference\tspecies\tn_samples\tn_studies\tpairs_compared\tpairs_written\n")
        for ref in done_ref:
            units = open(ref + ".units").read().split()
            n = sum(1 for _ in open(ref + ".fracn.tsv")) - 1  # samples actually in the MSA
            written = sum(1 for _ in open(ref + ".csv")) - 1
            acc = unver(os.path.basename(ref))
            n_st = len({study_of(unit_ind.get(u, "NA")) for u in units})
            fh.write(f"{os.path.basename(ref)}\t{species.get(acc, 'NA')}\t{n}\t{n_st}\t"
                     f"{n * (n - 1) // 2}\t{written}\n")

print(f"{len(done_ind)} individuals, {len(done_ref)} background references -> {out_path}")
for k in ("within_individual", "between_same_study", "between_study"):
    print(f"  {k}: {counts[k]} pairs")
if counts["dropped_over_threshold"]:
    print(f"  ({counts['dropped_over_threshold']} pairs above {args.snp_threshold} raw SNPs dropped)")
for sub, what in (("failed", "individuals"), (os.path.join("background", "failed"), "background references")):
    n_failed = len(glob.glob(os.path.join(args.outdir, sub, "*.log")))
    if n_failed:
        print(f"WARNING: {n_failed} {what} failed - see {os.path.join(args.outdir, sub)}/", file=sys.stderr)
