#!/usr/bin/env python3
"""Step 0 of a tracs_post run: create the run folder and its sample metadata.

Writes into --run-dir (e.g. results/tracs_distance/run-2026-11-06/):
  sources.tsv            copy of tracs_post/sources.tsv (later steps read this copy)
  run_config.tsv         max_frac_n, snp_threshold, min_sites, min_called, metalog file, date (later steps read this)
  sample_individual.tsv  run -> individual ("<study>::<subject>") and source name
  sample_times.tsv       run -> study, subject, body_site, days_from_first_collection, infant_age
                         (cmd_sample_times.tsv columns, + source)
  ena_run_sample.tsv     cache of the ENA run -> sample accession lookup for metalog sources

cmd sources: individuals and times come from cmd_sample_times.tsv (export_cmd_times.R, run it
first with the run folder as output). Units without a cMD record are reported and skipped, as
before.
metalog sources: each run's ENA sample accession is looked up in metalog's human_extended_wide
(column spire_sample_name). individual = study_code::subject_id (minus the "study_code." prefix),
days_from_first_collection = timepoint (days since the subject's first sample), body_site = stool
(environmental_package human-gut) or the package name otherwise. Any metalog run that cannot be
matched stops the script, so a new source can't drop out silently.

Re-running is safe; it refuses to change max_frac_n/snp_threshold/min_sites of an existing run unless --force.
"""
import argparse
import csv
import datetime
import glob
import gzip
import os
import shutil
import sys
import time
import urllib.request
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import sources as S

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run-dir", required=True)
p.add_argument("--max-frac-n", type=float, required=True,
               help="coverage cut used by background selection, longitudinal pairs and mxsure "
                    "(frac_N <= this in both samples; 0.5 = at least 50%% of the reference covered)")
p.add_argument("--snp-threshold", type=int, default=20000,
               help="tracs distance -D (raw SNP cap); 20000 since run-2026-11-12, earlier runs 5000")
p.add_argument("--right-truncation", type=int, default=20000,
               help="mxsure right truncation, in FILTERED SNPs (pairs at or above it are left out of the fit); "
                    "20000 since run-2026-11-12 (5000 before)")
p.add_argument("--min-called", type=float, default=0,
               help="longitudinal and unrelated pairs need frac_called (A/C/G/T share, from "
                    "base_composition.py) >= this in both samples")
p.add_argument("--min-sites", type=int, default=0,
               help="longitudinal and unrelated pairs need sites considered >= this. frac_N misses "
                    "'X' sites (reads, but no allele called), so pairs can share 0 called sites; "
                    "see README")
p.add_argument("--sources", default="tracs_post/sources.tsv")
p.add_argument("--metalog", default=None,
               help="human_extended_wide table (default: newest human_extended_wide_*.tsv.gz here)")
p.add_argument("--cmd-times", default=None, help="default: <run-dir>/cmd_sample_times.tsv")
p.add_argument("--allow-missing-times", action="store_true",
               help="metalog runs with no timepoint are kept with time NA instead of stopping")
p.add_argument("--force", action="store_true", help="overwrite an existing run_config.tsv with new values")
args = p.parse_args()

run_dir = args.run_dir
os.makedirs(run_dir, exist_ok=True)
srcs = S.read_sources(args.sources)
na = lambda v: v is None or v.strip() in ("", "NA", "nan")

# --- run_config.tsv -------------------------------------------------------------------
metalog_path = args.metalog
if any(s["metadata"] == "metalog" for s in srcs) and metalog_path is None:
    cands = sorted(glob.glob("human_extended_wide_*.tsv.gz"))
    if not cands:
        sys.exit("ERROR: metalog source listed but no human_extended_wide_*.tsv.gz found; pass --metalog")
    metalog_path = cands[-1]
config = {"max_frac_n": f"{args.max_frac_n:g}", "snp_threshold": str(args.snp_threshold),
          "min_sites": str(args.min_sites),
          "right_truncation": str(args.right_truncation),
          "min_called": f"{args.min_called:g}",
          "metalog": os.path.realpath(metalog_path) if metalog_path else "NA",
          "created": datetime.date.today().isoformat()}
cfg_path = os.path.join(run_dir, "run_config.tsv")
if os.path.exists(cfg_path):
    old = S.read_config(run_dir)
    old.setdefault("min_sites", "0")  # runs made before these existed used 0
    old.setdefault("min_called", "0")
    old.setdefault("right_truncation", "5000")
    changed = [k for k in ("max_frac_n", "snp_threshold", "min_sites", "min_called", "right_truncation") if old.get(k) != config[k]]
    if changed and not args.force:
        sys.exit(f"ERROR: {cfg_path} has different {changed} ({[old.get(k) for k in changed]}); "
                 "start a new run folder, or --force")
    config["created"] = old.get("created", config["created"])
with open(cfg_path, "w") as fh:
    fh.write("key\tvalue\n")
    for k, v in config.items():
        fh.write(f"{k}\t{v}\n")
shutil.copyfile(args.sources, os.path.join(run_dir, "sources.tsv"))

TIME_COLS = ["run", "study", "subject", "cmd_sample_id", "body_site", "days_from_first_collection",
             "infant_age", "visit_number", "source"]
ind_rows, time_rows = [], []  # (run, individual, source), dict
seen_run = {}


def add(run, individual, source, trow):
    if run in seen_run:
        if seen_run[run] != individual:
            sys.exit(f"ERROR: run {run} is in sources {seen_run[run]!r} and {individual!r}")
        return
    seen_run[run] = individual
    ind_rows.append((run, individual, source))
    if trow is not None:
        time_rows.append({**trow, "run": run, "source": source})


# --- cmd sources ----------------------------------------------------------------------
cmd_srcs = [s for s in srcs if s["metadata"] == "cmd"]
if cmd_srcs:
    cmd_times = args.cmd_times or os.path.join(run_dir, "cmd_sample_times.tsv")
    if not os.path.exists(cmd_times):
        sys.exit(f"ERROR: {cmd_times} missing. Run first:\n  module load r/4.4.1\n"
                 f"  Rscript tracs_post/export_cmd_times.R \"\" {cmd_times}")
    cmd = {}
    with open(cmd_times) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cmd.setdefault(r["run"], r)
    for s in cmd_srcs:
        runs, _ = S.unit_runs([s])
        all_runs = {r for rs in runs.values() for r in rs}
        hit = sorted(r for r in all_runs if r in cmd)
        for r in hit:
            row = cmd[r]
            add(r, f"{row['study']}::{row['subject']}", s["name"], {c: row.get(c, "NA") for c in TIME_COLS[:-1]})
        print(f"{s['name']} (cmd): {len(runs)} units, {len(all_runs)} runs, {len(hit)} with a cMD record "
              f"({len(all_runs) - len(hit)} without - skipped, as before)")


# --- metalog sources ------------------------------------------------------------------
def ena(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode()
        except Exception as e:  # transient ENA errors are common
            if i == tries - 1:
                raise
            print(f"  ENA retry ({e})", file=sys.stderr)
            time.sleep(10 * (i + 1))


ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&fields=run_accession,sample_accession,study_accession&accession="
met_srcs = [s for s in srcs if s["metadata"] == "metalog"]
if met_srcs:
    ena_cache = os.path.join(run_dir, "ena_run_sample.tsv")
    run_sample = {}
    if os.path.exists(ena_cache):
        with open(ena_cache) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                run_sample[r["run_accession"]] = r["sample_accession"]

    def lookup(runs):
        """fill run_sample for runs: one project-wide ENA query per project found."""
        missing = sorted(r for r in runs if r not in run_sample)
        tried = set()
        while missing:
            r0 = missing[0]
            lines = ena(ENA + r0).splitlines()
            if len(lines) < 2:
                sys.exit(f"ERROR: ENA has no record of run {r0}")
            proj = lines[1].split("\t")[2]
            if proj in tried:
                sys.exit(f"ERROR: run {r0} not in ENA filereport of its project {proj}")
            tried.add(proj)
            print(f"  ENA: fetching runs of {proj}", file=sys.stderr)
            for line in ena(ENA + proj).splitlines()[1:]:
                f = line.split("\t")
                if len(f) >= 2:
                    run_sample[f[0]] = f[1]
            missing = [r for r in missing if r not in run_sample]
        with open(ena_cache, "w") as fh:
            fh.write("run_accession\tsample_accession\n")
            for k in sorted(run_sample):
                fh.write(f"{k}\t{run_sample[k]}\n")

    met = defaultdict(list)
    with gzip.open(metalog_path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            met[r["spire_sample_name"]].append(r)

    for s in met_srcs:
        runs, _ = S.unit_runs([s])
        if not runs:
            sys.exit(f"ERROR: source {s['name']}: no units files match {s['units_glob']}")
        all_runs = {r for rs in runs.values() for r in rs}
        lookup(all_runs)
        bad, no_time, studies = [], [], Counter()
        for r in sorted(all_runs):
            rows = met.get(run_sample[r], [])
            if len({x["study_code"] for x in rows}) != 1:
                bad.append((r, run_sample[r], len(rows)))
                continue
            m = rows[0]
            study = m["study_code"]
            subj = m["subject_id"]
            subj = subj[len(study) + 1:] if subj.startswith(study + ".") else subj
            if na(subj):
                bad.append((r, run_sample[r], "no subject_id"))
                continue
            t = m.get("timepoint", "")
            try:
                t = f"{float(t):g}"
            except ValueError:
                no_time.append(r)
                t = "NA"
            pkg = m.get("environmental_package", "")
            body = "stool" if pkg == "human-gut" else (pkg or "NA")
            studies[study] += 1
            add(r, f"{study}::{subj}", s["name"], {
                "study": study, "subject": subj, "cmd_sample_id": m["sample_alias"], "body_site": body,
                "days_from_first_collection": t, "infant_age": "NA", "visit_number": "NA"})
        if bad:
            print(f"ERROR: source {s['name']}: {len(bad)} runs have no (or an ambiguous) metalog record:",
                  file=sys.stderr)
            for b in bad[:20]:
                print(f"  run {b[0]}  ENA sample {b[1]}  ({b[2]} metalog rows)", file=sys.stderr)
            sys.exit(1)
        if no_time and not args.allow_missing_times:
            sys.exit(f"ERROR: source {s['name']}: {len(no_time)} runs have no numeric timepoint "
                     f"(e.g. {no_time[:5]}); pass --allow-missing-times to keep them untimed")
        n_ind = len({seen_run[r] for r in all_runs})
        print(f"{s['name']} (metalog): {len(runs)} units, {len(all_runs)} runs, all matched; "
              f"{n_ind} individuals; studies {dict(studies)}" + (f"; {len(no_time)} untimed" if no_time else ""))

# --- a unit must not span individuals -------------------------------------------------
runs_all, _ = S.unit_runs(srcs)
for unit, rs in runs_all.items():
    inds = {seen_run[r] for r in rs if r in seen_run}
    if len(inds) > 1:
        sys.exit(f"ERROR: unit {unit} maps to several individuals: {sorted(inds)}")

with open(os.path.join(run_dir, "sample_individual.tsv"), "w") as fh:
    fh.write("run\tindividual\tsource\n")
    for r in sorted(ind_rows):
        fh.write("\t".join(r) + "\n")
with open(os.path.join(run_dir, "sample_times.tsv"), "w") as fh:
    w = csv.DictWriter(fh, TIME_COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for r in sorted(time_rows, key=lambda x: x["run"]):
        w.writerow({c: ("NA" if na(str(r.get(c, ""))) else r[c]) for c in TIME_COLS})
print(f"{len(ind_rows)} runs, {len({r[1] for r in ind_rows})} individuals -> {run_dir}/sample_individual.tsv")
print(f"config: {config}")
