#!/usr/bin/env Rscript
# Substitution rate per reference with mxsure, from the TRACS distances.
#   mixed     per_reference/<ref>.csv, rows with longitudinal == TRUE (same person,
#             same body site, sampled on different days; see build_per_reference.py)
#   unrelated background/per_ref/<ref>.csv, pairs between people in DIFFERENT studies,
#             frac_N <= MAX_FRAC_N for both samples, sites considered >= MIN_SITES
#             (same-study pairs are left out: some studies sampled families)
# Filtered SNP distances throughout, sites considered as mixed_sites, so lambda is in
# SNPs/site/year (lambda_per_genome = lambda x mean sites). CIs bootstrap whole individuals
# (mxsure_ci cluster_id = individual), since one person's pairs share samples and strains;
# CLUSTER_BOOTSTRAP=0 goes back to resampling pairs (runs up to run-2026-11-11 did that).
# The time randomisation test's own shuffled-data bootstraps still resample pairs.
# Time randomisation test (mxsure_timerandtest): the mixed set's times are shuffled
# TIMERAND_PERMUTATIONS times and refitted with TIMERAND_BOOTSTRAPS bootstraps each; if the
# rate is driven by time, the shuffled estimates fall well below the real one. timerand_fail
# = share of shuffled-data bootstrap rates >= the real estimate (above_estimate criterion);
# timerand_pass = no shuffled 95% CI reaches the real estimate. Plots -> mxsure/plots/.
#
# usage (from QC_pipeline/, module load r/4.4.1):
#   RUN_DIR=<run> Rscript tracs_post/mxsure_rates.R run <ref> [<ref> ...]   # -> mxsure/per_ref/<ref>.{tsv,rds}
#   RUN_DIR=<run> Rscript tracs_post/mxsure_rates.R combine                 # -> mxsure/rates.tsv
# A reference with an existing .tsv is skipped (delete it to redo).
# env: RUN_DIR (required; the run folder), BOOTSTRAPS (500), RIGHT_TRUNCATION (default: right_truncation in run_config.tsv, else 5000; filtered SNPs),
#      THRESHOLD_DAYS (365.25), MAX_FRAC_N (default: max_frac_n in <run>/run_config.tsv), MIN_SITES (default: min_sites in run_config.tsv, else 0), WORKERS (SLURM_CPUS_PER_TASK),
#      TIMERAND_PERMUTATIONS (10; 0 skips the test), TIMERAND_BOOTSTRAPS (100), CLUSTER_BOOTSTRAP (1)
suppressPackageStartupMessages({
  library(data.table)
  library(mxsure)
  library(future)
})
env <- function(k, d) { v <- Sys.getenv(k); if (v == "") d else v }
outdir     <- env("RUN_DIR", env("OUTDIR", ""))
if (outdir == "") stop("set RUN_DIR to the run folder (results/tracs_distance/run-YYYY-MM-DD)")
run_cfg    <- fread(file.path(outdir, "run_config.tsv"))[, setNames(value, key)]
bootstraps <- as.integer(env("BOOTSTRAPS", 500))
right_trunc <- as.numeric(env("RIGHT_TRUNCATION", if ("right_truncation" %in% names(run_cfg)) run_cfg[["right_truncation"]] else 5000))
thresh_days <- as.numeric(env("THRESHOLD_DAYS", 365.25))
max_frac_n <- as.numeric(env("MAX_FRAC_N", run_cfg[["max_frac_n"]]))
min_called <- as.numeric(env("MIN_CALLED", if ("min_called" %in% names(run_cfg)) run_cfg[["min_called"]] else 0))
min_sites  <- as.numeric(env("MIN_SITES", if ("min_sites" %in% names(run_cfg)) run_cfg[["min_sites"]] else 0))
workers    <- as.integer(env("WORKERS", env("SLURM_CPUS_PER_TASK", 1)))
cluster_boot <- env("CLUSTER_BOOTSTRAP", "1") != "0"
tr_perms   <- as.integer(env("TIMERAND_PERMUTATIONS", 10))
tr_boots   <- as.integer(env("TIMERAND_BOOTSTRAPS", 100))
mx_dir <- file.path(outdir, "mxsure", "per_ref")
plot_dir <- file.path(outdir, "mxsure", "plots")
dir.create(mx_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else ""

if (mode == "combine") {
  files <- Sys.glob(file.path(mx_dir, "*.tsv"))
  res <- rbindlist(lapply(files, fread), fill = TRUE)
  summ <- fread(file.path(outdir, "per_reference", "summary.tsv"))
  res <- merge(summ[, .(reference, species, max_time_days)], res,
               by = "reference", all.y = TRUE)
  setorder(res, -n_mixed)
  fwrite(res, file.path(outdir, "mxsure", "rates.tsv"), sep = "\t", na = "NA")
  cat(nrow(res), "references ->", file.path(outdir, "mxsure", "rates.tsv"), "\n")
  quit(save = "no")
}
if (mode != "run" || length(args) < 2) stop("usage: mxsure_rates.R run <ref> ... | combine")

ref_species <- fread(file.path(outdir, "per_reference", "summary.tsv"))[, setNames(species, reference)]
unit_study <- fread(file.path(outdir, "index.tsv"))[, .(unit, study = sub("::.*", "", individual))]
called <- NULL
if (min_called > 0) {  # base_composition.py tables, shared next to the run folders
  comp_files <- Sys.glob(file.path(dirname(normalizePath(outdir)), ".base_composition", "*.tsv"))
  if (!length(comp_files)) stop("min_called > 0 but no .base_composition tables")
  called <- rbindlist(lapply(comp_files, fread, select = c("unit", "reference", "frac_called")))
}
plan(multisession, workers = max(1, workers))

for (ref in args[-1]) {
  out_tsv <- file.path(mx_dir, paste0(ref, ".tsv"))
  if (file.exists(out_tsv)) { cat(ref, ": done already\n"); next }
  t0 <- Sys.time()

  mixed <- fread(file.path(outdir, "per_reference", paste0(ref, ".csv")))[longitudinal == TRUE]

  bg_path <- file.path(outdir, "background", "per_ref", ref)
  bg <- fread(paste0(bg_path, ".csv"))
  fn <- fread(paste0(bg_path, ".fracn.tsv"))[, .(sample, frac_N)]
  bg <- merge(bg, setnames(copy(fn), c("sampleA", "frac_N_A")), by = "sampleA")
  bg <- merge(bg, setnames(copy(fn), c("sampleB", "frac_N_B")), by = "sampleB")
  bg <- merge(bg, setnames(copy(unit_study), c("sampleA", "studyA")), by = "sampleA")
  bg <- merge(bg, setnames(copy(unit_study), c("sampleB", "studyB")), by = "sampleB")
  unrelated <- bg[studyA != studyB & frac_N_A <= max_frac_n & frac_N_B <= max_frac_n &
                    `sites considered` >= min_sites & !is.na(`filtered SNP distance`)]
  if (!is.null(called)) {
    cr <- called[reference == ref, .(unit, frac_called)]
    unrelated <- merge(unrelated, setnames(copy(cr), c("sampleA", "called_A")), by = "sampleA")
    unrelated <- merge(unrelated, setnames(copy(cr), c("sampleB", "called_B")), by = "sampleB")
    unrelated <- unrelated[called_A >= min_called & called_B >= min_called]
  }

  result <- tryCatch(
    mxsure_estimate(mixed$`filtered SNP distance`, unrelated$`filtered SNP distance`,
                    mixed$`date difference`, mixed$`sites considered`,
                    right_truncation = right_trunc, threshold_time = thresh_days, youden = TRUE),
    error = function(e) { message(ref, ": mxsure_estimate failed: ", conditionMessage(e)); NULL })
  ci <- NULL
  if (!is.null(result) && !anyNA(result$results$lambda) && bootstraps > 0) {
    ci <- tryCatch(
      mxsure_ci(mixed$`filtered SNP distance`, unrelated$`filtered SNP distance`,
                mixed$`date difference`, mixed$`sites considered`,
                right_truncation = right_trunc, bootstraps = bootstraps,
                cluster_id = if (cluster_boot) mixed$individual else NA,
                quiet = TRUE),
      error = function(e) { message(ref, ": mxsure_ci failed: ", conditionMessage(e)); NULL })
  }

  timerand <- NULL
  if (!is.null(ci) && is.data.frame(ci$raw_results) && tr_perms > 0) {
    timerand <- tryCatch(
      mxsure_timerandtest(mixed$`filtered SNP distance`, unrelated$`filtered SNP distance`,
                          mixed$`date difference`, mixed$`sites considered`,
                          right_truncation = right_trunc, bootstraps = tr_boots,
                          permutations = tr_perms, start_params = "Efficient",
                          original_result = result$results, ci_data = ci, quiet = TRUE,
                          title = paste(ref, sub(" strain=.*", "", ref_species[ref])), failure_criterion = "above_estimate"),
      error = function(e) { message(ref, ": mxsure_timerandtest failed: ", conditionMessage(e)); NULL })
    if (!is.null(timerand)) {
      ggplot2::ggsave(file.path(plot_dir, paste0(ref, "_timerand.pdf")), timerand$plot,
                      width = 6, height = 4)
    }
  }

  row <- data.table(
    reference = ref, n_mixed = nrow(mixed), n_mixed_individuals = uniqueN(mixed$individual),
    n_mixed_below_truncation = sum(mixed$`filtered SNP distance` < right_trunc),
    n_unrelated = nrow(unrelated),
    n_unrelated_below_truncation = sum(unrelated$`filtered SNP distance` < right_trunc),
    median_mixed_sites = median(mixed$`sites considered`),
    median_time_days = median(mixed$`date difference`))
  if (!is.null(result)) row <- cbind(row, as.data.table(result$results))
  if (!is.null(result$youden)) row <- cbind(row, as.data.table(result$youden))
  if (!is.null(ci)) {
    cis <- ci$confidence_intervals
    row[, `:=`(lambda_lowci = cis$lambda[1], lambda_highci = cis$lambda[2],
               k_lowci = cis$k[1], k_highci = cis$k[2],
               boot_samples = if (is.data.frame(ci$raw_results)) nrow(ci$raw_results) else NA_integer_)]
  }
  if (!is.null(timerand)) {
    o <- timerand$outcome
    row[, `:=`(timerand_permutations = o$n_permutations,
               timerand_pass = !o$any_overlapping_est,
               timerand_n_overlapping_est = o$n_overlapping_est,
               timerand_n_overlapping_lowci = o$n_overlapping_lowci,
               timerand_fail = o$failure_perc)]
  }
  row[, runtime_s := round(as.numeric(difftime(Sys.time(), t0, units = "secs")))]

  saveRDS(list(result = result, ci = ci, timerand = timerand, settings = list(
    bootstraps = bootstraps, timerand_permutations = tr_perms, timerand_bootstraps = tr_boots, right_truncation = right_trunc, threshold_days = thresh_days,
    max_frac_n = max_frac_n, min_sites = min_sites, min_called = min_called, cluster_bootstrap = cluster_boot)), file.path(mx_dir, paste0(ref, ".rds")))
  fwrite(row, out_tsv, sep = "\t", na = "NA")  # written last = done
  cat(sprintf("%s: n=%d (%d people) vs %d unrelated, lambda=%s, timerand fail=%s, %ds\n", ref,
              nrow(mixed), uniqueN(mixed$individual), nrow(unrelated),
              if (is.null(result)) "failed" else format(result$results$lambda, digits = 3),
              if (is.null(timerand)) "NA" else format(timerand$outcome$failure_perc, digits = 3),
              row$runtime_s))
}
