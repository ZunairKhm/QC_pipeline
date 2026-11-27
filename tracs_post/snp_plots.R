#!/usr/bin/env Rscript
# SNP-distance plots for chosen references of a finished run: mxsure snp_hist (within-person vs
# unrelated) and snp_over_time with the fitted rate, using the same pair sets as mxsure_rates.R.
#   RUN_DIR=<run> Rscript tracs_post/snp_plots.R <ref> [<ref> ...]
#   RUN_DIR=<run> Rscript tracs_post/snp_plots.R --top 10      # the 10 highest lambdas
# -> <run>/mxsure/plots/<ref>_snp.pdf (2 pages), and with --top (or OUT=<name>.pdf) all in one file there
suppressPackageStartupMessages({library(data.table); library(mxsure); library(ggplot2)})
env <- function(k, d) { v <- Sys.getenv(k); if (v == "") d else v }
outdir <- env("RUN_DIR", ""); if (outdir == "") stop("set RUN_DIR")
cfg <- fread(file.path(outdir, "run_config.tsv"))[, setNames(value, key)]
num <- function(k, d) if (k %in% names(cfg)) as.numeric(cfg[[k]]) else d
max_frac_n <- num("max_frac_n", 1); min_sites <- num("min_sites", 0); min_called <- num("min_called", 0)
right_trunc <- as.numeric(env("RIGHT_TRUNCATION", num("right_truncation", 5000)))
rates <- fread(file.path(outdir, "mxsure", "rates.tsv"))
args <- commandArgs(trailingOnly = TRUE)
top <- length(args) && args[1] == "--top"
refs <- if (top) rates[!is.na(lambda)][order(-lambda)][seq_len(as.integer(args[2])), reference] else args
unit_study <- fread(file.path(outdir, "index.tsv"))[, .(unit, study = sub("::.*", "", individual))]
called <- if (min_called > 0) rbindlist(lapply(Sys.glob(file.path(dirname(normalizePath(outdir)), ".base_composition", "*.tsv")),
                                                fread, select = c("unit", "reference", "frac_called"))) else NULL
plot_dir <- file.path(outdir, "mxsure", "plots"); dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
pages <- list()
for (ref in refs) {
  r <- rates[reference == ref]
  x <- readRDS(file.path(outdir, "mxsure", "per_ref", paste0(ref, ".rds")))
  mixed <- fread(file.path(outdir, "per_reference", paste0(ref, ".csv")))[longitudinal == TRUE]
  bgp <- file.path(outdir, "background", "per_ref", ref)
  bg <- fread(paste0(bgp, ".csv"))
  fn <- fread(paste0(bgp, ".fracn.tsv"))[, .(sample, frac_N)]
  bg <- merge(bg, setnames(copy(fn), c("sampleA", "frac_N_A")), by = "sampleA")
  bg <- merge(bg, setnames(copy(fn), c("sampleB", "frac_N_B")), by = "sampleB")
  bg <- merge(bg, setnames(copy(unit_study), c("sampleA", "studyA")), by = "sampleA")
  bg <- merge(bg, setnames(copy(unit_study), c("sampleB", "studyB")), by = "sampleB")
  un <- bg[studyA != studyB & frac_N_A <= max_frac_n & frac_N_B <= max_frac_n &
             `sites considered` >= min_sites & !is.na(`filtered SNP distance`)]
  if (!is.null(called)) {
    cr <- called[reference == ref, .(unit, frac_called)]
    un <- merge(un, setnames(copy(cr), c("sampleA", "called_A")), by = "sampleA")
    un <- merge(un, setnames(copy(cr), c("sampleB", "called_B")), by = "sampleB")
    un <- un[called_A >= min_called & called_B >= min_called]
  }
  res <- x$result$results
  ttl <- sprintf("%s (%s): %d within-person pairs from %d people, %d unrelated", sub(" strain=.*", "", r$species), ref,
                 nrow(mixed), uniqueN(mixed$individual), nrow(un))
  h <- snp_hist(mixed$`filtered SNP distance`, un$`filtered SNP distance`, snp_threshold = res$snp_threshold,
                limits = c(NA, min(right_trunc, max(c(mixed$`filtered SNP distance`, un$`filtered SNP distance`))) * 1.05),
                title = sprintf("Filtered SNP distance; dashed = mxsure threshold %.0f, est. FP %.2f", res$snp_threshold, res$estimated_fp)) +
    labs(subtitle = ttl)
  o <- snp_over_time(mixed$`filtered SNP distance`, un$`filtered SNP distance`, mixed$`date difference`,
                     mixed$`sites considered`, right_truncation = right_trunc, original_result = res,
                     ci_data = x$ci, threshold_time = 365,
                     title = sprintf("lambda %.2g /site/yr (%.0f SNPs/genome/yr), k %.2f, intercept %.0f, timerand_fail %s",
                                     res$lambda, r$lambda_per_genome, res$k, res$intercept, format(r$timerand_fail))) +
    labs(subtitle = ttl)
  f <- file.path(plot_dir, paste0(ref, "_snp.pdf"))
  pdf(f, width = 11, height = 6); print(h); print(o); dev.off()
  pages <- c(pages, list(h, o)); cat(ref, "->", f, "\n")
}
if (top || env("OUT", "") != "") {
  f <- file.path(plot_dir, env("OUT", "top_rates.pdf")); pdf(f, width = 11, height = 6); for (p in pages) print(p); dev.off(); cat("all ->", f, "\n")
}
