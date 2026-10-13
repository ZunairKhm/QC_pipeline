#!/usr/bin/env Rscript
# Sampling times for every run, from curatedMetagenomicData's sampleMetadata.
# The individuals in sample_individual.tsv are cMD study_name::subject_id, so the runs
# join on NCBI_accession (";"-separated runs per cMD sample).
#
#   time_days = days_from_first_collection; infant_age (days of life) is kept as a
#   fallback for studies without it (FerrettiP_2018, OlmMR_2017). build_per_reference.py
#   only compares two samples on the same time scale.
#
# usage (from QC_pipeline/, module load r/4.4.1):
#   Rscript tracs_post/export_cmd_times.R [sampleMetadata.rda] [out.tsv]
# With no .rda given, the current one is downloaded from the cMD GitHub repo.
args <- commandArgs(trailingOnly = TRUE)
rda <- if (length(args) >= 1) args[1] else NA
out <- if (length(args) >= 2) args[2] else "results/tracs_distance/cmd_sample_times.tsv"

if (is.na(rda)) {
  rda <- tempfile(fileext = ".rda")
  download.file("https://github.com/waldronlab/curatedMetagenomicData/raw/devel/data/sampleMetadata.rda",
                rda, mode = "wb", quiet = TRUE)
}
load(rda)  # -> sampleMetadata
m <- sampleMetadata
runs <- strsplit(as.character(m$NCBI_accession), ";", fixed = TRUE)
n <- lengths(runs)
out_df <- data.frame(
  run                        = unlist(runs),
  study                      = rep(m$study_name, n),
  subject                    = rep(m$subject_id, n),
  cmd_sample_id              = rep(m$sample_id, n),
  body_site                  = rep(m$body_site, n),
  days_from_first_collection = rep(m$days_from_first_collection, n),
  infant_age                 = rep(m$infant_age, n),
  visit_number               = rep(m$visit_number, n)
)
out_df <- out_df[!is.na(out_df$run) & out_df$run != "", ]
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.table(out_df, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(nrow(out_df), "runs from", nrow(m), "cMD samples ->", out, "\n")
