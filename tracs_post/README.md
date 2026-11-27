# TRACS combine + distance

Pairwise SNP distances and per-reference substitution rates (mxsure) from the TRACS align outputs.
The outputs are read straight out of the squashfs images under `results/tracs_align_squash/`.

There are two distance runs:

- **Within-individual:** every pair of samples from the same person, on every reference they share.
- **Between-person background:** for each reference, up to K samples, one per person, spread across
  studies, with K*(K-1)/2 ≤ `--max-pairs`. The runs use 1,000,000 pairs (K = 1414); the script default is
  10,000 (K = 141). This is the "unrelated" set for mxsure and for a "same strain" threshold.

Run everything from the `QC_pipeline/` directory.

## How it's organised

**Data sources: `tracs_post/sources.tsv`.** There is one line per source:
`name`, `units_glob`, `squash`, `metadata`.

| metadata | individuals and sampling times come from |
|---|---|
| `cmd` | curatedMetagenomicData `sampleMetadata` (`export_cmd_times.R`). individual = `study_name::subject_id`, time = `days_from_first_collection` (`infant_age` as a fallback) |
| `metalog` | metalog `human_extended_wide_*.tsv.gz`. Each run's ENA sample accession is matched to `spire_sample_name`. individual = `study_code::subject` (study prefix stripped), time = `timepoint` (days since the subject's first sample), body_site = `stool` for human-gut |

Current sources are `cmd` (batches 00–06), `olsson` (Olsson_2022_longitudinal) and `yaffe`
(Yaffe_2025_longitudinal_antibiotic). The helper code is in `tracs_post/sources.py`.

**Run folders: `results/tracs_distance/run-YYYY-MM-DD/`.** Each run is self-contained:

```
run-YYYY-MM-DD/
  run_config.tsv          max_frac_n, snp_threshold, metalog file, created   (prepare_run.py)
  sources.tsv             copy of tracs_post/sources.tsv at prepare time    (later steps read THIS copy)
  cmd_sample_times.tsv    export_cmd_times.R output
  ena_run_sample.tsv      ENA run -> sample accession cache (metalog sources)
  sample_individual.tsv   run -> individual, source
  sample_times.tsv        run -> study, subject, body_site, days_from_first_collection, infant_age, source
  index.tsv  chunks/  per_individual/  logs/  failed/           within-individual run
  background/{selection.tsv,chunks/,per_ref/,logs/,failed/}    background run
  all_distances.tsv.gz  background_summary.tsv                 merge
  per_reference/  mxsure/{per_ref/,plots/,logs/,rates.tsv}     mxsure
```

Every step takes `--run-dir` (shell scripts use `RUN_DIR`). The coverage cut `max_frac_n` is set once in
`run_config.tsv`, and the background selection, the longitudinal pairs and mxsure all read it from there.
The sbatch files refuse to run without `RUN_DIR`. The build scripts print the full `sbatch` line, with
`--export=ALL,RUN_DIR=...` and the log path already filled in.
`results/tracs_distance/.sqfs_listing/` caches the image listings for every run. Delete an image's entries
there if its contents change but its mtime doesn't.
`results/tracs_distance/.base_composition/<image>.tsv` holds the letter counts (A C G T N X other, frac_N,
frac_X, frac_called) of every consensus in every image, from `base_composition.py` (`build` -> sbatch ->
`combine`; `table` prints all of it). It is also shared across runs and rebuilt only for new images.

**Reusing work.** Within-individual distances don't depend on frac_N or on other people. So
`build_chunks.py --reuse-from <earlier run>` copies every finished individual whose unit list hasn't
changed, and only new or changed people run.
`build_background.py --reuse-from <earlier run>` does the same for each reference whose selection is
unchanged. It accepts several folders, checked in order (e.g. `--reuse-from results/tracs_distance/frac_n_0.75
$PREV` for a 0.75 run). Adding a study shifts most selections, though, so expect most references to rerun.

## Running a new run unattended (the usual way)

```bash
cd /home/zkhurram/scratch_directory/QC_pipeline
bash tracs_post/run_pipeline.sh results/tracs_distance/run-$(date +%F) \
     --max-frac-n 0.5 --min-called 0.5        # -D and mxsure truncation default to 20000
tail results/tracs_distance/run-<date>/pipeline.log     # ends "PIPELINE FINISHED" or "PIPELINE STOPPED"
```

`run_pipeline.sh` does step 0 below on the login node (ENA needs the internet), then submits
`pipeline_step.sbatch`. That job builds the chunks, submits the within-individual and background arrays
together, and resubmits itself after them (`afterany`). When nothing is left it merges, builds
`per_reference/`, submits the mxsure array, resubmits itself again and finally writes `mxsure/rates.tsv`.
Each stage is retried up to `MAX_ROUNDS` (3) times; finished work is skipped, so only failures rerun.
Options after the run folder go to `prepare_run.py`. Environment knobs (`MAX_PAIRS`, `WITHIN_TIME`,
`BG_TIME`, `THROTTLE`, `REUSE_FROM`, `RIGHT_TRUNCATION`, ...) are listed at the top of the two scripts.

- The background starts alongside the within-individual run because, with `min_called > 0`,
  `build_background.py` selects samples on frac_called from `.base_composition` and needs no frac_N.
  If new images were added, run `base_composition.py build`/`combine` first.
- `--snp-threshold` (default 20000) is the `-D` of both distance runs; both sbatch files read it from
  `run_config.tsv`, and merge/per_reference drop pairs above it. `-D` cuts on **raw** SNPs; the filtered
  distance mxsure uses is ~0.7x of raw for unrelated pairs. `--filter` time grows ~with raw distance squared,
  so the far pairs cost most: from run-2026-11-12's uncapped distances, `-D 20000` keeps 52% of unrelated
  pairs for ~13% of the filter time (`-D 5000`: 5% of pairs, which is what starved the old unrelated sets).
  `tests/raw_vs_filtered_by_ref.tsv` has the per-reference distributions.
- `--right-truncation` (default 20000; `right_truncation` in `run_config.tsv`) is mxsure's cut, on **filtered**
  SNPs. Change it without new distances: `RIGHT_TRUNCATION=<n>` (or edit run_config), delete
  `mxsure/per_ref/`, resubmit the mxsure array.
- Per-reference background cost at `-D 1000000` (run-2026-11-12, round 1): K=141 references took a median
  3.9 h and up to 11.4 h at 4 threads, and 15 of 118 tasks hit 12 h. At `-D 20000` expect ~1/8 of that.
- `MAX_PAIRS` stays 10000 per reference (141 people).
- To stop a run: `scancel` the `tracs-pipeline` job first, then the arrays.

## Running a new run by hand

```bash
cd /home/zkhurram/scratch_directory/QC_pipeline
RUN=results/tracs_distance/run-$(date +%F)
PREV=results/tracs_distance/run-2026-11-06          # most recent finished run, for --reuse-from
mkdir -p $RUN
module load r/4.4.1 squashfs/4.6.1

# 0. metadata + config
Rscript tracs_post/export_cmd_times.R "" $RUN/cmd_sample_times.tsv
python3 tracs_post/prepare_run.py --run-dir $RUN --max-frac-n 0.5 --min-sites 100000

# 1. within-individual
python3 tracs_post/build_chunks.py --run-dir $RUN --reuse-from $PREV    # prints the sbatch line
#    ... submit it, wait ...

# 2. background (after step 1 finishes, so samples can be filtered on frac_N)
python3 tracs_post/build_background.py --run-dir $RUN --reuse-from $PREV \
        --max-pairs 1000000 --seconds-per-chunk 7200                    # prints the sbatch line
#    ... submit it, wait ...

# 3. one table with both (unfiltered; frac_N columns included)
python3 tracs_post/merge_distances.py --run-dir $RUN

# 4. substitution rates per reference (mxsure)
python3 tracs_post/build_per_reference.py --run-dir $RUN                # prints # eligible
n=$(awk -F'\t' '$11=="TRUE"' $RUN/per_reference/summary.tsv | wc -l); mkdir -p $RUN/mxsure/logs
sbatch --array=0-$(( (n + 4) / 5 - 1 )) --export=ALL,RUN_DIR=$(realpath $RUN) \
       --output=$(realpath $RUN)/mxsure/logs/%A_%a.out tracs_post/mxsure_rates.sbatch
R_LIBS_USER=/software/projects/pawsey1172/zkhurram/R/4.4.1 RUN_DIR=$RUN Rscript tracs_post/mxsure_rates.R combine                 # -> $RUN/mxsure/rates.tsv
```

`prepare_run.py` won't change `max_frac_n`, `snp_threshold` or `min_sites` in an existing run folder. A different
threshold gets a new folder, and `--reuse-from` makes the within-individual part free.
Don't run a build script while that step's array job is still going, because it rewrites the chunk files
the tasks read.

**Check progress:** `squeue -u $USER`. Task logs go to `$RUN/logs/` or `$RUN/background/logs/`, and failures
leave `failed/<id>.log`. Re-running a build script picks up only what's missing or failed. Re-run it and
submit again until it prints `to run: 0`.

## Adding a new data source

1. **cMD batch:** run the pipeline with `--runtracsalign --run_label batchNN`, then squash it into
   `results/tracs_align_squash/tracs_align_batchNN.sqfs` (copy `squash_batch04.sh`, gzip compression).
   The `cmd` line in `sources.tsv` already matches it.
2. **Metalog study:** run the pipeline into its own results dir (e.g. `--outdir results/<name>`). That
   gives `results/<name>/metadata/units.batch*.tsv`. Then:
   - Squash its tracs_align into `results/tracs_align_squash/<name>/tracs_align.sqfs`.
   - Add one line to `sources.tsv`:
     `<name>  results/<name>/metadata/units.batch*.tsv  results/tracs_align_squash/<name>/*.sqfs  metalog`.

   Nothing study-specific is needed beyond that. `prepare_run.py` finds the ENA project from the runs, and
   finds the study_code and subjects in human_extended_wide. Use a newer human_extended_wide with `--metalog` if
   the study isn't in the current one.
   `prepare_run.py` stops if any metalog run can't be matched or has no numeric `timepoint`
   (`--allow-missing-times` keeps untimed runs, which still feed the background).
   Make sure the same samples aren't also in cMD under another study name: the same person would then count
   as two people. A run in two sources is an error, but a re-sequenced sample isn't caught.
3. Start a new run folder (above) with `--reuse-from` the last run.

## Outputs (in the run folder)

- `per_individual/<id>/`: `distances.csv`, `combined_metadata.csv` (frac_N), `units.txt` (= done).
- `background/per_ref/<ref>.{csv,fracn.tsv,units}`.
- `all_distances.tsv.gz`: one row per sample pair per reference. `comparison` is
  `within_individual` / `between_same_study` / `between_study`. Also includes SNPs, filtered SNPs,
  sites compared, SNPs per Mb, frac_N of both samples, and species. **Not** filtered on frac_N.
- `background_summary.tsv`: per reference, the samples, studies, and pairs compared vs written.
- `per_reference/<ref>.csv`: the within-individual pairs regrouped by reference. It has the tracs distance
  columns with `date difference` (days), plus individual, body sites, frac_N and `longitudinal`. The
  longitudinal pairs are the ones mxsure uses: same body site, different days, frac_N ≤ max_frac_n in both
  samples, raw SNPs ≤ the run's `snp_threshold`. `per_reference/summary.tsv` counts them, and `eligible` means ≥ 30
  longitudinal pairs and a background with pairs.
- `mxsure/rates.tsv`: per reference, lambda (SNPs/site/year), lambda_per_genome, k, the SNP threshold at
  1 year, the Youden threshold, and bootstrap CIs (resampling pairs). The unrelated set is background pairs
  between different studies, with frac_N ≤ max_frac_n in both samples. Per-reference fits and bootstraps are
  in `mxsure/per_ref/*.rds`.
  Time randomisation test (10 time shuffles x 100 bootstraps): `timerand_fail` is the share of shuffled-data
  bootstrap rates ≥ the real estimate, and `timerand_pass` means no shuffled 95% CI reaches it. Plots are in
  `mxsure/plots/<ref>_timerand.pdf`.

## Reading the results

- frac_N is the fraction of the reference that is N in a sample's TRACS consensus. `max_frac_n = 0.5`
  means at least 50% of the reference is covered. TRACS align itself already drops references below 25%
  coverage, so 0.75 is the loosest cut that makes sense.
- **"X" sites: frac_N is not the same as called coverage.** In the TRACS align consensus, `N` marks
  sites below the minimum coverage (or masked as coverage outliers), and `X` marks covered sites where no
  allele passed the posterior-frequency threshold (reads, but nothing called). Example:
  `results/tracs_distance/examples/sites0_SRR4052040_vs_SRR4052041_GCF_001025175.1.fasta`, where every
  posterior in SRR4052041 is 0. frac_N counts only `N`, so a sample can have frac_N 0.6 and no A/C/G/T at
  all. Its pairs then report `sites considered` = 0 and 0 SNPs. This is rare at frac_N ≤ 0.5 (0.1% of
  longitudinal pairs) but large at ≤ 0.75 (21%; 70% of between-study pairs). Pairs share either ~0 sites or
  > 500k, with little in between, so `min_sites` (run_config, default 0) removes them. It is applied to
  longitudinal pairs (build_per_reference) and unrelated pairs (mxsure). Use ≥ 100,000. The legacy runs and
  run-2026-11-06/09 used 0. `min_called` (run_config, `prepare_run.py --min-called`) filters on the
  per-sample called share instead; it needs the `.base_composition` tables. `base_composition.py` gives the per-sample frac_X / called share to filter
  samples directly (X is bimodal: most samples have ~0, broken ones 30%+).
- pairsnp (tracs distance) compares every position where both samples have something other than N or X.
  IUPAC ambiguity codes (R, Y, ... = several alleles kept) count as sites, and as SNPs when the allele
  sets don't overlap; they add ~1-4% to sites and up to ~5% to raw SNPs.
- Runs up to run-2026-11-11 used `--filter -D 5000` (`snp_threshold` in `run_config.tsv`), so pairs with
  **raw** SNP distance > 5000 are not written there (later runs: no cut, see above). `pairs_compared - pairs_written` in `background_summary.tsv` is the number of background pairs
  dropped.
- Filter on `sites_considered` and `frac_N` first, because pairs with little overlap show small or zero
  distances.
- Keep `between_same_study` apart from `between_study`, because some studies sampled families.
- Yaffe has ~16 closely spaced timepoints per person (daily around a ciprofloxacin course), so its people
  contribute many short-interval longitudinal pairs per reference.

## History

| folder | what |
|---|---|
| `results/tracs_distance/` (top level) | legacy layout, batches 00–06 cMD: `per_individual/`, `index.tsv`, `cmd_sample_times.tsv`, shared by the two folders below |
| `frac_n_0.75/`, `frac_n_0.5/` | legacy background / per_reference / mxsure / merge at frac_N ≤ 0.75 and ≤ 0.5 (Oct 2026) |
| `run-2026-11-06/` | first run-folder run: cMD + olsson + yaffe, frac_N ≤ 0.5. Within-individual reused from the legacy top level |
| `run-2026-11-09/` | same samples, frac_N ≤ 0.75, min_sites 0. Within-individual copied from run-2026-11-06; background reused from `frac_n_0.75/` then `run-2026-11-06/`. **~21% of longitudinal and ~70% of between-study pairs share 0 sites (see "X sites" above); treat its rates as contaminated** |
| `run-2026-11-10/` | same samples, frac_N ≤ 0.5 **and min_sites 100,000**. Everything except per_reference/mxsure copied from run-2026-11-06 |
| `run-2026-11-11/` | **min_called 0.5** only (frac_called ≥ 0.5 in both samples, from `.base_composition`; implies frac_N ≤ 0.5). Background reused from run-2026-11-10 (selected on frac_N, so fewer usable unrelated pairs). Rates ≈ run-2026-11-10 (median ratio 1.00); fixes some top lambdas, not those with estimated_fp > 0.5 |
| `run-2026-11-12/` | Started with no `-D` cut (1000000); after round 1 (3,094 of 3,114 background refs and all within-person done uncapped) switched to **`snp_threshold` 20000 / `right_truncation` 20000**, so the remaining refs ran at `-D 20000` and merge drops raw > 20000: equivalent to a `-D 20000` run. min_called 0.5, background reselected on frac_called with `MAX_PAIRS` 10000 (3,114 refs, 3.27M pairs), CIs by cluster bootstrap over individuals. First run started with `run_pipeline.sh`. Tests in `tests/background_D/`: Lawsonibacter went from 36 unrelated pairs (median 938 filtered SNPs) to 3,968 (median 7,722). Results: see `pipeline.log` / `mxsure/rates.tsv`. |

The legacy folders are read-only history; new runs never write there.

## Scratch purge / backup

/scratch purges files that haven't been accessed in a while. When backing up, keep each run folder whole,
plus `tracs_post/` (in git) and `results/tracs_distance/.sqfs_listing/` (optional, it rebuilds in minutes).
A run folder holds its own `sample_individual.tsv`, sources and config, so it can be restored and
extended without the top-level files. The cMD run→individual map can always be regenerated with
`export_cmd_times.R` + `prepare_run.py`; this was checked to reproduce the original map exactly.

## mxsure setup

mxsure and its dependencies are installed in `/software/projects/pawsey1172/zkhurram/R/4.4.1`
(`R_LIBS_USER`, which the sbatch sets). When installing R packages, take miniconda3 off `PATH` first.
Otherwise conda's headers get picked up and `igraph` (needed by phytools) fails to compile.

## Cost (Setonix, Oct 2026 tests)

- **Within:** about 27 s of task time per sample. That's about 360 core-hours for ~6k samples, ~2 h wall.
  With `--reuse-from`, only new or changed people cost anything.
- **Background:** about 16 stream-hours for batches 00–04 (~65 core-hours). Staging the fastas
  dominates, and `--filter` without `-D` would take hours per common reference.
- **Image listing:** the first `build_*` call after a new image lists it, which takes seconds to minutes
  per image. The listing is cached afterwards.
