# TRACS combine + distance

Pairwise SNP distances from the TRACS align outputs, read straight out of the squashfs images in
`results/tracs_align_squash/`. There are two runs:

- **Within-individual:** every pair of samples from the same person, on every reference they share.
- **Between-person background:** for each reference, up to 141 samples (≤ 10,000 pairs), one per person,
  spread across studies, for setting a "same strain" threshold.

Run everything from the `QC_pipeline/` directory. All outputs go to `results/tracs_distance/`.

## Inputs

| File | Used for |
|---|---|
| `sample_individual.tsv` | run → individual (`Study_YEAR::subject`). Must include every run you want used. |
| `results/metadata/units.batch*.tsv` | unit (= TRACS sample dir name) → runs, written by the pipeline for each `--run_label` |
| `results/tracs_align_squash/*.sqfs` | TRACS align outputs, one image per batch |

## Running it

```bash
cd /home/zkhurram/scratch_directory/QC_pipeline

# 1. within-individual
python3 tracs_post/build_chunks.py                     # prints the array range
sbatch --array=0-N%40 tracs_post/tracs_combine_distance.sbatch

# 2. background (after step 1 finishes, so samples can be filtered on frac_N)
python3 tracs_post/build_background.py --max-pairs 1000000 --seconds-per-chunk 7200   # frac_N <= 0.75 by default
sbatch --array=0-N%40 tracs_post/tracs_background.sbatch

# 3. one table with both
python3 tracs_post/merge_distances.py                  # -> all_distances.tsv.gz

# 4. substitution rates per reference (mxsure), after steps 1 and 2
module load r/4.4.1
Rscript tracs_post/export_cmd_times.R                  # cMD sampling times -> cmd_sample_times.tsv
python3 tracs_post/build_per_reference.py              # -> per_reference/, prints # eligible
n=$(awk -F'\t' '$11=="TRUE"' results/tracs_distance/per_reference/summary.tsv | wc -l)
sbatch --array=0-$(( (n + 4) / 5 - 1 )) tracs_post/mxsure_rates.sbatch
Rscript tracs_post/mxsure_rates.R combine              # -> mxsure/rates.tsv
```

Replace `N` with the number the build script prints. Don't run a build script while that run's array
job is still going: it rewrites the chunk files the tasks read.

**Check progress:** `squeue -u $USER`. Each task logs to `results/tracs_distance/logs/` or
`.../background/logs/`, and failures leave `failed/<id>.log`. Re-running a build script picks up only
what's missing or failed, so just re-run the build and submit again.

## Outputs

`results/tracs_distance/`
- `per_individual/<id>/` – `distances.csv`, `combined_metadata.csv` (frac_N), `units.txt` (= done)
- `background/per_ref/<ref>.{csv,fracn.tsv,units}`
- `all_distances.tsv.gz` – one row per sample pair per reference. `comparison` is
  `within_individual` / `between_same_study` / `between_study`. Also includes SNPs, filtered SNPs,
  sites compared, SNPs per Mb, frac_N of both samples, and species.
- `background_summary.tsv` – per reference: samples, studies, pairs compared vs written
- `per_reference/<ref>.csv` – the within-individual pairs regrouped by reference, named like
  `background/per_ref/<ref>.csv`. tracs distance columns with `date difference` filled in (days),
  plus individual, body sites, frac_N and `longitudinal` (the pairs mxsure uses: same body site,
  different days, frac_N <= 0.75 both, raw SNPs <= 5000). `per_reference/summary.tsv` counts them;
  `eligible` = >= 30 longitudinal pairs and a background with pairs.
- `cmd_sample_times.tsv` – run -> cMD study, subject, body site, days_from_first_collection,
  infant_age (fallback time scale where a study has no days_from_first_collection)
- `mxsure/rates.tsv` – per reference: lambda (SNPs/site/year), lambda_per_genome, k, SNP threshold at
  1 year, Youden threshold, bootstrap CIs (resampling pairs). Unrelated set = background pairs
  between different studies with frac_N <= 0.75. Per-reference fits and bootstraps in `mxsure/per_ref/*.rds`.
  Time randomisation test (10 time shuffles x 100 bootstraps): `timerand_fail` = share of shuffled-data
  bootstrap rates >= the real estimate, `timerand_pass` = no shuffled 95% CI reaches it. Plots in
  `mxsure/plots/<ref>_timerand.pdf`.

## Reading the results

- Both runs use `--filter -D 5000` (set with `SNP_THRESHOLD`), so pairs with **raw** SNP distance > 5000
  are not written. `merge_distances.py --snp-threshold 5000` applies the same cap to every row.
  That also covers individuals finished before `-D` was added.
  `pairs_compared - pairs_written` in `background_summary.tsv` is the number of background pairs dropped.
- Filter on `sites_considered` and `frac_N` first. Pairs with little overlap show small or zero
  distances.
- Keep `between_same_study` apart from `between_study`. Some studies sampled families.

## mxsure setup

mxsure and its dependencies are installed in `/software/projects/pawsey1172/zkhurram/R/4.4.1`
(`R_LIBS_USER`; the sbatch sets it). When installing R packages, take miniconda3 off `PATH` first:
otherwise conda's headers get picked up and `igraph` (needed by phytools) fails to compile.

## Adding a new batch later

1. Run the pipeline with `--runtracsalign --run_label batchNN` as usual. That writes
   `results/metadata/units.batchNN.tsv`.
2. Squash it into `results/tracs_align_squash/tracs_align_batchNN.sqfs` (copy `squash_batch04.sh`).
   Use gzip compression, as the existing images do.
3. Add the new runs to `sample_individual.tsv`. Runs that aren't in it are skipped.
4. Re-run the three steps above. What gets redone:
   - **Within-individual:** only people whose sample set changed. A person with a new
     timepoint is recomputed from scratch, which takes about a minute. Everyone else is skipped.
   - **Background:** by default, any reference whose selection changed is recomputed. Adding samples
     to a reference that's already at 141 reshuffles its selection. To keep finished references as
     they are and only run new or unfinished ones, use
     `build_background.py --keep-done --max-pairs 1000000`.
   - **Merge:** always rebuilt from everything on disk.

The `.sqfs` listings are cached in `results/tracs_distance/.sqfs_listing/` and refreshed when an image
changes. `--include-unsquashed` on either build script reads sample dirs that are still in
`results/tracs_align/`. Only use it once that batch has finished.

## Cost (Setonix, Oct 2026 tests)

- **Within:** about 27 s of task time per sample. About 360 core-hours for ~6k samples, ~2 h wall.
- **Background:** about 16 stream-hours for batches 00–04 (~65 core-hours). Staging the fastas
  dominates, and `--filter` without `-D` would take hours per common reference.
