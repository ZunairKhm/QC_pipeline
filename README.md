# Nextflow pipeline for preprocessing metagenomic reads and running through sylph profile and query. 
Edit nextflow.config to change work directory, and database locations.

for spartan
```bash
module load Nextflow/x.x.x
module load GCCcore/x.x.x
module load Apptainer/x.x.x
nextflow run main.nf --kingfisher accessions.txt
```

for pawsey
```bash
module load gcc-native/14.2
module load pawseyenv/2025.08
module load nextflow/25.04.6
module load singularity/4.1.0-slurm

nextflow run main.nf --kingfisher accessions.txt  -resume
```
## Notes:
- Many samples are split across multiple runs (in these accession lists, 10,523 runs = 6,160 samples). You still just pass an accession list; the pipeline runs `kingfisher annotate` first, groups runs by `sample_accession`, and merges each sample's runs into one unit. Single-run units skip the merge.
- A merged unit is named after its first run accession (lowest sorted), so output filenames stay accession-based and the id is stable across `-resume`. `results/metadata/units.tsv` maps each unit to its runs.
- Some samples mix PAIRED and SINGLE runs, which cannot be concatenated into one pair of files. In those samples the single-end runs are dropped and listed in `results/metadata/dropped_runs.tsv`. On these lists that drops 16% of runs but only 1.4% of bytes - though 28 samples lose more than half their data that way.
- Group on `sample_accession`, never the `number_of_runs_for_sample` metadata field - the latter is unreliable (reports 1 for samples that clearly have 8 runs).
- If some runs of a multi-run unit fail to download, the unit still merges from what arrived. `results/merge_logs/*.merge.tsv` records expected vs merged run counts; grep for PARTIAL to find incomplete units.
- `results/metadata/` holds the raw `run_metadata.tsv`, the derived `units.tsv` (unit -> runs) and `dropped_runs.tsv`.
- Nextflow's own options take a SINGLE dash (-resume, -with-report); a double dash makes it a pipeline parameter instead. `--resume` silently does nothing, so every task re-runs - use `-resume`.
- Primarily set up for pawsey currently but can be modified to work for spartan and other HPCs.s
- Publish directories need to be modified in each process (main.nf for deacon and kingfisher, and under modules for fastp, sylph profile, sylph query, and kraken). 
- Kraken2 can be run optionally with the flag --runkraken2
- TRACS align (strain-level transmission inference, https://github.com/gtonkinhill/tracs, docs: https://gthlab.au/tracs/#/alignment) can be run optionally with the flag --runtracsalign, aligning to the reference fasta file(s) given by --tracs_refseqs (no sourmash database is used, so alignment is direct to the given reference(s) with no on-the-fly genome downloads).
- deacon can be run instead of nohuman with the flag --filter_method="deacon"
- Unwanted steps can be commented out in the workflow section of main.nf.
- Currently set up to run sylph profile and sylph query with different databases (gtdb 95 and 99) however that can be easily changed.
- Currently set up for SLURM and Singularity (specifically for pawsey) but this can be modified in nextflow.config (look at commented out code). Remember to launch the pipeline with either sbatch or an interactive session (ideally on workflow nodes).
- nf-boost cleanup would've been very useful for deleting intermediate files but it does not consistently work.
- Maybe important to mention that I copied the available sylph-profile nf core module and swapped out "profile" for "query" to make a sylph-query module.
- For every accession specified, the pipeline will duplicate the files downloaded twice so consider not specifying more than ~250GB worth of accessions for each run for every 1TB of space available. 
- Each run writes `results/sample_manifest.tsv`: one row per unit with `runs` (how many runs it was built from) and `merge` (COMPLETE / PARTIAL / single_run), then PASS/FAIL per stage. Failed samples are skipped rather than killing the run, so this is how you spot them. Watch for `PARTIAL` - those units passed their stages but were merged from fewer runs than expected, so they are quietly short on data. When running batches into a shared --outdir, pass --run_label <name> so each batch writes its own `sample_manifest.<name>.tsv` instead of overwriting the previous batch's.
- Some utility scripts are provided: sbatch for running bracken (could be implemented into the pipeline but it took 5 mins to run across 100+ samples so its not too annoying); and some kraken and sylph utility scripts to merge the reports; and convert kraken to mpa format)
