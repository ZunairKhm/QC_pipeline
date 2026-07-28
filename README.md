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

nextflow run main.nf --kingfisher accessions.txt  --resume
```
## Notes:
- Primarily set up for pawsey currently but can be modified to work for spartan and other HPCs.s
- Publish directories need to be modified in each process (main.nf for deacon and kingfisher, and under modules for fastp, sylph profile, sylph query, and kraken). 
- Kraken2 can be run optionally with the flag --runkraken2
- nohuman can be run instead of deacon with the flag --filter_method="nohuman"
- Unwanted steps can be commented out in the workflow section of main.nf.
- Currently set up to run sylph profile and sylph query with different databases (gtdb 95 and 99) however that can be easily changed.
- Currently set up for SLURM and Singularity (specifically for pawsey) but this can be modified in nextflow.config (look at commented out code). Remember to launch the pipeline with either sbatch or an interactive session (ideally on workflow nodes).
- nf-boost cleanup would've been very useful for deleting intermediate files but it does not consistently work.
- Maybe important to mention that I copied the available sylph-profile nf core module and swapped out "profile" for "query" to make a sylph-query module.
- For every accession specified, the pipeline will duplicate the files downloaded twice so consider not specifying more than ~250GB worth of accessions for each run for every 1TB of space available. 
- Some utility scripts are provided: sbatch for running bracken (could be implemented into the pipeline but it took 5 mins to run across 100+ samples so its not too annoying); and some kraken and sylph utility scripts to merge the reports; and convert kraken to mpa format)
