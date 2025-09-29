Nextflow pipeline for preprocessing metagenomic reads and running through sylph profile and query. Edit nextflow.config to change work directory, deacon index location, and sylph database locations. 
```bash
module load Nextflow/x.x.x
module load GCCcore/x.x.x
module load Apptainer/x.x.x
nextflow run main.nf --kingfisher accessions.txt
```
Publish directories need to be modified in each process (main.nf for deacon and kingfisher, and under modules for fastp, sylph profile, and sylph query). 
Unwanted steps can be commented out in the workflow section of main.nf.
Currently set up to run sylph profile and sylph query with different databases (gtdb 95 and 99) however that can be easily changed.
Currently set up for SLURM and Apptainer but this can be modified in nextflow.config (look at commented out code). Remember to launch the pipeline with either sbatch or an interactive session.
nf-boost cleanup would've been very useful for deleting intermediate files but it does not consistently work.
Maybe important to mention that I copied the available sylph-profile nf core module and swapped out "profile" for "query" to make a sylph-query module.
