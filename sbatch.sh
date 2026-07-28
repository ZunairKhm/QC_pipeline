#!/bin/bash --login
# SLURM batch file for pawsey to run pipeline on accesstions.txt

# Partition for the job:
#SBATCH --partition=work

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="metagenomic-qcpipeline-head"

# The project ID which this job should run under:
#SBATCH --account=pawsey1172

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=02:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load gcc-native/14.2
module load pawseyenv/2025.08
module load nextflow/25.04.6
module load  singularity/4.1.0-slurm
unset SBATCH_EXPORT

# The job command(s):
nextflow run main.nf --kingfisher accessions.txt  --resume