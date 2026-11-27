#!/bin/bash --login
# Letter counts of every TRACS consensus (base_composition.py count), one array task per chunk.
#
# From the QC_pipeline directory:
#   python3 tracs_post/base_composition.py build      # prints the full sbatch line (exports COMP_DIR)
#   python3 tracs_post/base_composition.py combine    # after the job: parts -> per-image tables
# The counts describe the images, not a run, so they live in results/tracs_distance/.base_composition/.
# Each task streams ~3000 fastas out of one image with unsquashfs -cat (about a minute).

#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --job-name=tracs-basecomp
#SBATCH --account=pawsey1172
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --output=%x-%A_%a.out

if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

module load pawseyenv/2025.08
module load squashfs/4.6.1

PIPE=$(realpath "${SLURM_SUBMIT_DIR}")
COMP=${COMP_DIR:?set COMP_DIR (sbatch --export=ALL,COMP_DIR=<.base_composition dir>)}
CHUNK=$(printf '%s/chunks/chunk_%04d.tsv' "${COMP}" "${SLURM_ARRAY_TASK_ID}")
[ -s "${CHUNK}" ] || { echo "no chunk file ${CHUNK}"; exit 1; }

cd "${PIPE}"
python3 tracs_post/base_composition.py count --comp-dir "${COMP}" --chunk "${CHUNK}" \
        --processors "${SLURM_CPUS_PER_TASK}"
