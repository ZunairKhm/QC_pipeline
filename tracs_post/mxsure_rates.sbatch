#!/bin/bash --login
# mxsure substitution rates for every eligible reference, split over the array tasks
# (eligible = >= 30 longitudinal pairs and a finished background, from
# per_reference/summary.tsv). Bootstraps (CI and time randomisation test) run in
# parallel across the task's CPUs; the largest reference takes roughly an hour.
#
# From the QC_pipeline directory, once the background run has finished:
#   Rscript tracs_post/export_cmd_times.R           # once (module load r/4.4.1)
#   python3 tracs_post/build_per_reference.py       # prints how many are eligible
#   n=$(awk -F'\t' '$11=="TRUE"' results/tracs_distance/per_reference/summary.tsv | wc -l)
#   sbatch --array=0-$(( (n + 4) / 5 - 1 )) tracs_post/mxsure_rates.sbatch   # ~5 refs/task
#   Rscript tracs_post/mxsure_rates.R combine       # -> results/tracs_distance/mxsure/rates.tsv
# The array must start at 0 with no step: tasks pick refs by ID and SLURM_ARRAY_TASK_COUNT.
# Finished references are skipped, so resubmitting only runs what is missing.

#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --job-name=mxsure-rates
#SBATCH --account=pawsey1426
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28GB
#SBATCH --time=1-00:00:00
#SBATCH --output=results/tracs_distance/mxsure/logs/%A_%a.out

if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

module load r/4.4.1
# set unconditionally: the login environment already sets R_LIBS_USER to another path
export R_LIBS_USER=${MXSURE_R_LIB:-/software/projects/pawsey1172/zkhurram/R/4.4.1}

PIPE=$(realpath "${SLURM_SUBMIT_DIR}")
export OUTDIR=${OUTDIR:-${PIPE}/results/tracs_distance}
# eligible refs (summary.tsv is sorted largest first), dealt round-robin over the array
# tasks so each task gets a mix of large and small references
refs=$(awk -F'\t' -v n="${SLURM_ARRAY_TASK_COUNT}" -v i="${SLURM_ARRAY_TASK_ID}" \
       '$11=="TRUE" && (k++ % n) == i {print $1}' "${OUTDIR}/per_reference/summary.tsv")
if [ -z "${refs}" ]; then
    echo "no references for task ${SLURM_ARRAY_TASK_ID}"; exit 0
fi
echo "task ${SLURM_ARRAY_TASK_ID}: $(echo ${refs} | wc -w) references, ${SLURM_CPUS_PER_TASK} workers"
cd "${PIPE}" && Rscript tracs_post/mxsure_rates.R run ${refs}
