#!/bin/bash --login
# Between-person background distances: one array task per chunk built by
# build_background.py, each running PARALLEL references at a time
# (process_background_ref.sh).
#
# From the QC_pipeline directory, after build_chunks.py (and ideally after the
# within-individual run, so build_background.py can use frac_N):
#   python3 tracs_post/build_background.py          # prints the --array range to use
#   sbatch --array=0-N%40 tracs_post/tracs_background.sbatch
#   python3 tracs_post/merge_distances.py           # merges within + background
#
# Cost: --filter is very slow on far-apart between-person pairs (>2 CPU-hours for one
# 141-sample reference without -D), so SNP_THRESHOLD (default 5000) is passed as -D,
# which drops pairs with a RAW SNP distance above it before filtering (worst reference
# tested: 42s). Staging fastas then dominates, hence several references at once.

#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --job-name=tracs-background
#SBATCH --account=pawsey1172
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28GB
#SBATCH --time=02:00:00
#SBATCH --output=results/tracs_distance/background/logs/%A_%a.out

if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

module load pawseyenv/2025.08
module load singularity/4.1.0-slurm
module load squashfs/4.6.1

PIPE=$(realpath "${SLURM_SUBMIT_DIR}")
OUTDIR=${OUTDIR:-${PIPE}/results/tracs_distance}
BGDIR=${OUTDIR}/background
CHUNK=$(printf '%s/chunks/chunk_%04d.tsv' "${BGDIR}" "${SLURM_ARRAY_TASK_ID}")
PARALLEL=${PARALLEL:-4}                                   # references at once
THREADS=$(( SLURM_CPUS_PER_TASK / PARALLEL > 0 ? SLURM_CPUS_PER_TASK / PARALLEL : 1 ))

export TRACS_IMG=${PIPE}/singularity/zunairkhm-tracs-align-v3.0.img
export WORK_ROOT=/dev/shm/tracs_bg_${SLURM_JOB_ID}
export SNP_THRESHOLD=${SNP_THRESHOLD:-5000}
trap 'rm -rf "${WORK_ROOT}"' EXIT

if [ ! -s "${CHUNK}" ]; then
    echo "no chunk file ${CHUNK}"; exit 1
fi
mkdir -p "${BGDIR}/per_ref" "${BGDIR}/failed" "${WORK_ROOT}"
echo "chunk $(basename "${CHUNK}"): $(cut -f1 "${CHUNK}" | uniq | wc -l) references," \
     "$(wc -l < "${CHUNK}") fastas; ${PARALLEL} at a time x ${THREADS} threads; -D ${SNP_THRESHOLD}"

cut -f1 "${CHUNK}" | uniq | xargs -P "${PARALLEL}" -I{} bash -c '
    if bash "$0" "$1" "$2" "$3" "$4" 2> "$3/failed/$2.log"; then
        rm -f "$3/failed/$2.log"
    else
        echo "FAILED: $2 (see $3/failed/$2.log)"
        exit 1
    fi' "${PIPE}/tracs_post/process_background_ref.sh" "${CHUNK}" {} "${BGDIR}" "${THREADS}"
rc=$?

n_failed=$(for r in $(cut -f1 "${CHUNK}" | uniq); do [ -e "${BGDIR}/failed/${r}.log" ] && echo; done | wc -l)
echo "chunk $(basename "${CHUNK}") finished: ${n_failed} references failed"
exit $rc
