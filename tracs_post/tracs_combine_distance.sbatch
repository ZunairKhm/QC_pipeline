#!/bin/bash --login
# TRACS combine + distance per individual, one array task per chunk built by
# build_chunks.py. Each task runs PARALLEL individuals at once (see
# process_individual.sh for why), staging their TRACS align outputs out of the squashfs
# images into node memory, so nothing is unpacked onto scratch.
#
# From the QC_pipeline directory:
#   python3 tracs_post/build_chunks.py --run-dir $RUN   # prints the full sbatch line to use
#   (it exports RUN_DIR and sends logs to $RUN/logs/; RUN_DIR is required)
# Re-running build_chunks.py afterwards only re-chunks individuals that failed or
# whose sample set changed.
#
# KEEP_COMBINED=1 (sbatch --export=ALL,KEEP_COMBINED=1) also keeps each individual's
# combined MSAs as one tar, if distance will need re-running with other settings.

#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --job-name=tracs-combine-distance
#SBATCH --account=pawsey1172
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=28GB
#SBATCH --time=01:30:00
#SBATCH --output=%x-%A_%a.out

if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

module load pawseyenv/2025.08
module load singularity/4.1.0-slurm
module load squashfs/4.6.1

PIPE=$(realpath "${SLURM_SUBMIT_DIR}")
OUTDIR=${RUN_DIR:?set RUN_DIR (sbatch --export=ALL,RUN_DIR=<run folder>)}
CHUNK=$(printf '%s/chunks/chunk_%04d.tsv' "${OUTDIR}" "${SLURM_ARRAY_TASK_ID}")
PARALLEL=${PARALLEL:-4}                                   # individuals at once
THREADS=$(( SLURM_CPUS_PER_TASK / PARALLEL > 0 ? SLURM_CPUS_PER_TASK / PARALLEL : 1 ))

export TRACS_IMG=${PIPE}/singularity/zunairkhm-tracs-align-v3.0.img
# staged fastas + combined MSAs live in node memory (/dev/shm counts against --mem;
# staging is ~300 MB per sample and is freed before distance runs)
export WORK_ROOT=/dev/shm/tracs_${SLURM_JOB_ID}
export KEEP_COMBINED=${KEEP_COMBINED:-0}
# -D: SNP_THRESHOLD if exported, else the run's snp_threshold (run_config.tsv), else 5000
SNP_THRESHOLD=${SNP_THRESHOLD:-$(awk -F'\t' '$1=="snp_threshold"{print $2}' "${OUTDIR}/run_config.tsv" 2>/dev/null)}
export SNP_THRESHOLD=${SNP_THRESHOLD:-5000}
trap 'rm -rf "${WORK_ROOT}"' EXIT

if [ ! -s "${CHUNK}" ]; then
    echo "no chunk file ${CHUNK}"; exit 1
fi
mkdir -p "${OUTDIR}/per_individual" "${OUTDIR}/failed" "${WORK_ROOT}"
echo "chunk $(basename "${CHUNK}"): $(cut -f1 "${CHUNK}" | uniq | wc -l) individuals," \
     "$(wc -l < "${CHUNK}") samples; ${PARALLEL} at a time x ${THREADS} threads; -D ${SNP_THRESHOLD}"

cut -f1 "${CHUNK}" | uniq | xargs -P "${PARALLEL}" -I{} bash -c '
    if bash "$0" "$1" "$2" "$3" "$4" 2> "$3/failed/$2.log"; then
        rm -f "$3/failed/$2.log"
    else
        echo "FAILED: $2 (see $3/failed/$2.log)"
        exit 1
    fi' "${PIPE}/tracs_post/process_individual.sh" "${CHUNK}" {} "${OUTDIR}" "${THREADS}"
rc=$?

n_failed=$(for s in $(cut -f1 "${CHUNK}" | uniq); do [ -e "${OUTDIR}/failed/${s}.log" ] && echo; done | wc -l)
echo "chunk $(basename "${CHUNK}") finished: ${n_failed} individuals failed"
exit $rc
