#!/bin/bash --login
# Unattended driver for a whole run: started by run_pipeline.sh, then resubmits itself
# after each array job until the run is finished. Nothing needs doing in between.
#
#   STAGE=distances  build_chunks.py + build_background.py; anything left to do -> submit the
#                    within-individual and background arrays and this job again after them
#                    (afterany). Nothing left (or MAX_ROUNDS used up) -> merge_distances.py,
#                    build_per_reference.py, then STAGE=rates.
#   STAGE=rates      eligible references without an mxsure result -> submit the mxsure array
#                    and this job again after it. All done (or MAX_ROUNDS used up) ->
#                    mxsure_rates.R combine -> mxsure/rates.tsv, and "PIPELINE FINISHED".
# Each stage is retried up to MAX_ROUNDS times, because finished work is skipped and only
# failures (time-outs, node faults, the mxsure port clash) rerun. Progress: <run>/pipeline.log.
#
# env: RUN_DIR (required), STAGE (distances), ROUND (1), MAX_ROUNDS (3),
#      MAX_PAIRS (10000; background pairs per reference), SECONDS_PER_CHUNK (14400),
#      WITHIN_TIME / BG_TIME (12:00:00; array time limits), THROTTLE (60 tasks at once),
#      BG_EXTRA (extra build_background.py options),
#      plus anything the array jobs read (RIGHT_TRUNCATION, BOOTSTRAPS, ...), passed on.

#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --job-name=tracs-pipeline
#SBATCH --account=pawsey1172
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out

if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch (see run_pipeline.sh)"
   exit 1
fi
set -uo pipefail

module load squashfs/4.6.1

PIPE=$(realpath "${SLURM_SUBMIT_DIR}")
RUN=$(realpath "${RUN_DIR:?set RUN_DIR}")
STAGE=${STAGE:-distances}; ROUND=${ROUND:-1}; MAX_ROUNDS=${MAX_ROUNDS:-3}
MAX_PAIRS=${MAX_PAIRS:-10000}; SECONDS_PER_CHUNK=${SECONDS_PER_CHUNK:-14400}
WITHIN_TIME=${WITHIN_TIME:-12:00:00}; BG_TIME=${BG_TIME:-12:00:00}; THROTTLE=${THROTTLE:-60}
LOG=${RUN}/pipeline.log
cd "${PIPE}" || exit 1
mkdir -p "${RUN}/logs" "${RUN}/background/logs" "${RUN}/mxsure/logs"

say()  { echo "[$(date '+%F %T')] $*" | tee -a "${LOG}"; }
die()  { say "PIPELINE STOPPED (${STAGE}, round ${ROUND}): $*"; exit 1; }
nchunks() { find "$1" -maxdepth 1 -name 'chunk_*.tsv' 2>/dev/null | wc -l; }
again() {  # again <after job ids, colon separated> <stage> <round>
    sbatch --parsable --dependency=afterany:"$1" --export=ALL,RUN_DIR="${RUN}",STAGE="$2",ROUND="$3" \
           --output="${RUN}/logs/pipeline_%j.out" tracs_post/pipeline_step.sbatch \
        || die "could not resubmit the pipeline step"
}

if [ "${STAGE}" == "distances" ]; then
    python3 tracs_post/build_chunks.py --run-dir "${RUN}" > "${RUN}/logs/build_chunks.round${ROUND}.txt" \
        || die "build_chunks.py failed (logs/build_chunks.round${ROUND}.txt)"
    python3 tracs_post/build_background.py --run-dir "${RUN}" --max-pairs "${MAX_PAIRS}" \
            --seconds-per-chunk "${SECONDS_PER_CHUNK}" ${BG_EXTRA:-} > "${RUN}/logs/build_background.round${ROUND}.txt" \
        || die "build_background.py failed (logs/build_background.round${ROUND}.txt)"
    nw=$(nchunks "${RUN}/chunks"); nb=$(nchunks "${RUN}/background/chunks")
    if [ $((nw + nb)) -gt 0 ] && [ "${ROUND}" -le "${MAX_ROUNDS}" ]; then
        ids=""
        if [ "${nw}" -gt 0 ]; then
            j=$(sbatch --parsable --array=0-$((nw - 1))%${THROTTLE} --time="${WITHIN_TIME}" \
                --export=ALL,RUN_DIR="${RUN}" --output="${RUN}/logs/%A_%a.out" \
                tracs_post/tracs_combine_distance.sbatch) || die "within-individual sbatch failed"
            ids=${j}; say "distances round ${ROUND}: within-individual array ${j} (${nw} chunks)"
        fi
        if [ "${nb}" -gt 0 ]; then
            j=$(sbatch --parsable --array=0-$((nb - 1))%${THROTTLE} --time="${BG_TIME}" \
                --export=ALL,RUN_DIR="${RUN}" --output="${RUN}/background/logs/%A_%a.out" \
                tracs_post/tracs_background.sbatch) || die "background sbatch failed"
            ids=${ids:+${ids}:}${j}; say "distances round ${ROUND}: background array ${j} (${nb} chunks)"
        fi
        j=$(again "${ids}" distances $((ROUND + 1))); say "next step: job ${j}, after ${ids}"
        exit 0
    fi
    if [ $((nw + nb)) -gt 0 ]; then
        say "WARNING: still unfinished after ${MAX_ROUNDS} rounds (${nw} within chunks, ${nb} background chunks);" \
            "going on without them, see failed/ and background/failed/"
    fi
    say "distances complete; merging"
    python3 tracs_post/merge_distances.py --run-dir "${RUN}" > "${RUN}/logs/merge_distances.txt" \
        || die "merge_distances.py failed (logs/merge_distances.txt)"
    python3 tracs_post/build_per_reference.py --run-dir "${RUN}" > "${RUN}/logs/build_per_reference.txt" \
        || die "build_per_reference.py failed (logs/build_per_reference.txt)"
    STAGE=rates; ROUND=1
fi

if [ "${STAGE}" == "rates" ]; then
    summ=${RUN}/per_reference/summary.tsv
    n=$(awk -F'\t' '$11=="TRUE"' "${summ}" | wc -l)
    n_done=$(awk -F'\t' '$11=="TRUE" {print $1}' "${summ}" \
             | while read -r r; do [ -e "${RUN}/mxsure/per_ref/${r}.tsv" ] && echo; done | wc -l)
    [ "${n}" -gt 0 ] || die "no eligible references in per_reference/summary.tsv"
    if [ "${n_done}" -lt "${n}" ] && [ "${ROUND}" -le "${MAX_ROUNDS}" ]; then
        # the whole array every time: tasks pick references by task count, finished ones are skipped
        j=$(sbatch --parsable --array=0-$(( (n + 4) / 5 - 1 )) --export=ALL,RUN_DIR="${RUN}" \
            --output="${RUN}/mxsure/logs/%A_%a.out" tracs_post/mxsure_rates.sbatch) || die "mxsure sbatch failed"
        say "rates round ${ROUND}: mxsure array ${j} ($((n - n_done)) of ${n} references to do)"
        k=$(again "${j}" rates $((ROUND + 1))); say "next step: job ${k}, after ${j}"
        exit 0
    fi
    [ "${n_done}" -lt "${n}" ] && say "WARNING: $((n - n_done)) references have no mxsure result after ${MAX_ROUNDS} rounds"
    module load r/4.4.1
    export R_LIBS_USER=${MXSURE_R_LIB:-/software/projects/pawsey1172/zkhurram/R/4.4.1}
    RUN_DIR="${RUN}" Rscript tracs_post/mxsure_rates.R combine > "${RUN}/logs/mxsure_combine.txt" 2>&1 \
        || die "mxsure combine failed (logs/mxsure_combine.txt)"
    say "PIPELINE FINISHED: ${n_done} of ${n} eligible references -> ${RUN}/mxsure/rates.tsv"
fi
