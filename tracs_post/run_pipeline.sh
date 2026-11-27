#!/bin/bash
# Start a whole run and leave it: metadata + config here on the login node (prepare_run.py
# needs the internet for ENA), then pipeline_step.sbatch does within-individual + background
# distances -> merge -> per_reference -> mxsure -> rates.tsv, resubmitting itself between steps.
#
# From the QC_pipeline directory:
#   bash tracs_post/run_pipeline.sh results/tracs_distance/run-$(date +%F) \
#        --max-frac-n 0.5 --min-called 0.5 --snp-threshold 1000000 [other prepare_run.py options]
# Options after the run folder go to prepare_run.py. Environment (all optional, see
# pipeline_step.sbatch): REUSE_FROM=<earlier run> copies that run's finished individuals and
# unchanged background references first (only valid if it used the same snp_threshold);
# MAX_PAIRS, MAX_ROUNDS, WITHIN_TIME, BG_TIME, THROTTLE, RIGHT_TRUNCATION, BOOTSTRAPS, ...
# Watch: tail <run>/pipeline.log    (ends with "PIPELINE FINISHED" or "PIPELINE STOPPED")
set -euo pipefail
[ $# -ge 1 ] || { sed -n 2,13p "$0"; exit 1; }
RUN=$1; shift
[ -f tracs_post/pipeline_step.sbatch ] || { echo "run this from the QC_pipeline directory"; exit 1; }
mkdir -p "${RUN}/logs"
RUN=$(realpath "${RUN}")

if [ ! -s "${RUN}/cmd_sample_times.tsv" ]; then
    module load r/4.4.1
    Rscript tracs_post/export_cmd_times.R "" "${RUN}/cmd_sample_times.tsv"
fi
python3 tracs_post/prepare_run.py --run-dir "${RUN}" "$@"

if [ -n "${REUSE_FROM:-}" ]; then
    old=$(awk -F'\t' '$1=="snp_threshold"{print $2}' "${REUSE_FROM}/run_config.tsv")
    new=$(awk -F'\t' '$1=="snp_threshold"{print $2}' "${RUN}/run_config.tsv")
    [ "${old}" == "${new}" ] || { echo "REUSE_FROM used snp_threshold ${old}, this run ${new}: distances differ, not reusing"; exit 1; }
    python3 tracs_post/build_chunks.py --run-dir "${RUN}" --reuse-from "${REUSE_FROM}" > "${RUN}/logs/build_chunks.reuse.txt"
    python3 tracs_post/build_background.py --run-dir "${RUN}" --reuse-from "${REUSE_FROM}" \
            --max-pairs "${MAX_PAIRS:-10000}" > "${RUN}/logs/build_background.reuse.txt"
fi

j=$(sbatch --parsable --export=ALL,RUN_DIR="${RUN}",STAGE=distances,ROUND=1 \
    --output="${RUN}/logs/pipeline_%j.out" tracs_post/pipeline_step.sbatch)
echo "[$(date '+%F %T')] pipeline started: job ${j}" | tee -a "${RUN}/pipeline.log"
echo "watch: tail ${RUN}/pipeline.log"
