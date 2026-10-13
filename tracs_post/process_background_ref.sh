#!/bin/bash
# Between-person background for ONE reference. Called by tracs_background.sbatch.
#   1. unsquashfs just this reference's posterior_counts fasta for each selected sample
#   2. build_msa.py -> one MSA (+ per-sample frac_N)
#   3. tracs distance --filter -D SNP_THRESHOLD (default 5000): pairs further apart than
#      that are not written, which is what makes --filter affordable on between-person pairs
# Outputs in background/per_ref/: <ref>.csv, <ref>.fracn.tsv, <ref>.units (last = done).
#
# usage: process_background_ref.sh <chunk.tsv> <ref> <bgdir> <threads>
set -euo pipefail

chunk=$1; ref=$2; bgdir=$3; threads=$4
img=${TRACS_IMG:?TRACS_IMG not set}
here=$(dirname "$(realpath "$0")")
out="${bgdir}/per_ref"
work="${WORK_ROOT:?WORK_ROOT not set}/${ref}"
rm -rf "${work}"; mkdir -p "${work}/stage" "${out}"
trap 'rm -rf "${work}"' EXIT

awk -F'\t' -v r="${ref}" '$1==r' "${chunk}" > "${work}/rows.tsv"
cut -f2 "${work}/rows.tsv" | sort > "${work}/units.txt"
if [[ -f "${out}/${ref}.units" && -f "${out}/${ref}.csv" ]] && cmp -s "${work}/units.txt" "${out}/${ref}.units"; then
    exit 0
fi
t0=$SECONDS

for sqfs in $(awk -F'\t' '$5 ~ /\.sqfs$/ {print $5}' "${work}/rows.tsv" | sort -u); do
    awk -F'\t' -v q="${sqfs}" -v r="${ref}" '$5==q {print $2 "/" $2 "_posterior_counts_ref_" r ".fasta"}' \
        "${work}/rows.tsv" > "${work}/extract.txt"
    unsquashfs -q -n -f -d "${work}/stage" -ef "${work}/extract.txt" -processors 2 "${sqfs}" >/dev/null
done
t_stage=$((SECONDS - t0))

mapfile -t pairs < <(awk -F'\t' -v w="${work}/stage" -v r="${ref}" '{
    d = ($5 ~ /\.sqfs$/) ? w "/" $2 : $5
    print $2 "=" d "/" $2 "_posterior_counts_ref_" r ".fasta"}' "${work}/rows.tsv")
n_seq=$(python3 "${here}/build_msa.py" "${work}/${ref}.fasta.gz" "${work}/${ref}.fracn.tsv" "${pairs[@]}")
rm -rf "${work}/stage"

if [[ ${n_seq} -ge 2 ]]; then
    singularity exec "${img}" tracs distance --msa "${work}/${ref}.fasta.gz" --filter \
        -D "${SNP_THRESHOLD:-5000}" -t "${threads}" -o "${work}/${ref}.csv" --loglevel WARNING
else
    echo "sampleA,sampleB,date difference,SNP distance,transmission distance,expected K,filtered SNP distance,sites considered,MSA file" \
        > "${work}/${ref}.csv"
fi

rm -f "${out}/${ref}.units"
cp "${work}/${ref}.csv" "${work}/${ref}.fracn.tsv" "${out}/"
cp "${work}/units.txt" "${out}/${ref}.units"
echo "[$(date '+%F %T')] ${ref}: ${n_seq} samples, $(($(wc -l < "${out}/${ref}.csv") - 1)) pairs" \
     "(stage ${t_stage}s, total $((SECONDS - t0))s)"
