#!/bin/bash
# TRACS combine + distance for ONE individual. Called by tracs_combine_distance.sbatch, several at once.
#
# 1. stage: unsquashfs just this individual's posterior_counts fastas + sourmash hits
#    into node memory. Reading through a singularity squashfs mount ran at ~1 MB/s on
#    /scratch; unsquashfs straight off the image is several times faster, and running a
#    few individuals concurrently multiplies that again (latency-bound, not bandwidth).
# 2. tracs combine on the staged dirs, then tracs distance --filter -D SNP_THRESHOLD
#    (default 5000, same as the background run) on every reference
#    seen in >= 2 samples. --filter is single-threaded, hence parallel individuals.
# 3. keep only small outputs in per_individual/<safe_id>/ (units.txt written last =
#    done); with KEEP_COMBINED=1 also a tar of the combined MSAs.
#
# usage: process_individual.sh <chunk.tsv> <safe_id> <outdir> <threads>
set -euo pipefail

chunk=$1; safe=$2; outdir=$3; threads=$4
img=${TRACS_IMG:?TRACS_IMG not set}
out="${outdir}/per_individual/${safe}"
work="${WORK_ROOT:?WORK_ROOT not set}/${safe}"
rm -rf "${work}"; mkdir -p "${work}/stage"
trap 'rm -rf "${work}"' EXIT

awk -F'\t' -v s="${safe}" '$1==s {print $3}' "${chunk}" | sort > "${work}/units.txt"
if [[ -f "${out}/units.txt" && -f "${out}/distances.csv" ]] && cmp -s "${work}/units.txt" "${out}/units.txt"; then
    exit 0
fi
t0=$SECONDS

# --- stage ----------------------------------------------------------------------------
awk -F'\t' -v s="${safe}" '$1==s' "${chunk}" > "${work}/rows.tsv"
for sqfs in $(awk -F'\t' '$4 ~ /\.sqfs$/ {print $4}' "${work}/rows.tsv" | sort -u); do
    awk -F'\t' -v q="${sqfs}" '$4==q {print $3 "/*posterior_counts_ref_*.fasta"; print $3 "/*_sourmash_hits.csv"}' \
        "${work}/rows.tsv" > "${work}/extract.txt"
    unsquashfs -q -n -f -d "${work}/stage" -ef "${work}/extract.txt" -processors "${threads}" "${sqfs}" >/dev/null
done
# staged units read from node memory; units not squashed yet are read in place
awk -F'\t' -v w="${work}/stage" '{print ($4 ~ /\.sqfs$/) ? w "/" $3 : $4}' "${work}/rows.tsv" > "${work}/dirs.txt"
while read -r d; do
    if ! compgen -G "${d}/*posterior_counts_ref_*.fasta" >/dev/null; then
        echo "WARNING ${safe}: no posterior_counts fastas in ${d}" >&2
    fi
done < "${work}/dirs.txt"
t_stage=$((SECONDS - t0))

# --- combine + distance ---------------------------------------------------------------
singularity exec "${img}" tracs combine -i "${work}/dirs.txt" -o "${work}/combined" \
    -t "${threads}" --loglevel WARNING
rm -rf "${work}/stage"   # free memory before distance

: > "${work}/msas.txt"
for f in "${work}"/combined/*_combined.fasta.gz; do
    [[ -e "$f" ]] || continue
    if [[ $(zcat "$f" | grep -c '^>') -ge 2 ]]; then echo "$f" >> "${work}/msas.txt"; fi
done
if [[ -s "${work}/msas.txt" ]]; then
    singularity exec "${img}" tracs distance --msa $(cat "${work}/msas.txt") --filter \
        -D "${SNP_THRESHOLD:-5000}" \
        -t 1 -o "${work}/distances.csv" --loglevel WARNING
else
    echo "sampleA,sampleB,date difference,SNP distance,transmission distance,expected K,filtered SNP distance,sites considered,MSA file" \
        > "${work}/distances.csv"
fi

# --- keep outputs ---------------------------------------------------------------------
mkdir -p "${out}"
rm -f "${out}/units.txt"
cp "${work}/distances.csv" "${work}/combined/combined_metadata.csv" "${out}/"
if [[ "${KEEP_COMBINED:-0}" == 1 ]]; then
    tar -C "${work}/combined" -cf "${out}/combined.tar" .
fi
cp "${work}/units.txt" "${out}/units.txt"

echo "[$(date '+%F %T')] ${safe}: $(wc -l < "${work}/units.txt") samples," \
     "$(($(wc -l < "${out}/distances.csv") - 1)) pairs from $(wc -l < "${work}/msas.txt") refs" \
     "(stage ${t_stage}s, total $((SECONDS - t0))s)"
