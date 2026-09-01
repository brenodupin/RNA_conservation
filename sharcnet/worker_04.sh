#!/bin/bash

# Predict a common RNA structure for one cluster using mlocarna.
#
# Arguments:
#   container
#   source FASTA
#   cluster directory
#   cluster name
#   threads
#   folding temperature
#   timeout

set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo \
        "Usage: $0 CONTAINER SOURCE_FASTA CLUSTER_DIR CLUSTER THREADS TEMPERATURE TIMEOUT" \
        >&2
    exit 1
fi

rnatools=$1
source_fasta=$2
cluster_dir=$3
cluster=$4
threads=$5
fold_temperature=$6
locarna_timeout=$7

cluster_fasta="$cluster_dir/${cluster}_cluster.fasta"
safe_fasta="$cluster_dir/${cluster}_locarna.fasta"
id_map="$cluster_dir/${cluster}_id_map.tsv"
target_dir="$cluster_dir/${cluster}_locarnap"
result="$cluster_dir/${cluster}_result.stk"
result_tmp="$result.tmp"
timing="$cluster_dir/${cluster}_timing.tsv"
log="$cluster_dir/${cluster}_mlocarna.log"

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

[[ -s "$source_fasta" ]] || {
    echo "Cluster FASTA not found: $source_fasta" >&2
    exit 1
}

mkdir -p "$cluster_dir"

rm -rf "$target_dir"
rm -f "$result" "$result_tmp" "$timing"

trap 'rm -f "$result_tmp"' EXIT

cp "$source_fasta" "$cluster_fasta"

# LocARNA 1.9.2 shortens identifiers in some intermediate files. Long
# accession-and-coordinate identifiers may therefore collide. Give LocARNA
# short unique names and retain a mapping for restoring the final output.
awk -v map="$id_map" '
    /^>/ {
        original = substr($0, 2)
        sub(/[[:space:]].*/, "", original)

        safe = sprintf("S%06d", ++number)

        print safe "\t" original > map
        print ">" safe
        next
    }

    {
        gsub(/[[:space:]]/, "")
        print
    }
' "$source_fasta" > "$safe_fasta"

n_sequences=$(grep -c '^>' "$safe_fasta")

start_time=$(date '+%Y-%m-%d %H:%M:%S')
start_epoch=$(date +%s)
exit_code=0

timeout \
    --signal=TERM \
    --kill-after=5m \
    "$locarna_timeout" \
    apptainer exec \
        --bind "$cluster_dir:/work" \
        --pwd /work \
        --env LANG=C \
        --env LC_ALL=C \
        "$rnatools" \
        mlocarna \
        --probabilistic \
        --tgtdir "${cluster}_locarnap" \
        --moreverbose \
        --stockholm \
        --local-progressive \
        --threads="$threads" \
        --rnafold-temperature="$fold_temperature" \
        "${cluster}_locarna.fasta" \
        > "$log" 2>&1 || exit_code=$?

end_time=$(date '+%Y-%m-%d %H:%M:%S')
duration=$(( $(date +%s) - start_epoch ))

printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$cluster" \
    "$n_sequences" \
    "$start_time" \
    "$end_time" \
    "$duration" \
    "$exit_code" \
    > "$timing"

if ((exit_code != 0)); then
    echo \
        "mlocarna failed for $cluster: exit=$exit_code, sequences=$n_sequences" \
        >&2

    # xargs treats 255 specially and immediately abandons its remaining queue.
    # Return a normal failure so other independent clusters are attempted.
    exit 1
fi

stockholm="$target_dir/results/result.stk"

if [[ ! -s "$stockholm" ]]; then
    stockholm=$(
        find "$target_dir" \
            -type f \
            -name result.stk \
            -print \
            -quit \
            2>/dev/null || true
    )
fi

if [[ -z "$stockholm" || ! -s "$stockholm" ]]; then
    echo "mlocarna produced no result.stk for $cluster" >&2
    exit 1
fi

# Restore the original identifiers in sequence rows and per-sequence
# Stockholm annotations. Write through a temporary file so an interrupted
# conversion cannot leave a partial file that looks complete on resubmission.
awk '
    NR == FNR {
        split($0, fields, "\t")
        original[fields[1]] = fields[2]
        next
    }

    $1 in original {
        $1 = original[$1]
        print
        next
    }

    ($1 == "#=GS" || $1 == "#=GR") && $2 in original {
        $2 = original[$2]
        print
        next
    }

    {
        print
    }
' OFS='\t' "$id_map" "$stockholm" > "$result_tmp"

mv "$result_tmp" "$result"