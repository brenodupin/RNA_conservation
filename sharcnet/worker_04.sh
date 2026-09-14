#!/bin/bash
#SBATCH --job-name=rnac_04
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --array=1

# Predict a common RNA structure for one cluster using mlocarna. Runs as one
# task of a Slurm array. step_04a.sh and step_04b.sh are the normal way to
# submit this: they build a manifest for their own size range, then submit
# this script as an array over it with the resource requests appropriate for
# that range (--cpus-per-task/--mem/--time on their own sbatch call — the
# #SBATCH lines above are just placeholders for a direct invocation).
#
# The manifest is one line per cluster: "<member count>\t<cluster name>".
# Array task N processes line N. mlocarna's thread count is read from
# SLURM_CPUS_PER_TASK, so it always matches whatever --cpus-per-task the
# submitter actually requested.
#
# Completed clusters (an existing, non-empty *_result.stk) are skipped, so a
# task that already finished is a cheap no-op on resubmission, and a
# specific failed array index can be resubmitted directly without going
# through step_04a.sh/step_04b.sh or rebuilding the manifest, e.g.:
#
#   sbatch --array=3,17,412 \
#       --cpus-per-task=3 --mem=6G --time=1-01:00:00 \
#       --output=DATA_DIR/logs/step_04a_7_19_Sep_14/%A_%a.out \
#       SHARCNET_DIR/worker_04.sh MANIFEST DATA_DIR SHARCNET_DIR
#
# (match --cpus-per-task/--mem/--time to whatever the original run used, and
# reuse the --output folder step_04a.sh/step_04b.sh printed so the retried
# tasks land next to the rest of that run's logs. Find failed indices with:
#   sacct -j JOBID -X --format=JobID,State -P |
#       awk -F'|' '$2=="FAILED"{print $1}' | sed 's/.*_//' | paste -sd,)
#
# Pass --force as a 4th argument to redo a cluster even if it already has a
# result (step_04a.sh/step_04b.sh forward their own --force here too).
#
# Usage:
#   sbatch [--array=... --cpus-per-task=N --mem=... --time=...] \
#       worker_04.sh MANIFEST DATA_DIR SHARCNET_DIR [--force]

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: sbatch $0 MANIFEST DATA_DIR SHARCNET_DIR [--force]" >&2
    exit 1
fi

manifest=$1
data_dir=$2
script_dir=$3
force=${4:-}

source "$script_dir/info.sh"

module load "$apptainer_module"

: "${SLURM_ARRAY_TASK_ID:?This script must be submitted as a job array}"

[[ -s "$manifest" ]] || {
    echo "Missing manifest: $manifest" >&2
    exit 1
}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")

[[ -n "$line" ]] || {
    echo \
        "No manifest entry for array task $SLURM_ARRAY_TASK_ID" \
        "(manifest has $(wc -l < "$manifest") lines)" \
        >&2
    exit 1
}

IFS=$'\t' read -r count cluster <<< "$line"

threads=$SLURM_CPUS_PER_TASK

source_fasta="$step_02_dir/splits/${cluster}_cluster.fasta"
cluster_dir="$step_04_dir/$cluster"
cluster_fasta="$cluster_dir/${cluster}_cluster.fasta"
safe_fasta="$cluster_dir/${cluster}_locarna.fasta"
id_map="$cluster_dir/${cluster}_id_map.tsv"
target_dir="$cluster_dir/${cluster}_locarnap"
result="$cluster_dir/${cluster}_result.stk"
result_tmp="$result.tmp"
timing="$cluster_dir/${cluster}_timing.tsv"
log="$cluster_dir/${cluster}_mlocarna.log"

echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Cluster:    $cluster ($count members)"
echo "Threads:    $threads"

# The final Stockholm file is the completion marker, unless --force says to
# redo it anyway.
if [[ -s "$result" && "$force" != "--force" ]]; then
    echo "Already completed, skipping."
    exit 0
fi

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

echo "Completed cluster $cluster: $n_sequences sequences, result=$result"
rm -rf "$target_dir" "$cluster_fasta" "$safe_fasta" "$id_map"

date
