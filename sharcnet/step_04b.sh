#!/bin/bash
#SBATCH --job-name=rnac_04b_launch
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Build the manifest of clusters with step04b_cluster_min-step04b_cluster_max
# members (info.sh) that passed step 3, then submit worker_04.sh as a Slurm
# array over it -- one array task per cluster. This job itself is a plain,
# short-lived launcher: it builds the manifest, submits the real array,
# prints where to find it, and exits. It is not where the mlocarna work
# happens.
#
# Submit only after examining step 04a resource usage and timings --
# step04b_cpus/mem/time in info.sh are a scaled-up guess from 04a, not yet
# profiled against real 20-49-member clusters.
#
# Completed clusters are skipped by worker_04.sh, so resubmitting this whole
# step is cheap -- it re-submits an array of the same size, but tasks whose
# cluster already has a result finish immediately.
#
# To retry only specific failed clusters instead of the whole step, see the
# usage comment in worker_04.sh.
#
# Usage:
#   sharcnet/submit.sh 04b DATA_DIR

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR" >&2
    exit 1
fi

data_dir=$1
script_dir=$2

source "$script_dir/info.sh"

cluster_min=$step04b_cluster_min
cluster_max=$step04b_cluster_max

cluster_count="$step_02_dir/cluster_count.tsv"
passed_list="$step_03_dir/RNALalifold_passedList.txt"
worker="$script_dir/worker_04.sh"

[[ -s "$cluster_count" ]] || {
    echo "Missing cluster count: $cluster_count" >&2
    exit 1
}

[[ -s "$passed_list" ]] || {
    echo "Missing or empty passed list: $passed_list" >&2
    exit 1
}

[[ -x "$worker" ]] || {
    echo "Worker is not executable: $worker" >&2
    exit 1
}

run_tag="step_04b_${cluster_min}_${cluster_max}_$(date +%b_%d)"

manifest="$step_04_dir/${run_tag}.tsv"
manifest_tmp="$manifest.tmp"

mkdir -p "$step_04_dir"

awk -F '\t' -v min="$cluster_min" -v max="$cluster_max" '
    NR == FNR {
        passed[$1] = 1
        next
    }

    $2 in passed && $1 >= min && $1 <= max {
        print $1 "\t" $2
    }
' "$passed_list" "$cluster_count" |
    LC_ALL=C sort -t $'\t' -k1,1n -k2,2 \
    > "$manifest_tmp"

mv "$manifest_tmp" "$manifest"

cluster_total=$(wc -l < "$manifest")

if ((cluster_total == 0)); then
    echo \
        "No clusters with $cluster_min-$cluster_max members — nothing to submit." \
        >&2
    exit 1
fi

log_dir="$logs_dir/$run_tag"
mkdir -p "$log_dir"

submission=$(
    sbatch \
        --parsable \
        --job-name=rnac_04b \
        --array="1-${cluster_total}%${step04_jobs}" \
        --cpus-per-task="$step04b_cpus" \
        --mem="$step04b_mem" \
        --time="$step04b_time" \
        --mail-type=END,FAIL \
        --mail-user=bdupin@uwo.ca \
        --output="$log_dir/%A_%a.out" \
        "$worker" \
        "$manifest" \
        "$data_dir" \
        "$script_dir"
)

job_id=${submission%%;*}

printf 'Step 04b: mlocarna for clusters with %s-%s members.\n' \
    "$cluster_min" "$cluster_max"
printf '  Manifest:  %s (%s clusters)\n' "$manifest" "$cluster_total"
printf '  Array job: %s\n' "$job_id"
printf '  Array:     1-%s (%%%s concurrent)\n' "$cluster_total" "$step04_jobs"
printf \
    '  Resources: %s cpus, %s mem, %s time per task\n' \
    "$step04b_cpus" "$step04b_mem" "$step04b_time"
printf '  Logs:      %s/%s_*.out\n' "$log_dir" "$job_id"