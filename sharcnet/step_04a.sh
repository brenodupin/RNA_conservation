#!/bin/bash
#SBATCH --job-name=rnac_04a_launch
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Build the manifest of clusters with step04a_cluster_min-step04a_cluster_max
# members (info.sh) that passed step 3 and don't already have a result, then
# submit worker_04.sh as a Slurm array over it -- one array task per
# cluster. This job itself is a plain, short-lived launcher: it builds the
# manifest, submits the real array, prints where to find it, and exits. It
# is not where the mlocarna work happens.
#
# The manifest excludes clusters that already have a *_result.stk, so a
# rerun of this script only ever submits the clusters still left to do.
# Combined with the account's AssocMaxSubmitJobLimit (1000 total submitted
# jobs, array elements included) via step04_max_submit in info.sh: if more
# clusters remain than that, this submits only step04_max_submit of them
# and tells you to rerun once that batch has drained -- no need to track
# what's left by hand.
#
# To retry only specific failed clusters instead of the whole step, see the
# usage comment in worker_04.sh.
#
# Pass --force to redo every cluster in range, including ones that already
# have a result -- the manifest then includes them too, and worker_04.sh is
# told not to skip them.
#
# Usage:
#   sharcnet/submit.sh 04a DATA_DIR [--force]

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force]" >&2
    exit 1
fi

data_dir=$1
script_dir=$2
force=${3:-}

if [[ -n "$force" && "$force" != "--force" ]]; then
    echo "Unknown argument: $force (expected --force)" >&2
    exit 1
fi

source "$script_dir/info.sh"

cluster_min=$step04a_cluster_min
cluster_max=$step04a_cluster_max

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

run_tag="step_04a_${cluster_min}_${cluster_max}_$(date +%b_%d)"
[[ "$force" == "--force" ]] && run_tag="${run_tag}_force"

manifest="$step_04_dir/${run_tag}.tsv"
manifest_tmp="$manifest.tmp"

mkdir -p "$step_04_dir"

# Filter by size and passed-screening, then drop anything that already has
# a completed result -- so the array only ever covers remaining work.
# --force skips that last drop, so already-completed clusters are included
# (and re-run) too.
awk -F '\t' -v min="$cluster_min" -v max="$cluster_max" '
    NR == FNR {
        passed[$1] = 1
        next
    }

    $2 in passed && $1 >= min && $1 <= max {
        print $1 "\t" $2
    }
' "$passed_list" "$cluster_count" |
    LC_ALL=C sort -t $'\t' -k1,1n -k2,2 |
    while IFS=$'\t' read -r count cluster; do
        result="$step_04_dir/$cluster/${cluster}_result.stk"
        [[ "$force" != "--force" && -s "$result" ]] && continue
        printf '%s\t%s\n' "$count" "$cluster"
    done \
    > "$manifest_tmp"

mv "$manifest_tmp" "$manifest"

cluster_total=$(wc -l < "$manifest")

if ((cluster_total == 0)); then
    echo \
        "No remaining clusters with $cluster_min-$cluster_max members" \
        "-- nothing to submit." \
        >&2
    exit 0
fi

submit_count=$cluster_total
remaining_after=0

if ((cluster_total > step04_max_submit)); then
    submit_count=$step04_max_submit
    remaining_after=$((cluster_total - step04_max_submit))
fi

log_dir="$logs_dir/$run_tag"
mkdir -p "$log_dir"

worker_args=("$manifest" "$data_dir" "$script_dir")
[[ "$force" == "--force" ]] && worker_args+=(--force)

submission=$(
    sbatch \
        --parsable \
        --job-name=rnac_04a \
        --array="1-${submit_count}%${step04_jobs}" \
        --cpus-per-task="$step04a_cpus" \
        --mem="$step04a_mem" \
        --time="$step04a_time" \
        --mail-type=END,FAIL \
        --mail-user=bdupin@uwo.ca \
        --output="$log_dir/%A_%a.out" \
        "$worker" \
        "${worker_args[@]}"
)

job_id=${submission%%;*}

printf 'Step 04a: mlocarna for clusters with %s-%s members.\n' \
    "$cluster_min" "$cluster_max"
[[ "$force" == "--force" ]] && printf '  Mode:      --force (redoing completed clusters too)\n'
printf '  Manifest:  %s (%s remaining clusters)\n' "$manifest" "$cluster_total"
printf '  Array job: %s\n' "$job_id"
printf '  Array:     1-%s (%%%s concurrent)\n' "$submit_count" "$step04_jobs"
printf \
    '  Resources: %s cpus, %s mem, %s time per task\n' \
    "$step04a_cpus" "$step04a_mem" "$step04a_time"
printf '  Logs:      %s/%s_*.out\n' "$log_dir" "$job_id"

if ((remaining_after > 0)); then
    rerun_cmd="sharcnet/submit.sh 04a $data_dir"
    [[ "$force" == "--force" ]] && rerun_cmd="$rerun_cmd --force"

    printf '\n'
    printf \
        '  NOTE: %s clusters exceed the %s-submission cap (AssocMaxSubmitJobLimit=1000).\n' \
        "$cluster_total" "$step04_max_submit"
    printf \
        '        Submitted %s now; %s left for next time.\n' \
        "$submit_count" "$remaining_after"
    printf \
        '        Once this batch drains, rerun: %s\n' \
        "$rerun_cmd"
fi
