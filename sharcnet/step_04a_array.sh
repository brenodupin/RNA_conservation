#!/bin/bash
#SBATCH --job-name=rnac_04a
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=40G
#SBATCH --array=1-10
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Run mlocarna for clusters with 7-19 members that passed step 3, split
# across a Slurm job array. Completed clusters are skipped, so tasks may be
# resubmitted freely.
#
# Each array task selects its share of the manifest by stride rather than by
# contiguous block. The manifest is sorted by member count ascending, so
# contiguous blocks would hand task 1 only cheap clusters and the last task
# only expensive ones. Striding gives every task the same size mix.
#
# The number of array tasks is read from SLURM_ARRAY_TASK_COUNT, so changing
# the --array range above is enough to re-split the work.
#
# Usage:
#   sharcnet/submit.sh 04a_array DATA_DIR

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR" >&2
    exit 1
fi

data_dir=$1
script_dir=$2

source "$script_dir/info.sh"

module load "$apptainer_module"

: "${SLURM_ARRAY_TASK_ID:?This script must be submitted as a job array}"
: "${SLURM_ARRAY_TASK_COUNT:?This script must be submitted as a job array}"

task_id=$SLURM_ARRAY_TASK_ID
task_count=$SLURM_ARRAY_TASK_COUNT

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

# Every task writes only to its own slice file, so no task races another to
# build a shared manifest.

slice_dir="$step_04_dir/array_04a"
mkdir -p "$slice_dir"

slice="$slice_dir/task_${task_id}.tsv"

awk -F '\t' '
    NR == FNR {
        passed[$1] = 1
        next
    }

    $2 in passed && $1 >= 7 && $1 <= 19 {
        print $1 "\t" $2
    }
' "$passed_list" "$cluster_count" |
    LC_ALL=C sort -t $'\t' -k1,1n -k2,2 |
    awk -F '\t' \
        -v task="$task_id" \
        -v count="$task_count" \
        '(FNR - 1) % count == (task - 1) % count' \
        > "$slice"

required_cpus=$((locarna_jobs * locarna_threads))

if ((required_cpus > SLURM_CPUS_PER_TASK)); then
    echo \
        "Step 04 requires $required_cpus CPUs, but Slurm allocated $SLURM_CPUS_PER_TASK" \
        >&2
    exit 1
fi

echo "Array task:          $task_id of $task_count"
echo "Clusters in slice:   $(wc -l < "$slice")"
echo "Concurrent workers:  $locarna_jobs"
echo "Threads per worker:  $locarna_threads"
echo "Per-cluster timeout: $locarna_timeout"

while IFS=$'\t' read -r count cluster; do
    source_fasta="$step_02_dir/splits/${cluster}_cluster.fasta"
    cluster_dir="$step_04_dir/$cluster"
    result="$cluster_dir/${cluster}_result.stk"

    # The final Stockholm file is the completion marker.
    [[ -s "$result" ]] && continue

    printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
        "$rnatools" \
        "$source_fasta" \
        "$cluster_dir" \
        "$cluster" \
        "$locarna_threads" \
        "$fold_temperature" \
        "$locarna_timeout"
done < "$slice" |
    xargs -0 -r \
        -P "$locarna_jobs" \
        -n7 \
        "$worker"

date