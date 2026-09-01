#!/bin/bash
#SBATCH --job-name=rnac_04a
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=120000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Run mlocarna for clusters with 7-19 members that passed step 3.
# Completed clusters are skipped when the job is resubmitted.
#
# Usage:
#   sharcnet/submit.sh 04a DATA_DIR

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR" >&2
    exit 1
fi

data_dir=$1
script_dir=$2

source "$script_dir/info.sh"

module load "$apptainer_module"

cluster_count="$step_02_dir/cluster_count.tsv"
passed_list="$step_03_dir/RNALalifold_passedList.txt"
worker="$script_dir/worker_04.sh"
manifest="$step_04_dir/clusters_04a.tsv"

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

mkdir -p "$step_04_dir"

awk -F '\t' '
    NR == FNR {
        passed[$1] = 1
        next
    }

    $2 in passed && $1 >= 7 && $1 <= 19 {
        print $1 "\t" $2
    }
' "$passed_list" "$cluster_count" |
    LC_ALL=C sort -t $'\t' -k1,1n -k2,2 > "$manifest"

required_cpus=$((locarna_jobs * locarna_threads))

if ((required_cpus > SLURM_CPUS_PER_TASK)); then
    echo \
        "Step 04 requires $required_cpus CPUs, but Slurm allocated $SLURM_CPUS_PER_TASK" \
        >&2
    exit 1
fi

echo "Clusters selected: $(wc -l < "$manifest")"
echo "Concurrent workers: $locarna_jobs"
echo "Threads per worker: $locarna_threads"
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
done < "$manifest" |
    xargs -0 -r \
        -P "$locarna_jobs" \
        -n7 \
        "$worker"

date