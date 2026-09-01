#!/bin/bash
#SBATCH --job-name=rnac_03b
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Screen all Clustal Omega alignments with RNALalifold using parallel workers.
# Outputs are written to DATA_DIR/03_screen/<cluster>.
#
# Usage:
#   sharcnet/submit.sh 03b DATA_DIR

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR" >&2
    exit 1
fi

data_dir=$1
script_dir=$2

echo " step_03b.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

module load "$apptainer_module"

cluster_count="$step_02_dir/cluster_count.tsv"
worker="$script_dir/worker_03b.sh"

while IFS=$'\t' read -r count cluster; do
    cluster_dir="$step_03_dir/$cluster"
    alignment="${cluster}_aligned.aln"
    output="$cluster_dir/${cluster}_RNALalifold.out"

    printf '%s\0%s\0%s\0%s\0%s\0' \
        "$rnatools" \
        "$cluster_dir" \
        "$alignment" \
        "$output" \
        "$fold_temperature"
done < "$cluster_count" |
    xargs -0 -r \
        -P "$SLURM_CPUS_PER_TASK" \
        -n5 \
        "$worker"

date