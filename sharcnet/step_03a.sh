#!/bin/bash
#SBATCH --job-name=rnac_03a
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Align all clusters with Clustal Omega using parallel workers.
# Outputs are written to DATA_DIR/03_screen/<cluster>.
#
# Usage:
#   sharcnet/submit.sh 03a DATA_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 DATA_DIR" >&2
    exit 1
fi

data_dir=$1

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/info.sh"

module load "$clustalo_module"

cluster_count="$step_02_dir/cluster_count.tsv"
worker="$script_dir/worker_03a.sh"

mkdir -p "$step_03_dir"

while IFS=$'\t' read -r count cluster; do
    source_fasta="$step_02_dir/splits/${cluster}_cluster.fasta"
    cluster_dir="$step_03_dir/$cluster"

    cluster_fasta="$cluster_dir/${cluster}_cluster.fasta"
    alignment="$cluster_dir/${cluster}_aligned.aln"
    distance_matrix="$cluster_dir/${cluster}_distMat.csv"

    printf '%s\0%s\0%s\0%s\0' \
        "$source_fasta" \
        "$cluster_fasta" \
        "$alignment" \
        "$distance_matrix"
done < "$cluster_count" |
    xargs -0 -r \
        -P "$SLURM_CPUS_PER_TASK" \
        -n4 \
        "$worker"

date