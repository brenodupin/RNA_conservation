#!/bin/bash
#SBATCH --job-name=rnac_02
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Add reverse-complement windows, perform two-pass MMseqs2 clustering, and
# split multi-member clusters into individual FASTA files.
# Outputs are written to DATA_DIR/02_clusters.
#
# Usage:
#   sharcnet/submit.sh 02 DATA_DIR

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR" >&2
    exit 1
fi

data_dir=$1
script_dir=$2

echo " step_02.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

module load "$seqtk_module"
module load "$mmseqs_module"

mkdir -p "$step_02_dir"

unique="$step_01_dir/${windows_prefix}_processed_noNs_polyN_uniq.fasta"
reverse="$step_02_dir/${windows_prefix}_revComp.fasta"
both="$step_02_dir/${windows_prefix}_bothStrands.fasta"

pass_1="$step_02_dir/mmseqs_pass_1"
pass_2="$step_02_dir/mmseqs_pass_2"

tmp_1="$SLURM_TMPDIR/mmseqs_pass_1"
tmp_2="$SLURM_TMPDIR/mmseqs_pass_2"

cluster_count="$step_02_dir/cluster_count.tsv"
splits_dir="$step_02_dir/splits"

seqtk seq -r -l 0 "$unique" | sed '/^>/s/$/r/' > "$reverse"
cat "$unique" "$reverse" > "$both"

mmseqs easy-cluster \
    "$both" \
    "$pass_1" \
    "$tmp_1" \
    -c "$coverage" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --kmer-per-seq "$kmer_per_seq" \
    --min-seq-id "$pid_pass_1" \
    --cov-mode "$cov_mode" \
    --filter-hits "$filter_hits"

mmseqs easy-cluster \
    "${pass_1}_rep_seq.fasta" \
    "$pass_2" \
    "$tmp_2" \
    -c "$coverage" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --min-seq-id "$pid_pass_2" \
    --cov-mode "$cov_mode" \
    --filter-hits "$filter_hits"

cut -f1 "${pass_2}_cluster.tsv" |
    sort |
    uniq -c |
    awk '$1 > 1 {print $1 "\t" $2}' \
    > "$cluster_count"

rm -rf "$splits_dir"
mkdir -p "$splits_dir"

awk -v dir="$splits_dir" '
    NR == FNR {
        wanted[$2] = 1
        next
    }

    {
        if (previous ~ /^>/ && $0 ~ /^>/) {
            if (output != "") close(output)

            cluster = substr(previous, 2)
            sub(/[ \t].*/, "", cluster)

            output = cluster in wanted \
                ? dir "/" cluster "_cluster.fasta" \
                : ""
        }
        else if (previous ~ /^>/ && output != "") {
            print previous > output
            print $0 > output
        }

        previous = $0
    }
' "$cluster_count" "${pass_2}_all_seqs.fasta"

date