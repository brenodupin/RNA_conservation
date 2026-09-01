#!/bin/bash
#SBATCH --job-name=rnac_01
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Create overlapping windows from DATA_DIR/00_oneline, merge them, remove
# windows containing Ns or poly-nucleotide sequences, and remove duplicates.
# Outputs are written to DATA_DIR/01_windows.
#
# Usage:
#   sharcnet/submit.sh 01 DATA_DIR

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR" >&2
    exit 1
fi

data_dir=$1
script_dir=$2

echo " step_01.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

mkdir -p "$step_01_dir"

for fasta in "$step_00_dir"/*_oneLine.fasta; do
    name=$(basename "$fasta" _oneLine.fasta)

    perl "$step_01_scripts/createWindows.pl" \
        -f "$fasta" \
        -w "$window_size" \
        -p "$overlap" \
        -O "$step_01_dir/${name}_windows.fa"
done

combined="$step_01_dir/${windows_prefix}.fasta"
processed="$step_01_dir/${windows_prefix}_processed"
composition="$step_01_dir/${windows_prefix}_nuclComposition.tsv"

cat "$step_01_dir"/*_windows.fa > "$combined"

perl "$step_01_scripts/removeNs_polyN_windows.pl" "$combined" "$processed" > "$composition"

perl "$step_01_scripts/removeDuplicates.pl" "${processed}_noNs_polyN"

unique="${processed}_noNs_polyN_uniq.fasta"

paste - - < "$unique" | LC_ALL=C sort -k1,1 | tr '\t' '\n' > "${unique}.sorted"
mv "${unique}.sorted" "$unique"

date
