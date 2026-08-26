#!/bin/bash
#SBATCH --job-name=rnac_00
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Convert multiline FASTA files from DATA_DIR/input into one-line FASTA
# files in DATA_DIR/00_oneline.
#
# Usage:
#   sbatch sharcnet/step_00.sh DATA_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 DATA_DIR" >&2
    exit 1
fi

data_dir=$1

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/info.sh"

mkdir -p "$step_00_dir"

shopt -s nullglob

for fasta in "$input_dir"/*.fasta; do
    name=$(basename "$fasta")
    name=${name%.*}

    awk '
        /^>/ {
            if (seen) printf "\n"
            print
            seen = 1
            next
        }
        {
            gsub(/[[:space:]]/, "")
            printf "%s", $0
        }
        END {
            printf "\n"
        }
    ' "$fasta" > "$step_00_dir/${name}_oneLine.fasta"
done

date
