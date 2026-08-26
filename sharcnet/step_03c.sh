#!/bin/bash
#SBATCH --job-name=rnac_03c
#SBATCH --time=00:15:00
#SBATCH --mem=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Collect clusters whose RNALalifold output contains a predicted stem.
# The passing cluster names are written to RNALalifold_passedList.txt.
#
# Usage:
#   sharcnet/submit.sh 03c DATA_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 DATA_DIR" >&2
    exit 1
fi

data_dir=$1

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/info.sh"

passed_list="$step_03_dir/RNALalifold_passedList.txt"

shopt -s nullglob
outputs=("$step_03_dir"/*/*_RNALalifold.out)

{
    for output in "${outputs[@]}"; do
        if grep -qF '((' "$output"; then
            basename "$output" _RNALalifold.out
        fi
    done
} | LC_ALL=C sort > "$passed_list"

date