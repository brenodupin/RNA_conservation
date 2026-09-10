#!/bin/bash

# Submit an RNA conservation pipeline step to Slurm.
#
# Usage:
#   sharcnet/submit.sh STEP DATA_DIR
#
# Example:
#   sharcnet/submit.sh 00 sandbox/data_protists_plt

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 STEP DATA_DIR" >&2
    exit 1
fi

step=$1
data_dir=$2

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

step_script="$sharcnet_dir/step_${step}.sh"

# Array steps need %A_%a so tasks do not all write to the same log.
if grep -q '^#SBATCH --array' "$step_script"; then
    is_array=1
    log_pattern="$logs_dir/step_${step}_%A_%a.out"
else
    is_array=0
    log_pattern="$logs_dir/step_${step}_%j.out"
fi

mkdir -p "$logs_dir"

printf '\n'
printf 'Smith Lab RNA conservation.\n'
printf '  Date:   %s\n' "$(date)"
printf '  Step:   %s\n' "$step"
printf '  Data:   %s\n' "$data_dir"
printf '  Script: %s\n' "$step_script"
printf '\n'

submission=$(
    sbatch \
        --parsable \
        --output="$log_pattern" \
        "$step_script" \
        "$data_dir" \
        "$sharcnet_dir"
)

job_id=${submission%%;*}

if ((is_array)); then
    log_file="$logs_dir/step_${step}_${job_id}_*.out"
else
    log_file="$logs_dir/step_${step}_${job_id}.out"
fi

printf '  ✓ Submitted job %s\n' "$job_id"
printf '  → Log: %s\n\n' "$log_file"