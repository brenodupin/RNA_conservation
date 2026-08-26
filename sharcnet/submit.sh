#!/bin/bash

# Submit an RNA conservation pipeline step to Slurm.
#
# Usage:
#   sharcnet/submit.sh STEP PROJECT_DIR
#
# Example:
#   sharcnet/submit.sh 00 sandbox/data_protists_plt

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 STEP PROJECT_DIR" >&2
    exit 1
fi

step=$1
project_dir=$2

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

step_script="$sharcnet_dir/step_${step}.sh"
log_pattern="$logs_dir/step_${step}_%j.out"

mkdir -p "$logs_dir"

printf '\n'
printf 'Smith Lab RNA conservation.\n'
printf '  Date:     %s\n' "$(date)"
printf '  Step:     %s\n' "$step"
printf '  Project:  %s\n' "$project_dir"
printf '  Script:   %s\n' "$step_script"
printf '\n'

submission=$(
    sbatch \
        --parsable \
        --output="$log_pattern" \
        "$step_script" "$project_dir"
)

job_id=${submission%%;*}
log_file="$logs_dir/step_${step}_${job_id}.out"

printf '  ✓ Submitted job %s\n' "$job_id"
printf '  → Log: %s\n\n' "$log_file"
