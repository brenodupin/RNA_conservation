#!/bin/bash

# Submit an RNA conservation pipeline step to Slurm.
#
# Any extra arguments after DATA_DIR are forwarded to the step script as-is
# -- e.g. step_04a.sh/step_04b.sh accept a trailing --force to redo
# clusters that already have a result.
#
# Usage:
#   sharcnet/submit.sh STEP DATA_DIR [EXTRA_ARGS...]
#
# Example:
#   sharcnet/submit.sh 00 sandbox/data_protists_plt
#   sharcnet/submit.sh 04a data_protists_plt_July_31 --force

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 STEP DATA_DIR [EXTRA_ARGS...]" >&2
    exit 1
fi

step=$1
data_dir=$2
shift 2

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

step_script="$sharcnet_dir/step_${step}.sh"
log_pattern="$logs_dir/step_${step}_%j.out"

mkdir -p "$logs_dir"

printf '\n'
printf 'Smith Lab RNA conservation.\n'
printf '  Date:   %s\n' "$(date)"
printf '  Step:   %s\n' "$step"
printf '  Data:   %s\n' "$data_dir"
printf '  Script: %s\n' "$step_script"
if [[ $# -gt 0 ]]; then
    printf '  Extra:  %s\n' "$*"
fi
printf '\n'

submission=$(
    sbatch \
        --parsable \
        --output="$log_pattern" \
        "$step_script" \
        "$data_dir" \
        "$sharcnet_dir" \
        "$@"
)

job_id=${submission%%;*}
log_file="$logs_dir/step_${step}_${job_id}.out"

printf '  ✓ Submitted job %s\n' "$job_id"
printf '  → Log: %s\n\n' "$log_file"
