## Holds information about the data directory and other path variables. 
## This file is sourced by other scripts to access these variables.

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
container_dir="$sharcnet_dir/containers"
rnatools="$container_dir/rnatools_v2.1.sif"

repo_dir=$(dirname "$sharcnet_dir")
step_01_scripts="$repo_dir/step1_createWindows"

# Resolve the data directory from either:
#   1. The path supplied directly.
#   2. $SCRATCH/<name> on Alliance systems.

if [[ -d "$data_dir" ]]; then
    data_dir=$(realpath "$data_dir")
elif [[ -n "${SCRATCH:-}" && -d "$SCRATCH/$data_dir" ]]; then
    data_dir=$(realpath "$SCRATCH/$data_dir")
else
    echo "Data directory not found: $data_dir" >&2
    return 1
fi

input_dir="$data_dir/input"
step_00_dir="$data_dir/00_oneline"
step_01_dir="$data_dir/01_windows"
step_02_dir="$data_dir/02_clusters"
step_03_dir="$data_dir/03_screen"
step_04_dir="$data_dir/04_locarna"
step_05_dir="$data_dir/05_evaluation"

logs_dir="$data_dir/logs"

# step 1 parameters
window_size=250
overlap=75
windows_prefix="all_windows"

# step 2 parameters
pid_pass_1=0.95
pid_pass_2=0.50
coverage=0.80
kmer_per_seq=200
cov_mode=0
filter_hits=1

# step 4 parameters
locarna_jobs=10
locarna_threads=3
locarna_timeout="24h"

# multi-step parameters
fold_temperature=21

# Software modules
seqtk_module="seqtk/1.4"
mmseqs_module="mmseqs2/17-b804f"
clustalo_module="clustal-omega/1.2.4"
apptainer_module="apptainer/1.4.5"