## Holds information about the data directory and other path variables. 
## This file is sourced by other scripts to access these variables.

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
container_dir="$sharcnet_dir/containers"
rnatools="$container_dir/rnatools_v2.1.sif"

repo_dir=$(dirname "$sharcnet_dir")
step_01_scripts="$repo_dir/step1_createWindows"
step_05_scripts="$repo_dir/step5_evaluationOfRNAstructures"

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

# step 4 parameters (mlocarna via worker_04.sh, run as a Slurm array by
# step_04a.sh / step_04b.sh)
locarna_timeout="24h"
step04_jobs=20
step04_max_submit=900

# step 4a: small clusters
step04a_cluster_min=7
step04a_cluster_max=19
step04a_cpus=3
step04a_mem=6000M
step04a_time=1-01:00:00

# step 4b: larger clusters (scaled up from 4a; no direct profiling yet)
step04b_cluster_min=20
step04b_cluster_max=49
step04b_cpus=4
step04b_mem=16000M
step04b_time=1-01:00:00

# step 5 parameters (trimAlignment.pl and RNA-SCoRE via step_05a.sh, R-scape
# via step_05b.sh)

# RNA-SCoRE thresholds. Step 6 scores its homolog hits with RNA-SCoRE too, at
# its own thresholds, hence the step05 prefix.
step05_rnascore_mt=0.5      # fraction of the motif's base pairs a sequence must form
step05_rnascore_bp=0.75     # fraction of each stem's base pairs a sequence must form
step05_rnascore_gc=0.30     # fraction of a stem's base pairs that must be G:C or C:G
step05_rnascore_dupl=0      # 0 drops duplicate motif sequences, 1 keeps them

# RNA-SCoRE keeps one sequence out of each set of identical motifs and writes
# its rows in Perl's hash order, which is randomised per process: the rank is
# the same every time, the cleaned alignment is not. Pinning the seed makes
# reruns byte-identical, which is also what stops step_05b.sh from re-running
# R-scape on clusters whose alignment did not actually change. Any value works
# as long as it stays the same.
perl_hash_seed=42

# R-scape. Shared with steps 6 and 7, which run the same two-set test on their
# own alignments.
rscape_options="-s --cacofold"
rscape_seed=42
rscape_timeout="1h"

# Clusters RNA-SCoRE ranks at or above this are handed to step_05b.sh.
# High = 10 or more sequences passed, Mid = 7 to 9, Low = fewer (dropped).
#   allowed: High | Mid
rscape_min_rank=Mid

# multi-step parameters
fold_temperature=21

# Software modules
seqtk_module="seqtk/1.4"
mmseqs_module="mmseqs2/17-b804f"
clustalo_module="clustal-omega/1.2.4"
apptainer_module="apptainer/1.4.5"