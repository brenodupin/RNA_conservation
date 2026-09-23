#!/bin/bash
#SBATCH --job-name=rnac_06a
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000M
#SBATCH --array=1

# Build a covariance model from one cluster's seed alignment, calibrate it, and
# search the homolog database with it. Runs as one task of a Slurm array;
# step_06a.sh is the normal way to submit it, with the resources from info.sh
# on its own sbatch call (the #SBATCH lines above are only placeholders).
#
# The manifest is one line per cluster: "<cluster>\t<seed alignment>\t<md5>".
# Array task N processes line N. cmcalibrate and cmsearch use
# SLURM_CPUS_PER_TASK threads.
#
# Calibration is the slow part and depends only on the seed, so the calibrated
# model is kept and reused when a later run only changes the database or the
# search settings: <cluster>_cm.txt records what the model was built from.
#
# A cluster is done when it has a <cluster>_result_06a.txt marker recording the
# seed, the database and the search settings; a task whose marker already
# matches is a cheap no-op, so a failed index can be resubmitted directly:
#
#   sbatch --array=3,17 --cpus-per-task=8 --mem=8000M --time=12:00:00 \
#       --output=DATA_DIR/logs/step_06a_r1_Sep_23/%A_%a.out \
#       SHARCNET_DIR/worker_06a.sh MANIFEST DATA_DIR SHARCNET_DIR \
#       SEARCH_DB DB_MD5 'SETTINGS'
#
# (step_06a.sh prints the exact command, settings included, for its run.) A
# task that fails writes <cluster>_failed_06a.txt with the reason.
#
# Usage:
#   sbatch [--array=... --cpus-per-task=N --mem=... --time=...] \
#       worker_06a.sh MANIFEST DATA_DIR SHARCNET_DIR SEARCH_DB DB_MD5 SETTINGS [--force]

set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
    echo \
        "Usage: sbatch $0 MANIFEST DATA_DIR SHARCNET_DIR SEARCH_DB DB_MD5 SETTINGS [--force]" \
        >&2
    exit 1
fi

manifest=$1
data_dir=$2
script_dir=$3
search_db=$4
db_md5=$5
settings=$6
force=${7:-}

source "$script_dir/info.sh"

module load "$apptainer_module"

: "${SLURM_ARRAY_TASK_ID:?This script must be submitted as a job array}"

[[ -s "$manifest" ]] || {
    echo "Missing manifest: $manifest" >&2
    exit 1
}

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$manifest")

[[ -n "$line" ]] || {
    echo \
        "No manifest entry for array task $SLURM_ARRAY_TASK_ID" \
        "(manifest has $(wc -l < "$manifest") lines)" \
        >&2
    exit 1
}

IFS=$'\t' read -r cluster seed _ <<< "$line"

threads=${SLURM_CPUS_PER_TASK:-1}

round_dir="$step_06_dir/round_$step06_round"
cluster_dir="$round_dir/$cluster"
marker="$cluster_dir/${cluster}_result_06a.txt"
failed_note="$cluster_dir/${cluster}_failed_06a.txt"

seed_copy="${cluster}_seed.sto"
model="${cluster}.cm"
model_note="$cluster_dir/${cluster}_cm.txt"
hits_text="${cluster}_hits.txt"
hits_table="${cluster}_hits.tbl"

echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Cluster:    $cluster"
echo "Seed:       $seed"
echo "Database:   $search_db"
echo "Threads:    $threads"

mkdir -p "$cluster_dir"

fail() {
    echo "$1" >&2
    printf '%s\n' "$1" > "$failed_note"
    exit 1
}

# Anything else that goes wrong still leaves a note saying where.
trap 'fail "worker_06a.sh stopped at line $LINENO for $cluster, see the task log"' ERR

marker_field() {  # marker_field <field>
    awk -F '\t' -v field="$1" '$1 == field { print $2 }' "$marker"
}

[[ -s "$seed" ]] || fail "seed alignment not found: $seed"

seed_md5=$(md5sum < "$seed" | cut -d ' ' -f 1)

# Nothing changed since the last successful run: nothing to do.
if [[ "$force" != "--force" && -s "$marker" &&
    "$(marker_field seed_md5)" == "$seed_md5" &&
    "$(marker_field db_md5)" == "$db_md5" &&
    "$(marker_field settings)" == "$settings" ]]
then
    echo "Already completed, skipping."
    exit 0
fi

[[ -s "$rnatools" ]] || fail "container not found: $rnatools"
[[ -s "$search_db" ]] || fail "search database not found: $search_db"

rm -f "$marker" "$failed_note"
rm -f \
    "$cluster_dir/$hits_text" "$cluster_dir/$hits_text.tmp" \
    "$cluster_dir/$hits_table" "$cluster_dir/$hits_table.tmp"

cp "$seed" "$cluster_dir/$seed_copy"

# The cluster directory is the working directory; the database is mounted
# read-only next to it.
db_dir=$(dirname "$search_db")
db_name=$(basename "$search_db")

in_container() {
    apptainer exec \
        --bind "$cluster_dir:/work" \
        --bind "$db_dir:/db:ro" \
        --pwd /work \
        "$rnatools" \
        "$@"
}

# ---------------------------------------------------------------------------
# model: cmbuild + cmcalibrate, unless an identical one is already calibrated
# ---------------------------------------------------------------------------

model_key="seed_md5=$seed_md5 cmcalibrate=$step06_cmcalibrate_options"
build_seconds=0
calibrate_seconds=0

if [[ "$force" != "--force" && -s "$cluster_dir/$model" && -s "$model_note" &&
    "$(cat "$model_note")" == "$model_key" ]]
then
    echo "Reusing the calibrated model built from this seed."
else
    rm -f "$cluster_dir/$model" "$cluster_dir/$model.tmp" "$model_note"

    start=$SECONDS

    in_container cmbuild -F "$model.tmp" "$seed_copy" \
        > "$cluster_dir/${cluster}_cmbuild.log" 2>&1 ||
        fail "cmbuild failed for $cluster, see ${cluster}_cmbuild.log"

    build_seconds=$((SECONDS - start))
    start=$SECONDS

    # shellcheck disable=SC2086 # a list of cmcalibrate options, possibly empty
    in_container cmcalibrate --cpu "$threads" $step06_cmcalibrate_options "$model.tmp" \
        > "$cluster_dir/${cluster}_cmcalibrate.log" 2>&1 ||
        fail "cmcalibrate failed for $cluster, see ${cluster}_cmcalibrate.log"

    calibrate_seconds=$((SECONDS - start))

    mv "$cluster_dir/$model.tmp" "$cluster_dir/$model"
    printf '%s\n' "$model_key" > "$model_note"
fi

# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

start=$SECONDS

# shellcheck disable=SC2086 # a list of cmsearch options, possibly empty
in_container cmsearch \
    --cpu "$threads" \
    $step06_cmsearch_options \
    -E "$step06_cmsearch_evalue" \
    --tblout "$hits_table.tmp" \
    "$model" \
    "/db/$db_name" \
    > "$cluster_dir/$hits_text.tmp" 2> "$cluster_dir/${cluster}_cmsearch.err" ||
    fail "cmsearch failed for $cluster, see ${cluster}_cmsearch.err"

search_seconds=$((SECONDS - start))

mv "$cluster_dir/$hits_text.tmp" "$cluster_dir/$hits_text"
mv "$cluster_dir/$hits_table.tmp" "$cluster_dir/$hits_table"

[[ -s "$cluster_dir/${cluster}_cmsearch.err" ]] ||
    rm -f "$cluster_dir/${cluster}_cmsearch.err"

# The hits without cmsearch's comment lines, which carry the date: the md5 of
# what was found, so step 6b can tell a real change from a rerun that found the
# same thing.
# (grep finds nothing when there are no hits, which is a result, not an error.)
hits=$(grep -vc '^#' "$cluster_dir/$hits_table" || true)
hits_md5=$(
    { grep -v '^#' "$cluster_dir/$hits_table" || true; } | md5sum | cut -d ' ' -f 1
)

{
    printf 'cluster\t%s\n' "$cluster"
    printf 'round\t%s\n' "$step06_round"
    printf 'seed\t%s\n' "$seed"
    printf 'seed_md5\t%s\n' "$seed_md5"
    printf 'search_db\t%s\n' "$search_db"
    printf 'db_md5\t%s\n' "$db_md5"
    printf 'settings\t%s\n' "$settings"
    printf 'hits\t%s\n' "$hits"
    printf 'hits_md5\t%s\n' "$hits_md5"
    printf 'build_seconds\t%s\n' "$build_seconds"
    printf 'calibrate_seconds\t%s\n' "$calibrate_seconds"
    printf 'search_seconds\t%s\n' "$search_seconds"
    printf 'threads\t%s\n' "$threads"
    printf 'finished\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
} > "$marker.tmp"

mv "$marker.tmp" "$marker"

echo "Completed $cluster: $hits hits (build ${build_seconds}s," \
    "calibrate ${calibrate_seconds}s, search ${search_seconds}s)"

date
