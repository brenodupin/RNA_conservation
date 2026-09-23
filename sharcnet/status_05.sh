#!/bin/bash

# Where every cluster stands in step 5: what finished, how RNA-SCoRE ranked it,
# and how much covariation R-scape found.
#
# Everything comes from the per-cluster markers and the two failure lists, so
# this can be run at any time -- including while a step_05a.sh or step_05b.sh
# job is still going, where it shows the clusters done so far. step_05a.sh and
# step_05b.sh both refresh it when they finish.
#
# Writes DATA_DIR/05_evaluation/step_05_status.tsv:
#   cluster            cluster name
#   members            sequences in the cluster (step 2)
#   step_05a           done, failed, or pending
#   rank               High, Mid or Low
#   passed             sequences that passed RNA-SCoRE
#   step_05b           done, failed, pending, or - when the rank is too low
#   bpairs             base pairs in the motif structure
#   covarying          base pairs R-scape found covarying
#   percent_covarying  covarying as a percentage of bpairs
#
# Usage:
#   sharcnet/status_05.sh DATA_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 DATA_DIR" >&2
    exit 1
fi

data_dir=$1

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

[[ -d "$step_05_dir" ]] || {
    echo "Missing step 5 output: $step_05_dir -- run step 05a first." >&2
    exit 1
}

status="$step_05_dir/step_05_status.tsv"

# Clusters to report on: the ones step 4 folded, which is also what step_05a.sh
# works from.
mapfile -t clusters < <(
    find "$step_04_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result.stk' -printf '%P\n' |
        sed 's|/.*||' |
        LC_ALL=C sort -u
)

declare -A members=()

if [[ -s "$step_02_dir/cluster_count.tsv" ]]; then
    while IFS=$'\t' read -r count cluster; do
        members[$cluster]=$count
    done < "$step_02_dir/cluster_count.tsv"
fi

# Marker fields, one awk pass over each sub-step's markers.
declare -A rank=()
declare -A passed=()
declare -A bpairs=()
declare -A covarying=()
declare -A percent=()

read_markers() {  # read_markers <marker glob>
    find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f -name "$1" -print0 |
        xargs -0 -r awk -F '\t' '
            FNR == 1 {
                split(FILENAME, path, "/")
                cluster = path[length(path) - 1]
            }

            NF == 2 {
                print cluster "\t" $1 "\t" $2
            }
        '
}

while IFS=$'\t' read -r cluster field value; do
    case $field in
        rank)   rank[$cluster]=$value ;;
        passed) passed[$cluster]=$value ;;
    esac
done < <(read_markers '*_result_05a.txt')

while IFS=$'\t' read -r cluster field value; do
    case $field in
        bpairs)            bpairs[$cluster]=$value ;;
        covarying)         covarying[$cluster]=$value ;;
        percent_covarying) percent[$cluster]=$value ;;
    esac
done < <(read_markers '*_result_05b.txt')

declare -A failed_05a=()
declare -A failed_05b=()

if [[ -s "$step_05_dir/failed_05a.tsv" ]]; then
    while IFS=$'\t' read -r cluster _; do
        failed_05a[$cluster]=yes
    done < "$step_05_dir/failed_05a.tsv"
fi

if [[ -s "$step_05_dir/failed_05b.tsv" ]]; then
    while IFS=$'\t' read -r cluster _; do
        failed_05b[$cluster]=yes
    done < "$step_05_dir/failed_05b.tsv"
fi

done_05a=0
error_05a=0
pending_05a=0
done_05b=0
error_05b=0
pending_05b=0
not_ranked=0

{
    printf 'cluster\tmembers\tstep_05a\trank\tpassed'
    printf '\tstep_05b\tbpairs\tcovarying\tpercent_covarying\n'

    for cluster in "${clusters[@]}"; do
        cluster_rank=${rank[$cluster]:-}

        if [[ -n "$cluster_rank" ]]; then
            state_05a="done"
            done_05a=$((done_05a + 1))
        elif [[ -n "${failed_05a[$cluster]:-}" ]]; then
            state_05a="failed"
            error_05a=$((error_05a + 1))
        else
            state_05a="pending"
            pending_05a=$((pending_05a + 1))
        fi

        # Only High, or High and Mid, reach step 5b (rscape_min_rank).
        candidate=no

        if [[ "$cluster_rank" == High ]]; then
            candidate=yes
        elif [[ "$cluster_rank" == Mid && "$rscape_min_rank" == Mid ]]; then
            candidate=yes
        fi

        if [[ -n "${bpairs[$cluster]:-}" ]]; then
            state_05b="done"
            done_05b=$((done_05b + 1))
        elif [[ "$candidate" == no ]]; then
            state_05b="-"
            not_ranked=$((not_ranked + 1))
        elif [[ -n "${failed_05b[$cluster]:-}" ]]; then
            state_05b="failed"
            error_05b=$((error_05b + 1))
        else
            state_05b="pending"
            pending_05b=$((pending_05b + 1))
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$cluster" \
            "${members[$cluster]:--}" \
            "$state_05a" \
            "${cluster_rank:--}" \
            "${passed[$cluster]:--}" \
            "$state_05b" \
            "${bpairs[$cluster]:--}" \
            "${covarying[$cluster]:--}" \
            "${percent[$cluster]:--}"
    done
} > "$status.tmp"

mv "$status.tmp" "$status"

printf 'Step 05 status: %s clusters folded by step 4.\n' "${#clusters[@]}"
printf '  05a: %s done, %s failed, %s pending\n' \
    "$done_05a" "$error_05a" "$pending_05a"
printf '  05b: %s done, %s failed, %s pending, %s below %s\n' \
    "$done_05b" "$error_05b" "$pending_05b" "$not_ranked" "$rscape_min_rank"
printf '  Table: %s\n' "$status"
