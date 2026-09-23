#!/bin/bash
#SBATCH --job-name=rnac_05b
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Test the motifs step 05a kept for covariation, with R-scape's two-set test
# and CaCoFold. Outputs are written next to the alignment they came from, in
# DATA_DIR/05_evaluation/<cluster>.
#
# The clusters tested are the ones RNA-SCoRE ranked at or above
# rscape_min_rank (info.sh): High only, or High and Mid as upstream does.
# --list FILE runs a hand-picked set instead, one cluster name per line, blank
# lines and lines starting with # ignored -- copy high_mid_clusters.txt, cut it
# down after reading step_05_status.tsv, and pass it back here.
#
# A cluster is done when it has a <cluster>_result_05b.txt marker. As in 05a
# the marker records the md5 of the alignment it was built from, here
# <cluster>_motif_cleaned.sto, so rerunning 05a does not send unchanged
# clusters back through R-scape: with perl_hash_seed pinned, RNA-SCoRE writes
# the same cleaned alignment every time. --force redoes everything.
#
# A run without --list also clears the R-scape output of any cluster that is no
# longer ranked high enough, so refolding one in step 4 cannot leave a result
# behind that describes an alignment that no longer exists. That includes
# clusters tested earlier through --list.
#
# Per cluster, in DATA_DIR/05_evaluation/<cluster>:
#   <cluster>_rscape/               everything R-scape produced, including the
#                                   .power counts and the R2R drawings of the
#                                   given and the CaCoFold structure
#   <cluster>_rscape.log            R-scape's own output
#   <cluster>_rscape_covariation.tsv  the counts for this cluster
#   <cluster>_result_05b.txt        the marker
#
# And for the run as a whole, in DATA_DIR/05_evaluation:
#   rscape_covariation.tsv          every cluster's counts, most covariation
#                                   first -- step 6 keeps motifs with at least
#                                   5% of their base pairs covarying
#   failed_05b.tsv                  clusters with no result, and why
#   step_05_status.tsv              one row per cluster (see status_05.sh)
#
# Usage:
#   sharcnet/submit.sh 05b DATA_DIR [--force] [--list FILE]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force] [--list FILE]" >&2
    exit 1
fi

data_dir=$1
script_dir=$2
shift 2

force=no
list_file=

while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            force=yes
            ;;
        --list)
            [[ $# -ge 2 ]] || {
                echo "--list needs a file" >&2
                exit 1
            }
            list_file=$2
            shift
            ;;
        *)
            echo "Unknown argument: $1 (expected --force or --list FILE)" >&2
            exit 1
            ;;
    esac
    shift
done

echo " step_05b.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

rscape_worker="$script_dir/worker_rscape.sh"
status_script="$script_dir/status_05.sh"

[[ -d "$step_05_dir" ]] || {
    echo "Missing step 5 output: $step_05_dir -- run step 05a first." >&2
    exit 1
}

for script in "$rscape_worker" "$status_script"; do
    [[ -x "$script" ]] || {
        echo "Not executable: $script" >&2
        exit 1
    }
done

case $rscape_min_rank in
    High | Mid) ;;
    *)
        echo "rscape_min_rank must be High or Mid, got '$rscape_min_rank'" >&2
        exit 1
        ;;
esac

module load "$apptainer_module"

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

threads=${SLURM_CPUS_PER_TASK:-1}

settings=$(printf 'options=%s seed=%s' "$rscape_options" "$rscape_seed")

# ---------------------------------------------------------------------------
# which clusters to test
# ---------------------------------------------------------------------------
#
# From the 05a markers rather than from allClusters_evaluation.tsv: the markers
# are what 05a itself treats as the record of a finished cluster, and they
# cannot have been edited by hand between the two steps.

declare -A rank=()

while IFS=$'\t' read -r cluster field value; do
    [[ "$field" == rank ]] && rank[$cluster]=$value
done < <(
    find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_05a.txt' -print0 |
        xargs -0 -r awk -F '\t' '
            FNR == 1 {
                split(FILENAME, path, "/")
                cluster = path[length(path) - 1]
            }

            $1 == "rank" {
                print cluster "\t" $1 "\t" $2
            }
        '
)

candidates=()

for cluster in "${!rank[@]}"; do
    case ${rank[$cluster]} in
        High)
            candidates+=("$cluster")
            ;;
        Mid)
            [[ "$rscape_min_rank" == Mid ]] && candidates+=("$cluster")
            ;;
    esac
done

if ((${#candidates[@]} > 0)); then
    mapfile -t candidates < <(printf '%s\n' "${candidates[@]}" | LC_ALL=C sort)
fi

candidate_total=${#candidates[@]}

# What this run tests: every candidate, or the hand-picked list.
if [[ -n "$list_file" ]]; then
    [[ -s "$list_file" ]] || {
        echo "Cluster list not found or empty: $list_file" >&2
        exit 1
    }

    mapfile -t run_set < <(
        grep -vE '^[[:space:]]*(#.*)?$' "$list_file" |
            awk '{ print $1 }' |
            LC_ALL=C sort -u
    )
else
    run_set=("${candidates[@]}")
fi

if ((${#run_set[@]} == 0)); then
    echo \
        "No cluster ranked $rscape_min_rank or better in $step_05_dir" \
        "-- nothing to test." \
        >&2
    exit 0
fi

# Candidates and listed clusters together. The failure list and the covariation
# table are rebuilt from this, so a run with --list reports on the clusters it
# tested without dropping what a full run had already recorded.
mapfile -t considered < <(
    printf '%s\n' "${candidates[@]}" "${run_set[@]}" | LC_ALL=C sort -u
)

# ---------------------------------------------------------------------------
# what still needs doing
# ---------------------------------------------------------------------------

alignments=()
missing=()

for cluster in "${run_set[@]}"; do
    cleaned="$step_05_dir/$cluster/${cluster}_motif_cleaned.sto"

    if [[ -s "$cleaned" ]]; then
        alignments+=("$cleaned")
    else
        missing+=("$cluster")
    fi
done

declare -A input_md5=()

if ((${#alignments[@]} > 0)); then
    while read -r sum path; do
        cluster=${path%/*}
        cluster=${cluster##*/}
        input_md5[$cluster]=$sum
    done < <(printf '%s\0' "${alignments[@]}" | xargs -0 -r md5sum)
fi

declare -A marker_md5=()
declare -A marker_settings=()

while IFS=$'\t' read -r cluster field value; do
    case $field in
        input_md5)
            marker_md5[$cluster]=$value
            ;;
        settings)
            marker_settings[$cluster]=$value
            ;;
    esac
done < <(
    find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_05b.txt' -print0 |
        xargs -0 -r awk -F '\t' '
            FNR == 1 {
                split(FILENAME, path, "/")
                cluster = path[length(path) - 1]
            }

            $1 == "input_md5" || $1 == "settings" {
                print cluster "\t" $1 "\t" $2
            }
        '
)

todo=()
unchanged=0

for cluster in "${run_set[@]}"; do
    [[ -n "${input_md5[$cluster]:-}" ]] || continue

    if [[ "$force" != yes &&
        "${marker_md5[$cluster]:-}" == "${input_md5[$cluster]}" &&
        "${marker_settings[$cluster]:-}" == "$settings" ]]
    then
        unchanged=$((unchanged + 1))
        continue
    fi

    todo+=("$cluster")
done

todo_total=${#todo[@]}

if [[ -n "$list_file" ]]; then
    printf 'Step 05b: %s cluster(s) listed, %s already done, %s to test.\n' \
        "${#run_set[@]}" "$unchanged" "$todo_total"
    printf '  List:     %s (%s ranked %s or better in all)\n' \
        "$list_file" "$candidate_total" "$rscape_min_rank"
else
    printf 'Step 05b: %s cluster(s) ranked %s or better, %s already done, %s to test.\n' \
        "$candidate_total" "$rscape_min_rank" "$unchanged" "$todo_total"
fi

printf '  Workers:  %s\n' "$threads"
printf '  R-scape:  %s (seed %s, timeout %s)\n' \
    "$rscape_options" "$rscape_seed" "$rscape_timeout"

printf '\n'

if ((${#missing[@]} > 0)); then
    printf '  %s cluster(s) have no _motif_cleaned.sto and are skipped.\n\n' \
        "${#missing[@]}"
fi

# ---------------------------------------------------------------------------
# run R-scape
# ---------------------------------------------------------------------------

if ((todo_total > 0)); then
    for cluster in "${todo[@]}"; do
        printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" \
            "$step_05_dir/$cluster" \
            "${cluster}_motif_cleaned.sto" \
            "${cluster}_rscape" \
            "${cluster}_rscape" \
            "$rscape_seed" \
            "$rscape_timeout" \
            "$rscape_options"
    done |
        xargs -0 -r -P "$threads" -n8 "$rscape_worker" || true
fi

# ---------------------------------------------------------------------------
# markers and failures
# ---------------------------------------------------------------------------
#
# As in 05a, what happened is read back off disk: the covariation file is the
# last thing worker_rscape.sh writes, so a cluster that has one is done. The
# markers are written for what this run tested; the failure list then covers
# every cluster considered, which is how an earlier failure survives a run with
# --list. A candidate this run did not touch and that never ran is neither done
# nor failed, so it is left out of the list and shows up as pending in
# step_05_status.tsv.

failed="$step_05_dir/failed_05b.tsv"
finished_at=$(date '+%Y-%m-%d %H:%M:%S')

for cluster in "${todo[@]}"; do
    cluster_dir="$step_05_dir/$cluster"
    covariation="$cluster_dir/${cluster}_rscape_covariation.tsv"
    marker="$cluster_dir/${cluster}_result_05b.txt"

    if [[ ! -s "$covariation" ]]; then
        rm -f "$marker"
        continue
    fi

    read -r cluster_bpairs cluster_covarying cluster_percent < <(
        awk -F '\t' '{ printf "%s %s %s\n", $2, $5, $6 }' "$covariation"
    )

    {
        printf 'cluster\t%s\n' "$cluster"
        printf 'input\t%s\n' "$cluster_dir/${cluster}_motif_cleaned.sto"
        printf 'input_md5\t%s\n' "${input_md5[$cluster]}"
        printf 'settings\t%s\n' "$settings"
        printf 'bpairs\t%s\n' "$cluster_bpairs"
        printf 'covarying\t%s\n' "$cluster_covarying"
        printf 'percent_covarying\t%s\n' "$cluster_percent"
        printf 'finished\t%s\n' "$finished_at"
    } > "$marker.tmp"

    mv "$marker.tmp" "$marker"
done

declare -A attempted=()

for cluster in "${todo[@]}"; do
    attempted[$cluster]=yes
done

: > "$failed.tmp"

for cluster in "${considered[@]}"; do
    cluster_dir="$step_05_dir/$cluster"

    [[ -s "$cluster_dir/${cluster}_result_05b.txt" ]] && continue

    if [[ ! -s "$cluster_dir/${cluster}_motif_cleaned.sto" ]]; then
        printf '%s\t%s\n' \
            "$cluster" \
            "no _motif_cleaned.sto, nothing for R-scape to test" \
            >> "$failed.tmp"
    elif [[ -n "${attempted[$cluster]:-}" || -f "$cluster_dir/${cluster}_rscape.log" ]]; then
        printf '%s\t%s\n' \
            "$cluster" \
            "R-scape failed, see ${cluster}_rscape.log" \
            >> "$failed.tmp"
    fi
done

LC_ALL=C sort "$failed.tmp" > "$failed"
rm -f "$failed.tmp"

failed_total=$(wc -l < "$failed")

# ---------------------------------------------------------------------------
# clusters that are no longer candidates
# ---------------------------------------------------------------------------
#
# A cluster refolded by step 4 and rescored by 05a can come back ranked Low.
# Its old R-scape output describes an alignment that no longer exists, so it is
# removed rather than left to be read as a current result.
#
# Only when the candidates came from the markers: with --list every cluster
# outside that hand-picked list would otherwise look retired.

stale=0

if [[ -z "$list_file" ]]; then
    declare -A is_candidate=()

    for cluster in "${candidates[@]}"; do
        is_candidate[$cluster]=yes
    done

    while IFS= read -r -d '' marker; do
        cluster_dir=${marker%/*}
        cluster=${cluster_dir##*/}

        [[ -n "${is_candidate[$cluster]:-}" ]] && continue

        rm -rf "$cluster_dir/${cluster}_rscape"
        rm -f \
            "$marker" \
            "$cluster_dir/${cluster}_rscape.log" \
            "$cluster_dir/${cluster}_rscape_covariation.tsv"

        stale=$((stale + 1))
    done < <(
        find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_05b.txt' -print0
    )
fi

# ---------------------------------------------------------------------------
# collated table
# ---------------------------------------------------------------------------

covariation_table="$step_05_dir/rscape_covariation.tsv"

{
    printf 'cluster\tbpairs\texpected\texpected_sd\tcovarying\tpercent_covarying\n'

    for cluster in "${considered[@]}"; do
        covariation="$step_05_dir/$cluster/${cluster}_rscape_covariation.tsv"
        [[ -s "$step_05_dir/$cluster/${cluster}_result_05b.txt" ]] || continue
        [[ -s "$covariation" ]] || continue
        printf '%s\0' "$covariation"
    done |
        xargs -0 -r cat |
        awk -F '\t' 'BEGIN { OFS = "\t" } { sub(/_motif_cleaned\.sto$/, "", $1); print }' |
        LC_ALL=C sort -t $'\t' -k6,6gr -k1,1
} > "$covariation_table.tmp"

mv "$covariation_table.tmp" "$covariation_table"

tested=$(($(wc -l < "$covariation_table") - 1))

read -r with_covariation above_five < <(
    awk -F '\t' '
        NR == 1 {
            next
        }

        {
            if ($5 > 0) covarying++
            if ($6 >= 5) strong++
        }

        END {
            printf "%d %d\n", covarying, strong
        }
    ' "$covariation_table"
)

"$status_script" "$data_dir" > /dev/null ||
    echo "status_05.sh failed, step_05_status.tsv may be out of date" >&2

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

printf '\n'
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "ranked $rscape_min_rank or better"   "$candidate_total"
printf '%-42s %10s\n' "already done, unchanged"             "$unchanged"
printf '%-42s %10s\n' "tested this run"                     "$todo_total"
printf '%-42s %10s\n' "with R-scape results"                "$tested"
printf '%-42s %10s\n' "  any base pair covarying"           "$with_covariation"
printf '%-42s %10s\n' "  at least 5% covarying"             "$above_five"
printf '%-42s %10s\n' "failed"                              "$failed_total"

if ((stale > 0)); then
    printf '%-42s %10s\n' "removed, no longer ranked high enough" "$stale"
fi

printf '\n'
printf '  Covariation: %s\n' "$covariation_table"
printf '  Status:      %s\n' "$step_05_dir/step_05_status.tsv"

if ((failed_total > 0)); then
    printf '  Failures:    %s (%s clusters)\n' "$failed" "$failed_total"
fi

printf '\n'
date

# Slurm should mark the job failed when clusters did.
((failed_total == 0))
