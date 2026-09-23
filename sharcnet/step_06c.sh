#!/bin/bash
#SBATCH --job-name=rnac_06c
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Test the hits alignments step 06b kept for covariation, with R-scape's
# two-set test and CaCoFold -- the same test step 05b runs on the seeds, now on
# the seed plus its homologs. Works on the current round (step06_round in
# info.sh), in DATA_DIR/06_homolog/round_<n>/<cluster>.
#
# The clusters tested are the ones whose hits RNA-SCoRE ranked at or above
# step06_min_rank (info.sh). --list FILE runs a hand-picked set instead, one
# cluster name per line, blank lines and lines starting with # ignored.
#
# A cluster is done when it has a <cluster>_result_06c.txt marker recording
# the md5 of <cluster>_hits_cleaned.sto, so rerunning 06b does not send
# unchanged clusters back through R-scape. --force redoes everything. A run
# without --list also clears the R-scape output of any cluster no longer ranked
# high enough, including clusters tested earlier through --list.
#
# --optional-outputs also renders the cleaned hits alignment and CaCoFold's
# alignment as HTML (stockholm_to_html.pl), for comparing the two structures
# by eye as the step 6 README describes. The R2R drawings of both structures
# are always there, in <cluster>_hits_rscape/.
#
# Per cluster, in DATA_DIR/06_homolog/round_<n>/<cluster>:
#   <cluster>_hits_rscape/              everything R-scape produced
#   <cluster>_hits_rscape.log           R-scape's own output
#   <cluster>_hits_rscape_covariation.tsv  the counts for this cluster
#   <cluster>_result_06c.txt            the marker
#
# And for the round as a whole, in DATA_DIR/06_homolog/round_<n>:
#   hits_covariation.tsv      every cluster's counts, most covariation first
#   passed_covariation.txt    clusters with at least covariation_min_percent
#                             (info.sh) of their base pairs covarying -> step 7
#   failed_06c.tsv            clusters with no result, and why
#   step_06_status.tsv        one row per cluster (see status_06.sh)
#
# Usage:
#   sharcnet/submit.sh 06c DATA_DIR [--force] [--list FILE] [--optional-outputs]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force] [--list FILE] [--optional-outputs]" >&2
    exit 1
fi

data_dir=$1
script_dir=$2
shift 2

force=no
list_file=
optional_outputs=no

while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            force=yes
            ;;
        --optional-outputs)
            optional_outputs=yes
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
            echo "Unknown argument: $1 (expected --force, --list FILE or --optional-outputs)" >&2
            exit 1
            ;;
    esac
    shift
done

echo " step_06c.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

rscape_worker="$script_dir/worker_rscape.sh"
html_worker="$script_dir/worker_html.sh"
status_script="$script_dir/status_06.sh"
round_dir="$step_06_dir/round_$step06_round"

[[ -d "$round_dir" ]] || {
    echo "Missing round $step06_round of step 6: $round_dir -- run step 06a first." >&2
    exit 1
}

for script in "$rscape_worker" "$html_worker" "$status_script"; do
    [[ -x "$script" ]] || {
        echo "Not executable: $script" >&2
        exit 1
    }
done

case $step06_min_rank in
    High | Mid) ;;
    *)
        echo "step06_min_rank must be High or Mid, got '$step06_min_rank'" >&2
        exit 1
        ;;
esac

module load "$apptainer_module"

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

threads=${SLURM_CPUS_PER_TASK:-1}

settings=$(
    printf 'options=%s seed=%s optional_outputs=%s' \
        "$rscape_options" "$rscape_seed" "$optional_outputs"
)

# ---------------------------------------------------------------------------
# which clusters to test
# ---------------------------------------------------------------------------
#
# From the 06b markers rather than from hits_evaluation.tsv, for the same
# reason step 05b reads the 05a markers: they are the record of a finished
# cluster, and nobody edits them by hand.

declare -A rank=()

while IFS=$'\t' read -r cluster field value; do
    [[ "$field" == rank ]] && rank[$cluster]=$value
done < <(
    find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_06b.txt' -print0 |
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
            [[ "$step06_min_rank" == Mid ]] && candidates+=("$cluster")
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

# Nothing to test is not a reason to stop early: the cleanup below still has
# to clear results left by clusters that no longer qualify.
if ((${#run_set[@]} == 0)); then
    echo "No hits alignment ranked $step06_min_rank or better in $round_dir -- nothing to test."
fi

# Candidates and listed clusters together. The failure list and the covariation
# table are rebuilt from this, so a run with --list reports on the clusters it
# tested without dropping what a full run had already recorded.
mapfile -t considered < <(
    printf '%s\n' "${candidates[@]}" "${run_set[@]}" | grep -v '^$' | LC_ALL=C sort -u
)

# ---------------------------------------------------------------------------
# what still needs doing
# ---------------------------------------------------------------------------

alignments=()
missing=()

for cluster in "${run_set[@]}"; do
    cleaned="$round_dir/$cluster/${cluster}_hits_cleaned.sto"

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
    find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_06c.txt' -print0 |
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
    printf 'Step 06c, round %s: %s cluster(s) listed, %s already done, %s to test.\n' \
        "$step06_round" "${#run_set[@]}" "$unchanged" "$todo_total"
    printf '  List:     %s (%s ranked %s or better in all)\n' \
        "$list_file" "$candidate_total" "$step06_min_rank"
else
    printf 'Step 06c, round %s: %s cluster(s) ranked %s or better, %s already done, %s to test.\n' \
        "$step06_round" "$candidate_total" "$step06_min_rank" "$unchanged" "$todo_total"
fi

printf '  Workers:  %s\n' "$threads"
printf '  R-scape:  %s (seed %s, timeout %s)\n' \
    "$rscape_options" "$rscape_seed" "$rscape_timeout"
printf '\n'

if ((${#missing[@]} > 0)); then
    printf '  %s cluster(s) have no _hits_cleaned.sto and are skipped.\n\n' \
        "${#missing[@]}"
fi

# ---------------------------------------------------------------------------
# run R-scape, then the optional HTML views
# ---------------------------------------------------------------------------

if ((todo_total > 0)); then
    for cluster in "${todo[@]}"; do
        printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" \
            "$round_dir/$cluster" \
            "${cluster}_hits_cleaned.sto" \
            "${cluster}_hits_rscape" \
            "${cluster}_hits_rscape" \
            "$rscape_seed" \
            "$rscape_timeout" \
            "$rscape_options"
    done |
        xargs -0 -r -P "$threads" -n8 "$rscape_worker" || true

    for cluster in "${todo[@]}"; do
        cluster_dir="$round_dir/$cluster"
        rm -f "$cluster_dir/${cluster}_hits_cleaned.html"

        [[ "$optional_outputs" == yes ]] || continue
        [[ -s "$cluster_dir/${cluster}_hits_rscape_covariation.tsv" ]] || continue

        printf '%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" "$cluster_dir" \
            "${cluster}_hits_cleaned.sto" "${cluster}_hits_cleaned.html" \
            "$riboswitch_scripts"

        cacofold="${cluster}_hits_rscape/${cluster}_hits_rscape.cacofold.sto"

        if [[ -s "$cluster_dir/$cacofold" ]]; then
            printf '%s\0%s\0%s\0%s\0%s\0' \
                "$rnatools" "$cluster_dir" \
                "$cacofold" "${cacofold%.sto}.html" \
                "$riboswitch_scripts"
        fi
    done |
        xargs -0 -r -P "$threads" -n5 "$html_worker" || true
fi

# ---------------------------------------------------------------------------
# markers and failures
# ---------------------------------------------------------------------------
#
# As in step 05b: the covariation file is the last thing worker_rscape.sh
# writes, so a cluster that has one is done. The failure list covers every
# cluster considered, so a run with --list keeps what a full run recorded; a
# candidate that was never tested shows up as pending in step_06_status.tsv.

failed="$round_dir/failed_06c.tsv"
finished_at=$(date '+%Y-%m-%d %H:%M:%S')

for cluster in "${todo[@]}"; do
    cluster_dir="$round_dir/$cluster"
    covariation="$cluster_dir/${cluster}_hits_rscape_covariation.tsv"
    marker="$cluster_dir/${cluster}_result_06c.txt"

    if [[ ! -s "$covariation" ]]; then
        rm -f "$marker"
        continue
    fi

    read -r cluster_bpairs cluster_covarying cluster_percent < <(
        awk -F '\t' '{ printf "%s %s %s\n", $2, $5, $6 }' "$covariation"
    )

    {
        printf 'cluster\t%s\n' "$cluster"
        printf 'round\t%s\n' "$step06_round"
        printf 'input\t%s\n' "$cluster_dir/${cluster}_hits_cleaned.sto"
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
    cluster_dir="$round_dir/$cluster"

    [[ -s "$cluster_dir/${cluster}_result_06c.txt" ]] && continue

    if [[ ! -s "$cluster_dir/${cluster}_hits_cleaned.sto" ]]; then
        printf '%s\t%s\n' \
            "$cluster" \
            "no _hits_cleaned.sto, nothing for R-scape to test" \
            >> "$failed.tmp"
    elif [[ -n "${attempted[$cluster]:-}" || -f "$cluster_dir/${cluster}_hits_rscape.log" ]]; then
        printf '%s\t%s\n' \
            "$cluster" \
            "R-scape failed, see ${cluster}_hits_rscape.log" \
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
# A cluster whose search was redone in 06a and rescored in 06b can come back
# ranked too low; its old R-scape output describes an alignment that no longer
# exists. Only without --list, as in step 05b.

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

        rm -rf "$cluster_dir/${cluster}_hits_rscape"
        rm -f \
            "$marker" \
            "$cluster_dir/${cluster}_hits_rscape.log" \
            "$cluster_dir/${cluster}_hits_rscape_covariation.tsv" \
            "$cluster_dir/${cluster}_hits_cleaned.html"

        stale=$((stale + 1))
    done < <(
        find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_06c.txt' -print0
    )
fi

# ---------------------------------------------------------------------------
# collated tables
# ---------------------------------------------------------------------------

covariation_table="$round_dir/hits_covariation.tsv"
passed_covariation="$round_dir/passed_covariation.txt"

{
    printf 'cluster\tbpairs\texpected\texpected_sd\tcovarying\tpercent_covarying\n'

    for cluster in "${considered[@]}"; do
        covariation="$round_dir/$cluster/${cluster}_hits_rscape_covariation.tsv"
        [[ -s "$round_dir/$cluster/${cluster}_result_06c.txt" ]] || continue
        [[ -s "$covariation" ]] || continue
        printf '%s\0' "$covariation"
    done |
        xargs -0 -r cat |
        awk -F '\t' 'BEGIN { OFS = "\t" } { sub(/_hits_cleaned\.sto$/, "", $1); print }' |
        LC_ALL=C sort -t $'\t' -k6,6gr -k1,1
} > "$covariation_table.tmp"

mv "$covariation_table.tmp" "$covariation_table"

# The hand-off to step 7, rebuilt every run: copy it before editing it by hand.
awk -F '\t' -v cutoff="$covariation_min_percent" '
    NR > 1 && $6 != "" && $6 >= cutoff {
        print $1
    }
' "$covariation_table" |
    LC_ALL=C sort \
    > "$passed_covariation.tmp"

mv "$passed_covariation.tmp" "$passed_covariation"

tested=$(($(wc -l < "$covariation_table") - 1))
passed_total=$(wc -l < "$passed_covariation")

with_covariation=$(awk -F '\t' 'NR > 1 && $5 > 0 { n++ } END { print n + 0 }' "$covariation_table")

"$status_script" "$data_dir" > /dev/null ||
    echo "status_06.sh failed, step_06_status.tsv may be out of date" >&2

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

printf '\n'
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "hits ranked $step06_min_rank or better" "$candidate_total"
printf '%-42s %10s\n' "already done, unchanged"             "$unchanged"
printf '%-42s %10s\n' "tested this run"                     "$todo_total"
printf '%-42s %10s\n' "with R-scape results"                "$tested"
printf '%-42s %10s\n' "  any base pair covarying"           "$with_covariation"
printf '%-42s %10s\n' "  at least $covariation_min_percent% covarying" "$passed_total"
printf '%-42s %10s\n' "failed"                              "$failed_total"

if ((stale > 0)); then
    printf '%-42s %10s\n' "removed, no longer ranked high enough" "$stale"
fi

printf '\n'
printf '  Covariation: %s\n' "$covariation_table"
printf '  Step 7 list: %s (%s clusters)\n' "$passed_covariation" "$passed_total"
printf '  Status:      %s\n' "$round_dir/step_06_status.tsv"

if ((failed_total > 0)); then
    printf '  Failures:    %s (%s clusters)\n' "$failed" "$failed_total"
fi

printf '\n'
date

# Slurm should mark the job failed when clusters did.
((failed_total == 0))
