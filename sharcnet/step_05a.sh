#!/bin/bash
#SBATCH --job-name=rnac_05a
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Trim every cluster step 4 folded down to its motif region, then score that
# motif with RNA-SCoRE. Outputs are written to DATA_DIR/05_evaluation/<cluster>.
#
# The input is whatever step 4 has produced: every
# DATA_DIR/04_locarna/<cluster>/<cluster>_result.stk there is. Widening the
# cluster size range in step 4 and running this again therefore picks up the
# new clusters without any bookkeeping by hand.
#
# A cluster is done when it has a <cluster>_result_05a.txt marker. The marker
# records the md5 of the result.stk it was built from and the RNA-SCoRE
# settings it was built with, so a rerun redoes exactly the clusters whose
# input or settings changed and leaves the rest alone. --force redoes
# everything.
#
# --optional-outputs also keeps <cluster>_motif.fasta (for RNAz),
# <cluster>_motif.afa (for SQUARNA) and <cluster>_motif.aln (a clustal view,
# via esl-reformat in the container). Nothing in this pipeline reads them, and
# the setting is part of the marker, so turning it on redoes the clusters that
# were done without it.
#
# Per cluster, in DATA_DIR/05_evaluation/<cluster>:
#   <cluster>_motif.sto            the trimmed alignment
#   <cluster>_motif_cleaned.sto    the sequences that passed  -> step 5b, step 6
#   <cluster>_motif_evaluated.tsv  verdict per sequence       -> step 7
#   <cluster>_motif_rank.tsv       this cluster's rank line
#   <cluster>_result_05a.txt       the marker
#
# And for the run as a whole, in DATA_DIR/05_evaluation:
#   allClusters_evaluation.tsv     every rank line under one header
#   high_mid_clusters.txt          clusters ranked High or Mid -> step 5b, step 6
#   failed_05a.tsv                 clusters with no result, and why
#   step_05_status.tsv             one row per cluster (see status_05.sh)
#
# Usage:
#   sharcnet/submit.sh 05a DATA_DIR [--force] [--optional-outputs]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force] [--optional-outputs]" >&2
    exit 1
fi

data_dir=$1
script_dir=$2
shift 2

force=no
optional_outputs=no

while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            force=yes
            ;;
        --optional-outputs)
            optional_outputs=yes
            ;;
        *)
            echo "Unknown argument: $1 (expected --force or --optional-outputs)" >&2
            exit 1
            ;;
    esac
    shift
done

echo " step_05a.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

trim_pl="$step_05_scripts/trimAlignment.pl"
score_pl="$step_05_scripts/RNA-SCoRE.pl"
trim_worker="$script_dir/worker_trim.sh"
score_worker="$script_dir/worker_rnascore.sh"
status_script="$script_dir/status_05.sh"

[[ -d "$step_04_dir" ]] || {
    echo "Missing step 4 output: $step_04_dir" >&2
    exit 1
}

[[ -s "$trim_pl" ]] || {
    echo "trimAlignment.pl not found: $trim_pl" >&2
    exit 1
}

[[ -s "$score_pl" ]] || {
    echo \
        "RNA-SCoRE.pl not found: $score_pl" \
        "-- clone it from RodrigoReisLab/RNA-SCoRE into that directory" \
        >&2
    exit 1
}

for worker in "$trim_worker" "$score_worker" "$status_script"; do
    [[ -x "$worker" ]] || {
        echo "Not executable: $worker" >&2
        exit 1
    }
done

threads=${SLURM_CPUS_PER_TASK:-1}

# The container is only needed for the clustal view.
if [[ "$optional_outputs" == yes ]]; then
    module load "$apptainer_module"

    [[ -s "$rnatools" ]] || {
        echo "Container not found: $rnatools" >&2
        exit 1
    }
fi

mkdir -p "$step_05_dir"

# Recorded in every marker: a cluster whose settings differ from these is
# redone, so changing a threshold in info.sh does not leave old results behind.
settings=$(
    printf 'mt=%s t=%s gc=%s d=%s seed=%s optional_outputs=%s' \
        "$step05_rnascore_mt" \
        "$step05_rnascore_bp" \
        "$step05_rnascore_gc" \
        "$step05_rnascore_dupl" \
        "$perl_hash_seed" \
        "$optional_outputs"
)

# Clusters to consider: the ones step 4 actually folded. The directory name is
# the cluster name, and worker_04.sh names the result after it.
mapfile -t clusters < <(
    find "$step_04_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result.stk' -printf '%P\n' |
        sed 's|/.*||' |
        LC_ALL=C sort -u
)

cluster_total=${#clusters[@]}

if ((cluster_total == 0)); then
    echo "No *_result.stk under $step_04_dir -- run step 04 first." >&2
    exit 1
fi

results=()

for cluster in "${clusters[@]}"; do
    results+=("$step_04_dir/$cluster/${cluster}_result.stk")
done

declare -A input_md5=()

while read -r sum path; do
    cluster=${path%/*}
    cluster=${cluster##*/}
    input_md5[$cluster]=$sum
done < <(printf '%s\0' "${results[@]}" | xargs -0 -r md5sum)

# What the existing markers were built from. One awk pass over all of them
# rather than reading each marker in the loop below.
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
    find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_05a.txt' -print0 |
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

for cluster in "${clusters[@]}"; do
    if [[ "$force" != yes &&
        "${marker_md5[$cluster]:-}" == "${input_md5[$cluster]:-x}" &&
        "${marker_settings[$cluster]:-}" == "$settings" ]]
    then
        unchanged=$((unchanged + 1))
        continue
    fi

    todo+=("$cluster")
done

todo_total=${#todo[@]}

printf 'Step 05a: %s folded cluster(s), %s already done, %s to process.\n' \
    "$cluster_total" "$unchanged" "$todo_total"
printf '  Workers:    %s\n' "$threads"
printf '  Thresholds: %s\n' "$settings"
printf '\n'

# ---------------------------------------------------------------------------
# trim, then score what trimmed
# ---------------------------------------------------------------------------
#
# Two passes rather than one worker doing both: the tools stay one per script,
# so steps 6 and 7 can call worker_rnascore.sh with their own alignments.
# A cluster that failed to trim is dropped before scoring, which keeps
# RNA-SCoRE away from motifs with no sequences in them.

if ((todo_total > 0)); then
    # A cluster being redone starts without its previous RNA-SCoRE results:
    # if trimming fails this time, RNA-SCoRE never runs, and an old rank line
    # left in place would be recorded as this run's result.
    for cluster in "${todo[@]}"; do
        base="$step_05_dir/$cluster/${cluster}_motif"
        rm -f "${base}_rank.tsv" "${base}_cleaned.sto" "${base}_evaluated.tsv"
    done

    for cluster in "${todo[@]}"; do
        printf '%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$trim_pl" \
            "$step_04_dir/$cluster/${cluster}_result.stk" \
            "$step_05_dir/$cluster" \
            "$cluster" \
            "$optional_outputs" \
            "$rnatools"
    done |
        xargs -0 -r -P "$threads" -n6 "$trim_worker" || true

    scorable=()

    for cluster in "${todo[@]}"; do
        [[ -s "$step_05_dir/$cluster/${cluster}_motif.sto" ]] || continue
        scorable+=("$cluster")
    done

    if ((${#scorable[@]} > 0)); then
        for cluster in "${scorable[@]}"; do
            printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
                "$score_pl" \
                "$step_05_dir/$cluster" \
                "${cluster}_motif.sto" \
                "_motif.sto" \
                "$step05_rnascore_mt" \
                "$step05_rnascore_bp" \
                "$step05_rnascore_gc" \
                "$step05_rnascore_dupl" \
                "$perl_hash_seed"
        done |
            xargs -0 -r -P "$threads" -n9 "$score_worker" || true
    fi
fi

# ---------------------------------------------------------------------------
# markers and failures
# ---------------------------------------------------------------------------
#
# The workers report a failure on STDERR (which lands in this job's log) and
# never write to a shared file, so what happened is worked out here instead,
# from what is on disk. A cluster with a rank line is done and gets its marker;
# one without gets a line in failed_05a.tsv, and any marker from an earlier run
# is removed so a stale one cannot pass for a current result.

failed="$step_05_dir/failed_05a.tsv"
finished_at=$(date '+%Y-%m-%d %H:%M:%S')

: > "$failed.tmp"

for cluster in "${clusters[@]}"; do
    cluster_dir="$step_05_dir/$cluster"
    rank_file="$cluster_dir/${cluster}_motif_rank.tsv"
    marker="$cluster_dir/${cluster}_result_05a.txt"

    # Untouched clusters keep the marker they already have.
    if [[ -s "$marker" &&
        "${marker_md5[$cluster]:-}" == "${input_md5[$cluster]:-x}" &&
        "${marker_settings[$cluster]:-}" == "$settings" ]]
    then
        continue
    fi

    # A cluster whose step 4 result went missing between the listing and here
    # has no md5 to record, so whatever it still has on disk is stale.
    if [[ -s "$rank_file" && -n "${input_md5[$cluster]:-}" ]]; then
        read -r rank passed < <(
            awk -F '\t' '
                {
                    split($NF, confidence, " ")
                    printf "%s %s\n", confidence[1], $7
                }
            ' "$rank_file"
        )

        {
            printf 'cluster\t%s\n' "$cluster"
            printf 'input\t%s\n' "$step_04_dir/$cluster/${cluster}_result.stk"
            printf 'input_md5\t%s\n' "${input_md5[$cluster]}"
            printf 'settings\t%s\n' "$settings"
            printf 'rank\t%s\n' "$rank"
            printf 'passed\t%s\n' "$passed"
            printf 'finished\t%s\n' "$finished_at"
        } > "$marker.tmp"

        mv "$marker.tmp" "$marker"
        continue
    fi

    if [[ -z "${input_md5[$cluster]:-}" ]]; then
        reason="no ${cluster}_result.stk in $step_04_dir/$cluster"
    elif [[ -s "$cluster_dir/${cluster}_motif.sto.invalid" ]]; then
        reason="empty motif, trimAlignment.pl wrote no sequence rows"
    elif [[ ! -s "$cluster_dir/${cluster}_motif.sto" ]]; then
        reason="trimAlignment.pl failed"
    else
        reason="RNA-SCoRE failed"
    fi

    rm -f "$marker"
    printf '%s\t%s\n' "$cluster" "$reason" >> "$failed.tmp"
done

LC_ALL=C sort "$failed.tmp" > "$failed"
rm -f "$failed.tmp"

failed_total=$(wc -l < "$failed")

# A marker whose step 4 result is gone describes a cluster that no longer
# exists; without this, step 5b would keep treating it as a candidate. Only the
# marker goes -- the files stay for anyone who wants to look at them.
declare -A folded=()

for cluster in "${clusters[@]}"; do
    folded[$cluster]=yes
done

retired=0

while IFS= read -r -d '' marker; do
    cluster=${marker%/*}
    cluster=${cluster##*/}

    [[ -n "${folded[$cluster]:-}" ]] && continue

    rm -f "$marker"
    retired=$((retired + 1))
done < <(
    find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_05a.txt' -print0
)

# ---------------------------------------------------------------------------
# collated tables
# ---------------------------------------------------------------------------
#
# Rebuilt from the per-cluster files on every run, so they always describe what
# is on disk now, and clusters whose step 4 result has since been removed drop
# out of them. Same reasoning as step 4's timings.tsv.

evaluation="$step_05_dir/allClusters_evaluation.tsv"
high_mid="$step_05_dir/high_mid_clusters.txt"

{
    printf 'clusterFile\tNseqs\tNuniqSeqs\tssConsensus\tss_consLen'
    printf '\tTotal_basepairs\tNseqsPassingEval\tConfidenceOnStructure\n'

    for cluster in "${clusters[@]}"; do
        rank_file="$step_05_dir/$cluster/${cluster}_motif_rank.tsv"
        [[ -s "$step_05_dir/$cluster/${cluster}_result_05a.txt" ]] || continue
        [[ -s "$rank_file" ]] || continue
        printf '%s\0' "$rank_file"
    done |
        xargs -0 -r cat |
        LC_ALL=C sort
} > "$evaluation.tmp"

mv "$evaluation.tmp" "$evaluation"

# The hand-off to step 5b and step 6: cluster names only, High before Mid, so
# it can be copied and cut down by hand and fed back with --list.
awk -F '\t' '
    NR == 1 {
        next
    }

    {
        split($NF, confidence, " ")

        if (confidence[1] == "High") print "1\t" $1
        else if (confidence[1] == "Mid") print "2\t" $1
    }
' "$evaluation" |
    LC_ALL=C sort -t $'\t' -k1,1 -k2,2 |
    cut -f2 \
    > "$high_mid.tmp"

mv "$high_mid.tmp" "$high_mid"

read -r scored high mid low < <(
    awk -F '\t' '
        NR == 1 {
            next
        }

        {
            split($NF, confidence, " ")
            scored++

            if (confidence[1] == "High") high++
            else if (confidence[1] == "Mid") mid++
            else low++
        }

        END {
            printf "%d %d %d %d\n", scored, high, mid, low
        }
    ' "$evaluation"
)

"$status_script" "$data_dir" > /dev/null ||
    echo "status_05.sh failed, step_05_status.tsv may be out of date" >&2

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

printf '\n'
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "folded by step 4"                   "$cluster_total"
printf '%-42s %10s\n' "already done, unchanged"            "$unchanged"
printf '%-42s %10s\n' "processed this run"                 "$todo_total"
printf '%-42s %10s\n' "scored by RNA-SCoRE"                "$scored"
printf '%-42s %10s\n' "  ranked High"                      "$high"
printf '%-42s %10s\n' "  ranked Mid"                       "$mid"
printf '%-42s %10s\n' "  ranked Low (dropped)"             "$low"
printf '%-42s %10s\n' "failed"                             "$failed_total"

if ((retired > 0)); then
    printf '%-42s %10s\n' "markers cleared, step 4 result gone" "$retired"
fi

printf '\n'
printf '  Evaluation:  %s\n' "$evaluation"
printf '  Status:      %s\n' "$step_05_dir/step_05_status.tsv"
printf '  Step 5b list: %s (%s clusters)\n' "$high_mid" "$((high + mid))"

if ((failed_total > 0)); then
    printf '  Failures:    %s (%s clusters)\n' "$failed" "$failed_total"
fi

printf '\n'
date

# Slurm should mark the job failed when clusters did, so the mail says so and
# sacct can be trusted. Everything above has already been written.
((failed_total == 0))
