#!/bin/bash
#SBATCH --job-name=rnac_06b
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Align the homologs step 06a found to their covariance model, score the hits
# alignment with RNA-SCoRE, and check whether the search found the seed again.
# Works on every cluster of the current round (step06_round in info.sh) that
# has a finished 06a search, in DATA_DIR/06_homolog/round_<n>/<cluster>.
#
# Per cluster:
#   1. cmsearch_reformatv1_1.pl + cmalign   -> <cluster>_hits.fna, _hits.sto
#   2. RNA-SCoRE with the step06_rnascore_* thresholds
#                                           -> <cluster>_hits_cleaned.sto,
#                                              _hits_evaluated.tsv, _hits_rank.tsv
#   3. seed recovery: which of the seed's sequences a hit overlaps (same genome,
#      same strand, inside the step 1 window the seed sequence came from),
#      whether the best hit is one of them, and how many hits are new. This is
#      the check the step 6 README asks for by hand: a model that does not find
#      its own seed is not a reliable one.
#
# A cluster is done when it has a <cluster>_result_06b.txt marker. The marker
# records what 06a found (the md5 of the seed and of the hits table) and the
# RNA-SCoRE settings, so a rerun redoes exactly the clusters whose search or
# settings changed. --force redoes everything. A search with no hits at all is
# recorded as done, ranked "none".
#
# --optional-outputs also renders <cluster>_hits.html (stockholm_to_html.pl).
#
# And for the round as a whole, in DATA_DIR/06_homolog/round_<n>:
#   hits_evaluation.tsv     every rank line under one header
#   seed_recovery.tsv       hits, seed sequences found again, best hit, new hits
#   high_mid_hits.txt       clusters whose hits ranked High or Mid -> step 6c
#   failed_06b.tsv          clusters with no result, and why
#   step_06_status.tsv      one row per cluster (see status_06.sh)
#
# Usage:
#   sharcnet/submit.sh 06b DATA_DIR [--force] [--optional-outputs]

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

echo " step_06b.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

score_pl="$step_05_scripts/RNA-SCoRE.pl"
align_worker="$script_dir/worker_cmalign.sh"
score_worker="$script_dir/worker_rnascore.sh"
html_worker="$script_dir/worker_html.sh"
status_script="$script_dir/status_06.sh"
round_dir="$step_06_dir/round_$step06_round"

[[ -d "$round_dir" ]] || {
    echo "Missing round $step06_round of step 6: $round_dir -- run step 06a first." >&2
    exit 1
}

[[ -s "$score_pl" ]] || {
    echo "RNA-SCoRE.pl not found: $score_pl" >&2
    exit 1
}

for script in "$align_worker" "$score_worker" "$html_worker" "$status_script"; do
    [[ -x "$script" ]] || {
        echo "Not executable: $script" >&2
        exit 1
    }
done

module load "$apptainer_module"

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

threads=${SLURM_CPUS_PER_TASK:-1}

settings=$(
    printf 'mt=%s t=%s gc=%s d=%s seed=%s optional_outputs=%s' \
        "$step06_rnascore_mt" \
        "$step06_rnascore_bp" \
        "$step06_rnascore_gc" \
        "$step06_rnascore_dupl" \
        "$perl_hash_seed" \
        "$optional_outputs"
)

# ---------------------------------------------------------------------------
# what 06a found, and what 06b already did with it
# ---------------------------------------------------------------------------

read_markers() {  # read_markers <marker glob> <field>...
    local glob=$1
    shift

    find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name "$glob" -print0 |
        xargs -0 -r awk -F '\t' -v fields="$*" '
            BEGIN {
                n = split(fields, list, " ")
                for (i = 1; i <= n; i++) wanted[list[i]] = 1
            }

            FNR == 1 {
                split(FILENAME, path, "/")
                cluster = path[length(path) - 1]
            }

            $1 in wanted {
                print cluster "\t" $1 "\t" $2
            }
        '
}

declare -A seed_md5=()
declare -A hits_md5=()
declare -A hits=()

while IFS=$'\t' read -r cluster field value; do
    case $field in
        seed_md5) seed_md5[$cluster]=$value ;;
        hits_md5) hits_md5[$cluster]=$value ;;
        hits)     hits[$cluster]=$value ;;
    esac
done < <(read_markers '*_result_06a.txt' seed_md5 hits_md5 hits)

# What 06b's result depends on from 06a: the seed (for the recovery check) and
# the hits themselves.
declare -A found=()

for cluster in "${!seed_md5[@]}"; do
    found[$cluster]="seed=${seed_md5[$cluster]} hits=${hits_md5[$cluster]:-}"
done

mapfile -t clusters < <(printf '%s\n' "${!found[@]}" | grep -v '^$' | LC_ALL=C sort)

cluster_total=${#clusters[@]}

if ((cluster_total == 0)); then
    echo "No finished 06a search in $round_dir -- nothing to do." >&2
    exit 1
fi

declare -A marker_input=()
declare -A marker_settings=()

while IFS=$'\t' read -r cluster field value; do
    case $field in
        input)    marker_input[$cluster]=$value ;;
        settings) marker_settings[$cluster]=$value ;;
    esac
done < <(read_markers '*_result_06b.txt' input settings)

is_current() {  # is_current <cluster>
    [[ "${marker_input[$1]:-}" == "${found[$1]}" &&
        "${marker_settings[$1]:-}" == "$settings" ]]
}

todo=()
unchanged=0

for cluster in "${clusters[@]}"; do
    if [[ "$force" != yes ]] && is_current "$cluster"; then
        unchanged=$((unchanged + 1))
        continue
    fi

    todo+=("$cluster")
done

todo_total=${#todo[@]}

printf 'Step 06b, round %s: %s searched cluster(s), %s already done, %s to process.\n' \
    "$step06_round" "$cluster_total" "$unchanged" "$todo_total"
printf '  Workers:    %s\n' "$threads"
printf '  Thresholds: %s\n' "$settings"
printf '\n'

# ---------------------------------------------------------------------------
# align, then score what aligned
# ---------------------------------------------------------------------------

if ((todo_total > 0)); then
    # A cluster being redone starts without its previous RNA-SCoRE results, so
    # a failed realignment cannot leave an old rank line to be recorded.
    for cluster in "${todo[@]}"; do
        base="$round_dir/$cluster/${cluster}_hits"
        rm -f "${base}_rank.tsv" "${base}_cleaned.sto" "${base}_evaluated.tsv"
    done

    for cluster in "${todo[@]}"; do
        ((${hits[$cluster]:-0} > 0)) || continue

        printf '%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" \
            "$round_dir/$cluster" \
            "${cluster}.cm" \
            "${cluster}_hits.txt" \
            "${cluster}_hits" \
            "$riboswitch_scripts"
    done |
        xargs -0 -r -P "$threads" -n6 "$align_worker" || true

    for cluster in "${todo[@]}"; do
        [[ -s "$round_dir/$cluster/${cluster}_hits.sto" ]] || continue

        printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$score_pl" \
            "$round_dir/$cluster" \
            "${cluster}_hits.sto" \
            "_hits.sto" \
            "$step06_rnascore_mt" \
            "$step06_rnascore_bp" \
            "$step06_rnascore_gc" \
            "$step06_rnascore_dupl" \
            "$perl_hash_seed"
    done |
        xargs -0 -r -P "$threads" -n9 "$score_worker" || true

    for cluster in "${todo[@]}"; do
        rm -f "$round_dir/$cluster/${cluster}_hits.html"

        [[ "$optional_outputs" == yes ]] || continue
        [[ -s "$round_dir/$cluster/${cluster}_hits.sto" ]] || continue

        printf '%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" \
            "$round_dir/$cluster" \
            "${cluster}_hits.sto" \
            "${cluster}_hits.html" \
            "$riboswitch_scripts"
    done |
        xargs -0 -r -P "$threads" -n5 "$html_worker" || true
fi

# ---------------------------------------------------------------------------
# seed recovery
# ---------------------------------------------------------------------------
#
# Seed sequences name the step 1 window they came from:
#   <accession>_<counter>_<start>_<end>[r]-<motif start>_<motif end>[r]
# with start 0-based and end exclusive, r for the reverse complement. Seeds
# from a later round are named like hits, <accession>/<from>-<to>_<strand>,
# and when step 7's hmmsearch realigned them they carry one more
# /<from>-<to> (the part of that hit it aligned), which is ignored. A seed
# counts as recovered when a hit on the same accession and strand overlaps
# that window or hit: the motif coordinates trimAlignment.pl appends are
# counted in alignment columns, not genome positions, so the window is the
# only exact location a step 5 name carries.
#
# Prints: seed sequences, recovered, best hit is a seed (yes/no/-), new hits.

seed_recovery() {  # seed_recovery <seed.sto> <hits.tbl>
    awk '
        function locate(id,    rest, w, a, b, t) {
            # the /<from>-<to> hmmsearch adds after a hit name
            if (id ~ /\/[0-9]+-[0-9]+_-?1\/[0-9]+-[0-9]+$/) sub(/\/[0-9]+-[0-9]+$/, "", id)

            if (match(id, /_[0-9]+_[0-9]+_[0-9]+r?(-[0-9]+_[0-9]+r?)?$/)) {
                where_acc = substr(id, 1, RSTART - 1)
                rest = substr(id, RSTART + 1)
                sub(/-.*/, "", rest)
                split(rest, w, "_")
                where_strand = (w[3] ~ /r$/) ? "-" : "+"
                sub(/r$/, "", w[3])
                where_lo = w[2] + 1
                where_hi = w[3] + 0
                return 1
            }

            if (match(id, /\/[0-9]+-[0-9]+(_-?1)?$/)) {
                where_acc = substr(id, 1, RSTART - 1)
                rest = substr(id, RSTART + 1)
                where_strand = (rest ~ /_-1$/) ? "-" : "+"
                sub(/_-?1$/, "", rest)
                split(rest, w, "-")
                a = w[1] + 0
                b = w[2] + 0
                if (a > b) { t = a; a = b; b = t; where_strand = "-" }
                where_lo = a
                where_hi = b
                return 1
            }

            return 0
        }

        FNR == NR {
            if ($0 ~ /^#/ || $0 ~ /^\// || NF < 2) next
            seeds++
            if (!locate($1)) next
            placed++
            acc[placed] = where_acc
            strand[placed] = where_strand
            lo[placed] = where_lo
            hi[placed] = where_hi
            next
        }

        /^#/ {
            next
        }

        {
            rank++
            a = $8 + 0
            b = $9 + 0
            if (a > b) { t = a; a = b; b = t }

            is_seed = 0

            for (s = 1; s <= placed; s++) {
                if (acc[s] == $1 && strand[s] == $10 && lo[s] <= b && hi[s] >= a) {
                    recovered[s] = 1
                    is_seed = 1
                }
            }

            if (rank == 1) top = is_seed ? "yes" : "no"
            if (!is_seed) new++
        }

        END {
            for (s in recovered) found++
            printf "%d %d %s %d\n", seeds, found, (rank ? top : "-"), new
        }
    ' "$1" "$2"
}

# ---------------------------------------------------------------------------
# markers and failures
# ---------------------------------------------------------------------------
#
# As in step 5, what happened is read back off disk: a cluster with a rank line
# (or a search that found nothing) is done and gets its marker; one without
# gets a line in failed_06b.tsv, and any marker from an earlier run is removed.

failed="$round_dir/failed_06b.tsv"
finished_at=$(date '+%Y-%m-%d %H:%M:%S')

: > "$failed.tmp"

for cluster in "${clusters[@]}"; do
    cluster_dir="$round_dir/$cluster"
    marker="$cluster_dir/${cluster}_result_06b.txt"
    rank_file="$cluster_dir/${cluster}_hits_rank.tsv"

    if [[ "$force" != yes ]] && is_current "$cluster" && [[ -s "$marker" ]]; then
        continue
    fi

    rank=
    passed=0

    if ((${hits[$cluster]:-0} == 0)); then
        rank=none
    elif [[ -s "$rank_file" ]]; then
        read -r rank passed < <(
            awk -F '\t' '
                {
                    split($NF, confidence, " ")
                    printf "%s %s\n", confidence[1], $7
                }
            ' "$rank_file"
        )
    fi

    if [[ -z "$rank" ]]; then
        if [[ ! -s "$cluster_dir/${cluster}_hits.sto" ]]; then
            reason="cmsearch_reformat or cmalign failed, see ${cluster}_hits_cmalign.log"
        else
            reason="RNA-SCoRE failed"
        fi

        rm -f "$marker"
        printf '%s\t%s\n' "$cluster" "$reason" >> "$failed.tmp"
        continue
    fi

    read -r seed_seqs recovered top_hit new_hits < <(
        seed_recovery "$cluster_dir/${cluster}_seed.sto" "$cluster_dir/${cluster}_hits.tbl"
    )

    {
        printf 'cluster\t%s\n' "$cluster"
        printf 'round\t%s\n' "$step06_round"
        printf 'input\t%s\n' "${found[$cluster]}"
        printf 'settings\t%s\n' "$settings"
        printf 'hits\t%s\n' "${hits[$cluster]:-0}"
        printf 'rank\t%s\n' "$rank"
        printf 'passed\t%s\n' "$passed"
        printf 'seed_seqs\t%s\n' "$seed_seqs"
        printf 'seeds_recovered\t%s\n' "$recovered"
        printf 'top_hit_is_seed\t%s\n' "$top_hit"
        printf 'new_hits\t%s\n' "$new_hits"
        printf 'finished\t%s\n' "$finished_at"
    } > "$marker.tmp"

    mv "$marker.tmp" "$marker"
done

LC_ALL=C sort "$failed.tmp" > "$failed"
rm -f "$failed.tmp"

failed_total=$(wc -l < "$failed")

# A 06b marker whose 06a search is no longer done (it was redone and failed, or
# is queued again) describes hits that are being replaced: clear it, and step
# 06c then clears the R-scape result built on it.
retired=0

while IFS= read -r -d '' marker; do
    cluster=${marker%/*}
    cluster=${cluster##*/}

    [[ -n "${found[$cluster]:-}" ]] && continue

    rm -f "$marker"
    retired=$((retired + 1))
done < <(
    find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_06b.txt' -print0
)

# ---------------------------------------------------------------------------
# collated tables
# ---------------------------------------------------------------------------

evaluation="$round_dir/hits_evaluation.tsv"
recovery="$round_dir/seed_recovery.tsv"
high_mid="$round_dir/high_mid_hits.txt"

{
    printf 'clusterFile\tNseqs\tNuniqSeqs\tssConsensus\tss_consLen'
    printf '\tTotal_basepairs\tNseqsPassingEval\tConfidenceOnStructure\n'

    for cluster in "${clusters[@]}"; do
        rank_file="$round_dir/$cluster/${cluster}_hits_rank.tsv"
        [[ -s "$round_dir/$cluster/${cluster}_result_06b.txt" ]] || continue
        [[ -s "$rank_file" ]] || continue
        printf '%s\0' "$rank_file"
    done |
        xargs -0 -r cat |
        LC_ALL=C sort
} > "$evaluation.tmp"

mv "$evaluation.tmp" "$evaluation"

{
    printf 'cluster\thits\trank\tseed_seqs\tseeds_recovered\tpercent_recovered'
    printf '\ttop_hit_is_seed\tnew_hits\n'

    read_markers '*_result_06b.txt' \
        hits rank seed_seqs seeds_recovered top_hit_is_seed new_hits |
        awk -F '\t' '
            {
                value[$1, $2] = $3
                seen[$1] = 1
            }

            END {
                for (c in seen) {
                    seeds = value[c, "seed_seqs"] + 0
                    found = value[c, "seeds_recovered"] + 0

                    percent = seeds > 0 ? 100 * found / seeds : 0

                    printf "%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\n", \
                        c, value[c, "hits"], value[c, "rank"], \
                        seeds, found, percent, \
                        value[c, "top_hit_is_seed"], value[c, "new_hits"]
                }
            }
        ' |
        LC_ALL=C sort
} > "$recovery.tmp"

mv "$recovery.tmp" "$recovery"

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

read -r no_hits all_seeds top_seed < <(
    awk -F '\t' '
        NR == 1 {
            next
        }

        {
            if ($3 == "none") none++
            if ($4 > 0 && $5 == $4) all++
            if ($7 == "yes") top++
        }

        END {
            printf "%d %d %d\n", none, all, top
        }
    ' "$recovery"
)

"$status_script" "$data_dir" > /dev/null ||
    echo "status_06.sh failed, step_06_status.tsv may be out of date" >&2

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

printf '\n'
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "searched by step 06a"               "$cluster_total"
printf '%-42s %10s\n' "already done, unchanged"            "$unchanged"
printf '%-42s %10s\n' "processed this run"                 "$todo_total"
printf '%-42s %10s\n' "no hits at all"                     "$no_hits"
printf '%-42s %10s\n' "every seed sequence found again"    "$all_seeds"
printf '%-42s %10s\n' "best hit is a seed sequence"        "$top_seed"
printf '%-42s %10s\n' "hits scored by RNA-SCoRE"           "$scored"
printf '%-42s %10s\n' "  ranked High"                      "$high"
printf '%-42s %10s\n' "  ranked Mid"                       "$mid"
printf '%-42s %10s\n' "  ranked Low (dropped)"             "$low"
printf '%-42s %10s\n' "failed"                             "$failed_total"

if ((retired > 0)); then
    printf '%-42s %10s\n' "markers cleared, 06a search not done" "$retired"
fi

printf '\n'
printf '  Evaluation:    %s\n' "$evaluation"
printf '  Seed recovery: %s\n' "$recovery"
printf '  Status:        %s\n' "$round_dir/step_06_status.tsv"
printf '  Step 6c list:  %s (%s clusters)\n' "$high_mid" "$((high + mid))"

if ((failed_total > 0)); then
    printf '  Failures:      %s (%s clusters)\n' "$failed" "$failed_total"
fi

printf '\n'
date

((failed_total == 0))
