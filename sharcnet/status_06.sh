#!/bin/bash

# Where every cluster of the current step 6 round stands: the search, the hits
# RNA-SCoRE kept, whether the seed was found again, and the covariation.
#
# Everything comes from the per-cluster markers, the 06a failure notes and the
# 06b/06c failure lists, so this can be run at any time, including while the
# 06a array is still going. step_06b.sh and step_06c.sh both refresh it when
# they finish.
#
# Writes DATA_DIR/06_homolog/round_<n>/step_06_status.tsv:
#   cluster            cluster name
#   step_06a           done, failed, or pending (queued or running)
#   hits               homologs cmsearch reported
#   step_06b           done, failed, pending, or - before 06a is done
#   rank               RNA-SCoRE on the hits: High, Mid, Low, or none (no hits)
#   passed             hits that passed RNA-SCoRE
#   seed_seqs          sequences in the seed alignment
#   seeds_recovered    of those, how many a hit landed on
#   top_hit_is_seed    whether the best hit is one of them
#   new_hits           hits outside the seed
#   step_06c           done, failed, pending, or - when the rank is too low
#   bpairs             base pairs in the hits alignment's structure
#   covarying          base pairs R-scape found covarying
#   percent_covarying  covarying as a percentage of bpairs
#   passes             yes when at or above covariation_min_percent -> step 7
#
# Usage:
#   sharcnet/status_06.sh DATA_DIR

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 DATA_DIR" >&2
    exit 1
fi

data_dir=$1

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

round_dir="$step_06_dir/round_$step06_round"

[[ -d "$round_dir" ]] || {
    echo "Missing round $step06_round of step 6: $round_dir -- run step 06a first." >&2
    exit 1
}

status="$round_dir/step_06_status.tsv"

mapfile -t clusters < <(
    find "$round_dir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | LC_ALL=C sort
)

# Every marker field, one awk pass over all three sub-steps' markers.
declare -A value=()

while IFS=$'\t' read -r cluster step field content; do
    value[$cluster,$step,$field]=$content
done < <(
    find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_06?.txt' -print0 |
        xargs -0 -r awk -F '\t' '
            FNR == 1 {
                n = split(FILENAME, path, "/")
                cluster = path[n - 1]
                step = path[n]
                sub(/.*_result_/, "", step)
                sub(/\.txt$/, "", step)
            }

            NF == 2 {
                print cluster "\t" step "\t" $1 "\t" $2
            }
        '
)

declare -A failed_06b=()
declare -A failed_06c=()

for step in 06b 06c; do
    [[ -s "$round_dir/failed_$step.tsv" ]] || continue

    while IFS=$'\t' read -r cluster _; do
        if [[ $step == 06b ]]; then failed_06b[$cluster]=yes; else failed_06c[$cluster]=yes; fi
    done < "$round_dir/failed_$step.tsv"
done

declare -A count=()

tally() {  # tally <key>
    count[$1]=$((${count[$1]:-0} + 1))
}

{
    printf 'cluster\tstep_06a\thits\tstep_06b\trank\tpassed\tseed_seqs'
    printf '\tseeds_recovered\ttop_hit_is_seed\tnew_hits\tstep_06c\tbpairs'
    printf '\tcovarying\tpercent_covarying\tpasses\n'

    for cluster in "${clusters[@]}"; do
        if [[ -n "${value[$cluster,06a,hits]:-}" ]]; then
            state_a="done"
        elif [[ -s "$round_dir/$cluster/${cluster}_failed_06a.txt" ]]; then
            state_a="failed"
        else
            state_a="pending"
        fi

        rank=${value[$cluster,06b,rank]:-}

        if [[ -n "$rank" ]]; then
            state_b="done"
        elif [[ "$state_a" != "done" ]]; then
            state_b="-"
        elif [[ -n "${failed_06b[$cluster]:-}" ]]; then
            state_b="failed"
        else
            state_b="pending"
        fi

        candidate=no
        [[ "$rank" == High ]] && candidate=yes
        [[ "$rank" == Mid && "$step06_min_rank" == Mid ]] && candidate=yes

        percent=${value[$cluster,06c,percent_covarying]:-}
        passes="-"

        if [[ -n "${value[$cluster,06c,bpairs]:-}" ]]; then
            state_c="done"
            passes=no
            awk -v p="$percent" -v c="$covariation_min_percent" \
                'BEGIN { exit !(p != "" && p + 0 >= c + 0) }' && passes=yes
        elif [[ "$candidate" == no ]]; then
            state_c="-"
        elif [[ -n "${failed_06c[$cluster]:-}" ]]; then
            state_c="failed"
        else
            state_c="pending"
        fi

        tally "a_$state_a"
        tally "b_$state_b"
        tally "c_$state_c"
        tally "passes_$passes"

        printf '%s' "$cluster"
        printf '\t%s' \
            "$state_a" \
            "${value[$cluster,06a,hits]:--}" \
            "$state_b" \
            "${rank:--}" \
            "${value[$cluster,06b,passed]:--}" \
            "${value[$cluster,06b,seed_seqs]:--}" \
            "${value[$cluster,06b,seeds_recovered]:--}" \
            "${value[$cluster,06b,top_hit_is_seed]:--}" \
            "${value[$cluster,06b,new_hits]:--}" \
            "$state_c" \
            "${value[$cluster,06c,bpairs]:--}" \
            "${value[$cluster,06c,covarying]:--}" \
            "${percent:--}" \
            "$passes"
        printf '\n'
    done
} > "$status.tmp"

mv "$status.tmp" "$status"

printf 'Step 06 status, round %s: %s clusters.\n' "$step06_round" "${#clusters[@]}"
printf '  06a: %s done, %s failed, %s pending\n' \
    "${count[a_done]:-0}" "${count[a_failed]:-0}" "${count[a_pending]:-0}"
printf '  06b: %s done, %s failed, %s pending\n' \
    "${count[b_done]:-0}" "${count[b_failed]:-0}" "${count[b_pending]:-0}"
printf '  06c: %s done, %s failed, %s pending, %s not eligible (hits not ranked %s or better)\n' \
    "${count[c_done]:-0}" "${count[c_failed]:-0}" "${count[c_pending]:-0}" \
    "${count[c_-]:-0}" "$step06_min_rank"
printf '  At or above %s%% covarying: %s\n' "$covariation_min_percent" "${count[passes_yes]:-0}"
printf '  Table: %s\n' "$status"
