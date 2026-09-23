#!/bin/bash
#SBATCH --job-name=rnac_07
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Check that each motif's sequences really fit its structure: realign the motif
# around one representative sequence with HMMER, then see whether the
# structure still covaries on that alignment. Works on step 6 round
# step07_round (info.sh) and writes DATA_DIR/07_compat/round_<n>/<cluster>.
#
# Why: LocARNA aligns sequences and structure together, and can pull a
# sequence into a structured column it does not really belong in, which shows
# up as covariation that is not there. An HMM of a single sequence knows
# nothing about structure, so its alignment is a sequence-only second opinion.
#
# Per cluster:
#   1. representative  one sequence RNA-SCoRE passed in step 6, picked by
#                      upstream's findRepresentative.sh rule (below), or by
#                      hand through --reps
#   2. worker_hmmer.sh HMM of the representative, hmmsearch of all the
#                      motif's sequences with it, and the representative's
#                      structure placed on the resulting alignment
#                      -> <cluster>_hmmer.sto, <cluster>_hmmer_withss.sto
#   3. R-scape twice   --cacofold alone on <cluster>_hmmer.sto: a structure
#                      from covariation only        -> <cluster>_hmmer_cacofold/
#                      rscape_options on _withss.sto: does the carried-over
#                      structure still covary       -> <cluster>_hmmer_withss_rscape/
#   4. comparison      both structures sit on the same alignment columns, so
#                      the base pairs they share are counted directly
#
# Which clusters: those whose step 6c result reached covariation_min_percent
# (06_homolog/round_<n>/passed_covariation.txt). --list FILE runs a hand-picked
# set instead, one cluster per line.
#
# Representative rule (findRepresentative.sh, with its column numbers updated
# for the current RNA-SCoRE output and only Passed sequences considered):
# rank the sequences by base pairs formed, keep the top 30, and among those
# with every stem at 99% or more take the one with the best stem percentages,
# first stem first; if none is at 99% everywhere, the one with most base pairs.
# Ties go to the first name in sort order. --reps FILE sets any cluster's
# representative by hand, one "<cluster>\t<sequence name>" per line.
#
# A cluster is done when it has a <cluster>_result_07.txt marker recording the
# md5 of its input alignment, the representative and the settings; a rerun
# redoes only what changed. --force redoes everything.
#
# For the round as a whole, in DATA_DIR/07_compat/round_<n>:
#   compatibility.tsv       per cluster: sequences kept by hmmsearch, the
#                           covariation of both structures, the base pairs
#                           they share, and whether covariation was kept
#   representatives.tsv     the representative of each cluster, and whether it
#                           was picked by the rule or by hand
#   next_round_seeds.tsv    "<cluster>\t<alignment>" for every cluster whose
#                           carried-over structure still reaches
#                           covariation_min_percent, seeded with that
#                           alignment (<cluster>_hmmer_withss.sto). Edit it --
#                           swap in <cluster>_hmmer_cacofold.sto, the same
#                           alignment with CaCoFold's structure, where that
#                           looks better; drop what should not go on -- and
#                           pass it to step 6 as the next round:
#                             step06_round=<n+1>, submit.sh 06a DATA_DIR --seeds FILE
#   failed_07.tsv           clusters with no result, and why
#
# Usage:
#   sharcnet/submit.sh 07 DATA_DIR [--force] [--list FILE] [--reps FILE]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force] [--list FILE] [--reps FILE]" >&2
    exit 1
fi

data_dir=$1
script_dir=$2
shift 2

force=no
list_file=
reps_file=

while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            force=yes
            ;;
        --list | --reps)
            [[ $# -ge 2 ]] || {
                echo "$1 needs a file" >&2
                exit 1
            }
            if [[ $1 == --list ]]; then list_file=$2; else reps_file=$2; fi
            shift
            ;;
        *)
            echo "Unknown argument: $1 (expected --force, --list FILE or --reps FILE)" >&2
            exit 1
            ;;
    esac
    shift
done

echo " step_07.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

hmmer_worker="$script_dir/worker_hmmer.sh"
rscape_worker="$script_dir/worker_rscape.sh"
addstruct_pl="$step_07_scripts/addStructToHMMERaln.pl"
in_round="$step_06_dir/round_$step07_round"
out_round="$step_07_dir/round_$step07_round"

[[ -d "$in_round" ]] || {
    echo "Missing step 6 round $step07_round: $in_round -- run step 06 first." >&2
    exit 1
}

for script in "$hmmer_worker" "$rscape_worker"; do
    [[ -x "$script" ]] || {
        echo "Not executable: $script" >&2
        exit 1
    }
done

[[ -s "$addstruct_pl" ]] || {
    echo "addStructToHMMERaln.pl not found: $addstruct_pl" >&2
    exit 1
}

case $step07_input in
    cleaned | cacofold) ;;
    *)
        echo "step07_input must be cleaned or cacofold, got '$step07_input'" >&2
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
    printf 'input=%s hmmsearch=%s rscape=%s seed=%s' \
        "$step07_input" "$step07_hmmsearch_options" "$rscape_options" "$rscape_seed"
)

# R-scape writes its alignments with / | : % ' , ; and . in sequence names
# turned into _, so NC_000932.1/5906-6060_-1 comes back as
# NC_000932_1_5906-6060_-1. Its alignments are used here as an input
# (step07_input=cacofold) and handed on as a possible next-round seed, and
# both need the real names: the representative is looked up by name, and step
# 6 reads the accession and coordinates out of them. restore_names maps every
# name back through the alignment R-scape was given.
restore_names() {  # restore_names <alignment R-scape read> <R-scape's alignment> <output>
    awk '
        function mangled(name) {
            gsub(/[\/|:%'"'"',;.]/, "_", name)
            return name
        }

        FNR == NR {
            if ($0 !~ /^#/ && $0 !~ /^\/\// && NF >= 2) original[mangled($1)] = $1
            next
        }

        $0 !~ /^#/ && $0 !~ /^\/\// && NF >= 2 && ($1 in original) {
            $0 = original[$1] substr($0, length($1) + 1)
        }

        ($1 == "#=GS" || $1 == "#=GR") && ($2 in original) {
            at = index($0, $2)
            $0 = substr($0, 1, at - 1) original[$2] substr($0, at + length($2))
        }

        { print }
    ' "$1" "$2" > "$3.tmp"

    mv "$3.tmp" "$3"
}

input_of() {  # input_of <cluster> -> the step 6 alignment step 7 starts from
    if [[ "$step07_input" == cleaned ]]; then
        printf '%s' "$in_round/$1/${1}_hits_cleaned.sto"
    else
        printf '%s' "$in_round/$1/${1}_hits_rscape/${1}_hits_rscape.cacofold.sto"
    fi
}

read_markers() {  # read_markers <directory> <marker glob> <field>...
    local dir=$1 glob=$2
    shift 2

    find "$dir" -mindepth 2 -maxdepth 2 -type f -name "$glob" -print0 |
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

# ---------------------------------------------------------------------------
# which clusters
# ---------------------------------------------------------------------------

clusters=()

if [[ -n "$list_file" ]]; then
    [[ -s "$list_file" ]] || {
        echo "Cluster list not found or empty: $list_file" >&2
        exit 1
    }

    mapfile -t clusters < <(
        grep -vE '^[[:space:]]*(#.*)?$' "$list_file" | awk '{ print $1 }' | LC_ALL=C sort -u
    )
else
    while IFS=$'\t' read -r cluster _ percent; do
        awk -v p="$percent" -v c="$covariation_min_percent" \
            'BEGIN { exit !(p != "" && p + 0 >= c + 0) }' || continue
        clusters+=("$cluster")
    done < <(read_markers "$in_round" '*_result_06c.txt' percent_covarying)

    if ((${#clusters[@]} > 0)); then
        mapfile -t clusters < <(printf '%s\n' "${clusters[@]}" | LC_ALL=C sort)
    fi
fi

cluster_total=${#clusters[@]}

if ((cluster_total == 0)); then
    echo \
        "No cluster in $in_round reached $covariation_min_percent% covariation" \
        "-- nothing to check."
    exit 0
fi

mkdir -p "$out_round"

# ---------------------------------------------------------------------------
# representatives
# ---------------------------------------------------------------------------

declare -A manual_rep=()

if [[ -n "$reps_file" ]]; then
    [[ -s "$reps_file" ]] || {
        echo "Representative list not found or empty: $reps_file" >&2
        exit 1
    }

    while IFS=$'\t' read -r cluster name _; do
        [[ -z "$cluster" || "$cluster" == \#* || -z "${name:-}" ]] && continue
        manual_rep[$cluster]=$name
    done < "$reps_file"
fi

pick_representative() {  # pick_representative <evaluated.tsv>
    awk -F '\t' '
        NR > 1 && $12 == "Passed" {
            print $10 "\t" $1 "\t" $9
        }
    ' "$1" |
        LC_ALL=C sort -t $'\t' -k1,1gr -k2,2 |
        head -n 30 |
        awk -F '\t' '
            # Stem by stem, first difference decides, as compare_values() did.
            function better(a, b,    x, y, n, i) {
                n = split(a, x, ",")
                split(b, y, ",")
                for (i = 1; i <= n; i++) {
                    if (x[i] == "") continue
                    if (x[i] + 0 > y[i] + 0) return 1
                    if (x[i] + 0 < y[i] + 0) return 0
                }
                return 0
            }

            {
                if (NR == 1) top = $2

                n = split($3, stem, ",")
                all = 1

                for (i = 1; i <= n; i++) {
                    if (stem[i] != "" && stem[i] + 0 < 99) {
                        all = 0
                        break
                    }
                }

                if (all && (best == "" || better($3, best_stems))) {
                    best = $2
                    best_stems = $3
                }
            }

            END {
                print (best != "" ? best : top)
            }
        '
}

declare -A rep=()
declare -A rep_source=()
declare -A input=()
missing=()

for cluster in "${clusters[@]}"; do
    path=$(input_of "$cluster")

    if [[ ! -s "$path" ]]; then
        missing+=("$cluster")
        continue
    fi

    input[$cluster]=$path

    if [[ -n "${manual_rep[$cluster]:-}" ]]; then
        rep[$cluster]=${manual_rep[$cluster]}
        rep_source[$cluster]=manual
    elif [[ -s "$in_round/$cluster/${cluster}_hits_evaluated.tsv" ]]; then
        rep[$cluster]=$(pick_representative "$in_round/$cluster/${cluster}_hits_evaluated.tsv")
        rep_source[$cluster]=rule
    fi
done

# ---------------------------------------------------------------------------
# what still needs doing
# ---------------------------------------------------------------------------

declare -A md5_of=()
declare -A input_md5=()

if ((${#input[@]} > 0)); then
    while read -r sum path; do
        md5_of[$path]=$sum
    done < <(printf '%s\0' "${input[@]}" | xargs -0 -r md5sum)

    for cluster in "${!input[@]}"; do
        input_md5[$cluster]=${md5_of[${input[$cluster]}]:-}
    done
fi

declare -A marker_md5=()
declare -A marker_rep=()
declare -A marker_settings=()

while IFS=$'\t' read -r cluster field value; do
    case $field in
        input_md5)      marker_md5[$cluster]=$value ;;
        representative) marker_rep[$cluster]=$value ;;
        settings)       marker_settings[$cluster]=$value ;;
    esac
done < <(read_markers "$out_round" '*_result_07.txt' input_md5 representative settings)

is_current() {  # is_current <cluster>
    [[ "${marker_md5[$1]:-}" == "${input_md5[$1]:-x}" &&
        "${marker_rep[$1]:-}" == "${rep[$1]:-x}" &&
        "${marker_settings[$1]:-}" == "$settings" ]]
}

todo=()
unchanged=0

for cluster in "${clusters[@]}"; do
    [[ -n "${input[$cluster]:-}" && -n "${rep[$cluster]:-}" ]] || continue

    if [[ "$force" != yes ]] && is_current "$cluster"; then
        unchanged=$((unchanged + 1))
        continue
    fi

    todo+=("$cluster")
done

todo_total=${#todo[@]}

printf 'Step 07, round %s: %s cluster(s), %s already done, %s to check.\n' \
    "$step07_round" "$cluster_total" "$unchanged" "$todo_total"
printf '  Input:    %s alignment from %s\n' "$step07_input" "$in_round"
printf '  Workers:  %s\n' "$threads"
printf '  Settings: %s\n' "$settings"
printf '\n'

# ---------------------------------------------------------------------------
# realign, then R-scape on both structures
# ---------------------------------------------------------------------------

if ((todo_total > 0)); then
    for cluster in "${todo[@]}"; do
        dir="$out_round/$cluster"
        mkdir -p "$dir"

        # Start clean: nothing from an earlier representative or input may be
        # read as this run's result.
        rm -rf "$dir/${cluster}_hmmer_cacofold" "$dir/${cluster}_hmmer_withss_rscape"
        rm -f "$dir/${cluster}_result_07.txt" "$dir/${cluster}"_hmmer_*covariation.tsv

        if [[ "$step07_input" == cacofold ]]; then
            restore_names \
                "$in_round/$cluster/${cluster}_hits_cleaned.sto" \
                "${input[$cluster]}" \
                "$dir/${cluster}_input.sto"
        else
            cp "${input[$cluster]}" "$dir/${cluster}_input.sto"
        fi
    done

    for cluster in "${todo[@]}"; do
        printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" \
            "$out_round/$cluster" \
            "${cluster}_input.sto" \
            "${rep[$cluster]}" \
            "${cluster}_hmmer" \
            "$hmmer_bin" \
            "$easel_bin" \
            "$addstruct_pl" \
            "$step07_hmmsearch_options"
    done |
        xargs -0 -r -P "$threads" -n9 "$hmmer_worker" || true

    for cluster in "${todo[@]}"; do
        dir="$out_round/$cluster"
        [[ -s "$dir/${cluster}_hmmer_withss.sto" ]] || continue

        printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" "$dir" "${cluster}_hmmer.sto" \
            "${cluster}_hmmer_cacofold" "${cluster}_hmmer_cacofold" \
            "$rscape_seed" "$rscape_timeout" "--cacofold"

        printf '%s\0%s\0%s\0%s\0%s\0%s\0%s\0%s\0' \
            "$rnatools" "$dir" "${cluster}_hmmer_withss.sto" \
            "${cluster}_hmmer_withss_rscape" "${cluster}_hmmer_withss_rscape" \
            "$rscape_seed" "$rscape_timeout" "$rscape_options"
    done |
        xargs -0 -r -P "$threads" -n8 "$rscape_worker" || true

    # CaCoFold's alignment with the real names back: the alternative seed to
    # swap into next_round_seeds.tsv where that structure looks better.
    for cluster in "${todo[@]}"; do
        dir="$out_round/$cluster"
        folded="$dir/${cluster}_hmmer_cacofold/${cluster}_hmmer_cacofold.cacofold.sto"
        rm -f "$dir/${cluster}_hmmer_cacofold.sto"
        [[ -s "$folded" ]] || continue
        restore_names "$dir/${cluster}_hmmer.sto" "$folded" "$dir/${cluster}_hmmer_cacofold.sto"
    done
fi

# ---------------------------------------------------------------------------
# per-cluster numbers
# ---------------------------------------------------------------------------

# power_counts <.power file> -> "bpairs covarying percent", or "- - -"
power_counts() {
    [[ -s "$1" ]] || {
        echo "- - -"
        return
    }

    awk '
        /^# BPAIRS [0-9]+$/ { bpairs = $NF }
        /^# BPAIRS observed to covary [0-9]+$/ { covarying = $NF }

        END {
            percent = bpairs > 0 ? 100 * covarying / bpairs : 0
            printf "%d %d %.2f\n", bpairs, covarying, percent
        }
    ' "$1"
}

# shared_pairs <withss.sto> <cacofold.sto> -> "shared", or "-" when the two do
# not have the same columns. WUSS: <([{ open and >)]} close; a pseudoknot is an
# upper-case letter closed by its lower-case partner.
shared_pairs() {
    awk '
        function pairs(ss, set,    i, c, key, n) {
            delete depth
            for (i = 1; i <= length(ss); i++) {
                c = substr(ss, i, 1)

                if (index("<([{", c)) { key = "b"; opens = 1 }
                else if (index(">)]}", c)) { key = "b"; opens = 0 }
                else if (c ~ /[A-Z]/) { key = c; opens = 1 }
                else if (c ~ /[a-z]/) { key = toupper(c); opens = 0 }
                else continue

                if (opens) {
                    stack[key, ++depth[key]] = i
                } else if (depth[key] > 0) {
                    set[stack[key, depth[key]--] "," i] = 1
                    n++
                }
            }
            return n
        }

        FNR == 1 { file++ }

        $1 == "#=GC" && $2 == "SS_cons" {
            ss[file] = ss[file] $3
        }

        END {
            if (length(ss[1]) == 0 || length(ss[1]) != length(ss[2])) {
                print "-"
                exit
            }

            pairs(ss[1], given)
            pairs(ss[2], folded)

            for (p in given) if (p in folded) shared++
            print shared + 0
        }
    ' "$1" "$2"
}

failed="$out_round/failed_07.tsv"
finished_at=$(date '+%Y-%m-%d %H:%M:%S')

: > "$failed.tmp"

for cluster in "${missing[@]}"; do
    printf '%s\t%s\n' "$cluster" "no $step07_input alignment in $in_round/$cluster" >> "$failed.tmp"
done

for cluster in "${clusters[@]}"; do
    [[ -n "${input[$cluster]:-}" ]] || continue

    dir="$out_round/$cluster"
    marker="$dir/${cluster}_result_07.txt"

    if [[ -z "${rep[$cluster]:-}" ]]; then
        printf '%s\t%s\n' "$cluster" "no representative: no ${cluster}_hits_evaluated.tsv" >> "$failed.tmp"
        continue
    fi

    if [[ "$force" != yes ]] && is_current "$cluster" && [[ -s "$marker" ]]; then
        continue
    fi

    if [[ ! -s "$dir/${cluster}_hmmer_withss.sto" ]]; then
        printf '%s\t%s\n' "$cluster" "realignment failed, see ${cluster}_hmmer.log" >> "$failed.tmp"
        rm -f "$marker"
        continue
    fi

    if [[ ! -s "$dir/${cluster}_hmmer_cacofold_covariation.tsv" ||
        ! -s "$dir/${cluster}_hmmer_withss_rscape_covariation.tsv" ]]
    then
        printf '%s\t%s\n' "$cluster" "R-scape failed, see ${cluster}_hmmer_*.log" >> "$failed.tmp"
        rm -f "$marker"
        continue
    fi

    seqs_in=$(grep -c '^>' "$dir/${cluster}_hmmer_seqs.fa" || true)

    # hmmsearch names aligned rows <name>/<from>-<to>, and can align two
    # domains of one sequence, so count the names without that suffix.
    seqs_kept=$(
        awk '$0 !~ /^#/ && $0 !~ /^\/\// && NF >= 2 { sub(/\/[0-9]+-[0-9]+$/, "", $1); print $1 }' \
            "$dir/${cluster}_hmmer.sto" | LC_ALL=C sort -u | wc -l
    )

    read -r given_bpairs given_covarying given_percent < <(
        power_counts "$dir/${cluster}_hmmer_withss_rscape/${cluster}_hmmer_withss_rscape.power"
    )

    read -r folded_bpairs folded_covarying folded_percent < <(
        power_counts "$dir/${cluster}_hmmer_cacofold/${cluster}_hmmer_cacofold.cacofold.power"
    )

    shared=$(
        shared_pairs \
            "$dir/${cluster}_hmmer_withss.sto" \
            "$dir/${cluster}_hmmer_cacofold/${cluster}_hmmer_cacofold.cacofold.sto" \
            2> /dev/null || echo "-"
    )

    keeps=no
    awk -v p="$given_percent" -v c="$covariation_min_percent" \
        'BEGIN { exit !(p != "-" && p + 0 >= c + 0) }' && keeps=yes

    {
        printf 'cluster\t%s\n' "$cluster"
        printf 'round\t%s\n' "$step07_round"
        printf 'input\t%s\n' "${input[$cluster]}"
        printf 'input_md5\t%s\n' "${input_md5[$cluster]}"
        printf 'representative\t%s\n' "${rep[$cluster]}"
        printf 'representative_from\t%s\n' "${rep_source[$cluster]}"
        printf 'settings\t%s\n' "$settings"
        printf 'seqs_in\t%s\n' "$seqs_in"
        printf 'seqs_kept\t%s\n' "$seqs_kept"
        printf 'given_bpairs\t%s\n' "$given_bpairs"
        printf 'given_covarying\t%s\n' "$given_covarying"
        printf 'given_percent\t%s\n' "$given_percent"
        printf 'cacofold_bpairs\t%s\n' "$folded_bpairs"
        printf 'cacofold_covarying\t%s\n' "$folded_covarying"
        printf 'cacofold_percent\t%s\n' "$folded_percent"
        printf 'shared_bpairs\t%s\n' "$shared"
        printf 'keeps_covariation\t%s\n' "$keeps"
        printf 'finished\t%s\n' "$finished_at"
    } > "$marker.tmp"

    mv "$marker.tmp" "$marker"
done

LC_ALL=C sort "$failed.tmp" > "$failed"
rm -f "$failed.tmp"

failed_total=$(wc -l < "$failed")

# ---------------------------------------------------------------------------
# collated tables
# ---------------------------------------------------------------------------
#
# Over the clusters of this run's selection only, so a cluster step 6 no longer
# passes (or that a --list left out) is not reported, though its folder stays.

compatibility="$out_round/compatibility.tsv"
representatives="$out_round/representatives.tsv"
next_seeds="$out_round/next_round_seeds.tsv"

declare -A value=()

while IFS=$'\t' read -r cluster field content; do
    value[$cluster,$field]=$content
done < <(
    read_markers "$out_round" '*_result_07.txt' \
        representative representative_from seqs_in seqs_kept \
        given_bpairs given_covarying given_percent \
        cacofold_bpairs cacofold_covarying cacofold_percent \
        shared_bpairs keeps_covariation
)

done_total=0
kept_total=0

{
    printf 'cluster\tseqs_in\tseqs_kept\tgiven_bpairs\tgiven_covarying\tgiven_percent'
    printf '\tcacofold_bpairs\tcacofold_covarying\tcacofold_percent\tshared_bpairs'
    printf '\tkeeps_covariation\n'

    for cluster in "${clusters[@]}"; do
        [[ -n "${value[$cluster,seqs_in]:-}" ]] || continue
        done_total=$((done_total + 1))
        [[ "${value[$cluster,keeps_covariation]}" == yes ]] && kept_total=$((kept_total + 1))

        printf '%s' "$cluster"
        for field in seqs_in seqs_kept given_bpairs given_covarying given_percent \
            cacofold_bpairs cacofold_covarying cacofold_percent shared_bpairs keeps_covariation
        do
            printf '\t%s' "${value[$cluster,$field]}"
        done
        printf '\n'
    done
} > "$compatibility.tmp"

mv "$compatibility.tmp" "$compatibility"

{
    printf 'cluster\trepresentative\tpicked_by\n'

    for cluster in "${clusters[@]}"; do
        [[ -n "${value[$cluster,representative]:-}" ]] || continue
        printf '%s\t%s\t%s\n' \
            "$cluster" "${value[$cluster,representative]}" "${value[$cluster,representative_from]}"
    done
} > "$representatives.tmp"

mv "$representatives.tmp" "$representatives"

# The hand-off to step 6's next round, rebuilt every run: copy it before
# editing it by hand.
for cluster in "${clusters[@]}"; do
    [[ "${value[$cluster,keeps_covariation]:-}" == yes ]] || continue
    printf '%s\t%s\n' "$cluster" "$out_round/$cluster/${cluster}_hmmer_withss.sto"
done > "$next_seeds.tmp"

mv "$next_seeds.tmp" "$next_seeds"

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

printf '\n'
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "selected from step 6"                "$cluster_total"
printf '%-42s %10s\n' "already done, unchanged"             "$unchanged"
printf '%-42s %10s\n' "checked this run"                    "$todo_total"
printf '%-42s %10s\n' "with a result"                       "$done_total"
printf '%-42s %10s\n' "  covariation kept (>= $covariation_min_percent%)" "$kept_total"
printf '%-42s %10s\n' "failed"                              "$failed_total"
printf '\n'
printf '  Compatibility:   %s\n' "$compatibility"
printf '  Representatives: %s\n' "$representatives"
printf '  Next round:      %s (%s clusters)\n' "$next_seeds" "$kept_total"

if ((failed_total > 0)); then
    printf '  Failures:        %s (%s clusters)\n' "$failed" "$failed_total"
fi

printf '\n'
date

((failed_total == 0))
