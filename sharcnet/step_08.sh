#!/bin/bash
#SBATCH --job-name=rnac_08
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Check the seeds step 7 hands to the next round against Rfam, and drop the
# ones that are a known RNA family. Works on round step08_round (info.sh):
# reads 07_compat/round_<n> and writes 08_rfam/round_<n>, including the seed
# list for step 6 round n+1.
#
# Why: the pipeline finds every conserved structure in the genomes, known ones
# included (tRNAs, rRNA domains, RNase P, tmRNA, ...). Rediscovering them shows
# the pipeline works, and they are reported here, but another homolog search
# round on a tRNA only spends time. Removing them keeps the next round on the
# candidates.
#
# Seeds. By default 07_compat/round_<n>/next_round_seeds.tsv. --seeds FILE
# takes a curated list instead, same format, "<cluster>\t<alignment>" per line
# (blank lines and lines starting with # ignored), so clusters added by hand,
# dropped, or swapped to their CaCoFold alignment come through as edited. A
# relative alignment path is read relative to the seed file's directory.
#
# Per cluster: the alignment's sequences, ungapped, are scanned with cmscan
# against Rfam (rfam_dir in info.sh) using Rfam's own gathering thresholds
# (--cut_ga) and clan competition, so every hit reported is a family Rfam
# itself would assign. A cluster is "known" when at least
# step08_min_fraction of its sequences hit the same family, "partial" when
# some hit but fewer, and "none" when none hit. Only known clusters are
# dropped; partial ones are worth a look by eye (rfam_by_cluster.tsv).
#
# Only the seeds are scanned, not every homolog: a seed already names the
# family, and thousands of hits per cluster would not change that answer.
#
# A cluster is done when it has a <cluster>_result_08.txt marker recording
# the md5 of its seed alignment, of Rfam.cm and the cmscan options, so a rerun
# scans exactly the clusters whose seed, Rfam release or options changed.
# Changing step08_min_fraction needs no rescan: known/partial/none is decided
# when the tables are written. --force rescans everything.
#
# Per cluster, in DATA_DIR/08_rfam/round_<n>/<cluster>:
#   <cluster>_seqs.fa          the seed's sequences, ungapped
#   <cluster>_rfam.tsv         its Rfam hits (one row per sequence and family)
#   <cluster>_result_08.txt    the marker
#
# And for the round as a whole, in DATA_DIR/08_rfam/round_<n>:
#   cmscan.tblout, cmscan.txt  this run's raw cmscan output
#   rfam_hits.tsv              every cluster's Rfam hits under one header
#   rfam_by_cluster.tsv        one row per cluster: top family, how many of
#                              its sequences hit it, status
#   known.txt                  clusters that are a known family
#   novel.txt                  the rest (partial and none)
#   seeds_round_<n+1>.tsv      the input seeds minus known.txt -> step 6
#                                step06_round=<n+1>,
#                                submit.sh 06a DATA_DIR --seeds FILE
#   failed_08.tsv              clusters with no result, and why
#
# rfam_dir must hold Rfam.cm (pressed: Rfam.cm.i1m and friends), Rfam.clanin
# and a VERSION file with "release\t<release>" and "md5\t<md5 of Rfam.cm>"
# lines. Nothing is downloaded here: compute nodes may have no internet.
#
# Usage:
#   sharcnet/submit.sh 08 DATA_DIR [--force] [--seeds FILE]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force] [--seeds FILE]" >&2
    exit 1
fi

data_dir=$1
script_dir=$2
shift 2

force=no
seeds_file=

while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            force=yes
            ;;
        --seeds)
            [[ $# -ge 2 ]] || {
                echo "$1 needs a file" >&2
                exit 1
            }
            seeds_file=$2
            shift
            ;;
        *)
            echo "Unknown argument: $1 (expected --force or --seeds FILE)" >&2
            exit 1
            ;;
    esac
    shift
done

echo " step_08.sh: data_dir = $data_dir, script_dir = $script_dir"

source "$script_dir/info.sh"

in_round="$step_07_dir/round_$step08_round"
out_round="$step_08_dir/round_$step08_round"
next_round=$((step08_round + 1))

[[ -n "$seeds_file" ]] || seeds_file="$in_round/next_round_seeds.tsv"

[[ -s "$seeds_file" ]] || {
    echo "Seed list not found or empty: $seeds_file -- run step 07 first, or pass --seeds FILE." >&2
    exit 1
}

seeds_file=$(realpath "$seeds_file")
seeds_dir=$(dirname "$seeds_file")

for file in Rfam.cm Rfam.cm.i1m Rfam.clanin VERSION; do
    [[ -s "$rfam_dir/$file" ]] || {
        echo "Rfam is not set up: $rfam_dir/$file missing (rfam_dir in info.sh)." >&2
        exit 1
    }
done

rfam_md5=$(awk -F '\t' '$1 == "md5" { print $2 }' "$rfam_dir/VERSION")
rfam_release=$(awk -F '\t' '$1 == "release" { print $2 }' "$rfam_dir/VERSION")

[[ -n "$rfam_md5" ]] || {
    echo "No md5 line in $rfam_dir/VERSION." >&2
    exit 1
}

awk -v f="$step08_min_fraction" 'BEGIN { exit !(f + 0 > 0 && f + 0 <= 1) }' || {
    echo "step08_min_fraction must be in (0, 1], got '$step08_min_fraction'" >&2
    exit 1
}

module load "$apptainer_module"

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

threads=${SLURM_CPUS_PER_TASK:-1}
settings="rfam_md5=$rfam_md5 cmscan=$step08_cmscan_options"

mkdir -p "$out_round"

# ---------------------------------------------------------------------------
# which clusters, and their seed alignments
# ---------------------------------------------------------------------------

clusters=()
declare -A seed=()
failures=()

while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^[[:space:]]*(#.*)?$ ]] && continue

    IFS=$'\t' read -r cluster path _ <<< "$line"

    if [[ -z "${path:-}" ]]; then
        failures+=("$cluster"$'\t'"no alignment path in $seeds_file")
        continue
    fi

    [[ "$path" == /* ]] || path="$seeds_dir/$path"

    if [[ -n "${seed[$cluster]:-}" ]]; then
        echo "Cluster listed twice in $seeds_file: $cluster" >&2
        exit 1
    fi

    clusters+=("$cluster")
    seed[$cluster]=$path
done < "$seeds_file"

cluster_total=${#clusters[@]}

if ((cluster_total == 0)); then
    echo "No clusters in $seeds_file -- nothing to check."
    exit 0
fi

declare -A seed_md5=()
missing=()

for cluster in "${clusters[@]}"; do
    if [[ ! -s "${seed[$cluster]}" ]]; then
        missing+=("$cluster")
        continue
    fi
    seed_md5[$cluster]=$(md5sum "${seed[$cluster]}" | awk '{ print $1 }')
done

# ---------------------------------------------------------------------------
# what still needs scanning
# ---------------------------------------------------------------------------

marker_field() {  # marker_field <marker> <field>
    awk -F '\t' -v f="$2" '$1 == f { print $2; exit }' "$1" 2> /dev/null || true
}

is_current() {  # is_current <cluster>
    local marker="$out_round/$1/${1}_result_08.txt"
    [[ -s "$marker" &&
        "$(marker_field "$marker" seed_md5)" == "${seed_md5[$1]:-x}" &&
        "$(marker_field "$marker" settings)" == "$settings" ]]
}

todo=()
unchanged=0

for cluster in "${clusters[@]}"; do
    [[ -n "${seed_md5[$cluster]:-}" ]] || continue

    if [[ "$force" != yes ]] && is_current "$cluster"; then
        unchanged=$((unchanged + 1))
        continue
    fi

    todo+=("$cluster")
done

todo_total=${#todo[@]}

printf 'Step 08, round %s: %s cluster(s), %s already done, %s to scan.\n' \
    "$step08_round" "$cluster_total" "$unchanged" "$todo_total"
printf '  Seeds:    %s\n' "$seeds_file"
printf '  Rfam:     %s (%s), md5 %s\n' "$rfam_dir" "$rfam_release" "$rfam_md5"
printf '  Workers:  %s\n' "$threads"
printf '  Settings: %s\n' "$settings"
printf '\n'

# ---------------------------------------------------------------------------
# sequences, then one cmscan over all of them
# ---------------------------------------------------------------------------
#
# Each sequence is named <cluster>|<name> in the combined FASTA, so the hits
# split back to their cluster afterwards. Stockholm can wrap a sequence over
# several blocks; the pieces are joined by name. Gap characters (. - _ ~) go.

extract_sequences() {  # extract_sequences <cluster> <alignment>
    awk -v cluster="$1" '
        /^#/ || /^\/\// || NF < 2 { next }
        {
            if (!($1 in seq)) order[++n] = $1
            seq[$1] = seq[$1] $2
        }
        END {
            for (i = 1; i <= n; i++) {
                s = seq[order[i]]
                gsub(/[-._~]/, "", s)
                if (s == "") continue
                print ">" cluster "|" order[i]
                print toupper(s)
            }
        }
    ' "$2"
}

scan_fasta="$out_round/cmscan_input.fa"
scan_tbl="$out_round/cmscan.tblout"
scan_txt="$out_round/cmscan.txt"
scan_ok=yes
declare -A seqs_of=()

if ((todo_total > 0)); then
    : > "$scan_fasta"

    for cluster in "${todo[@]}"; do
        dir="$out_round/$cluster"
        mkdir -p "$dir"
        rm -f "$dir/${cluster}_result_08.txt" "$dir/${cluster}_rfam.tsv"

        extract_sequences "$cluster" "${seed[$cluster]}" > "$dir/${cluster}_seqs.fa"
        seqs_of[$cluster]=$(grep -c '^>' "$dir/${cluster}_seqs.fa" || true)
        cat "$dir/${cluster}_seqs.fa" >> "$scan_fasta"
    done

    rm -f "$scan_tbl" "$scan_txt"

    if [[ -s "$scan_fasta" ]]; then
        # shellcheck disable=SC2086
        if ! apptainer exec \
            --bind "$out_round:/work" \
            --bind "$rfam_dir:/rfam:ro" \
            --pwd /work \
            "$rnatools" \
            cmscan \
            --cpu "$threads" \
            $step08_cmscan_options \
            --fmt 2 \
            --clanin /rfam/Rfam.clanin \
            --oclan \
            --tblout "$(basename "$scan_tbl")" \
            -o "$(basename "$scan_txt")" \
            /rfam/Rfam.cm \
            "$(basename "$scan_fasta")"
        then
            echo "cmscan failed, see $scan_txt" >&2
            scan_ok=no
        fi
    fi
fi

# ---------------------------------------------------------------------------
# per-cluster hits and markers
# ---------------------------------------------------------------------------
#
# cmscan --fmt 2 tblout columns: 1 idx, 2 family, 3 accession, 4 query,
# 6 clan, 10-11 seq from/to, 12 strand, 17 score, 18 E-value, 19 inc, 20 olp.
# Kept: significant hits (inc "!") that are not a lower-scoring overlap of a
# better one (olp "="). Per sequence and family only the best hit counts.

failed="$out_round/failed_08.tsv"
finished_at=$(date '+%Y-%m-%d %H:%M:%S')

: > "$failed.tmp"

for entry in "${failures[@]}"; do
    printf '%s\n' "$entry" >> "$failed.tmp"
done

for cluster in "${missing[@]}"; do
    printf '%s\t%s\n' "$cluster" "seed alignment not found: ${seed[$cluster]}" >> "$failed.tmp"
done

if ((todo_total > 0)) && [[ "$scan_ok" == yes ]]; then
    hits_tmp="$out_round/hits.tmp"

    if [[ -s "$scan_tbl" ]]; then
        awk '
            /^#/ { next }
            $19 != "!" || $20 == "=" { next }
            {
                split($4, q, "|")
                cluster = q[1]
                name = substr($4, length(cluster) + 2)
                key = cluster SUBSEP name SUBSEP $3
                if (key in best && best[key] + 0 >= $17 + 0) next
                best[key] = $17
                row[key] = cluster "\t" name "\t" $2 "\t" $3 "\t" $6 "\t" \
                    $10 "\t" $11 "\t" $12 "\t" $17 "\t" $18
            }
            END { for (k in row) print row[k] }
        ' "$scan_tbl" | LC_ALL=C sort -t $'\t' -k1,1 -k2,2 -k9,9gr > "$hits_tmp"
    else
        : > "$hits_tmp"
    fi

    for cluster in "${todo[@]}"; do
        dir="$out_round/$cluster"
        marker="$dir/${cluster}_result_08.txt"
        seqs=${seqs_of[$cluster]:-0}

        if ((seqs == 0)); then
            printf '%s\t%s\n' "$cluster" "no sequences in ${seed[$cluster]}" >> "$failed.tmp"
            continue
        fi

        awk -F '\t' -v c="$cluster" '$1 == c' "$hits_tmp" > "$dir/${cluster}_rfam.tsv"

        # family counts: sequences per family, best score and E-value
        read -r seqs_hit top_family top_accession top_clan top_seqs top_score top_evalue others < <(
            awk -F '\t' '
                { hit[$2] = 1 }
                {
                    fam = $4
                    name[fam] = $3; clan[fam] = $5; count[fam]++
                    if (!(fam in score) || $9 + 0 > score[fam] + 0) { score[fam] = $9; evalue[fam] = $10 }
                }
                END {
                    for (s in hit) n++
                    for (f in count) {
                        if (top == "" || count[f] > count[top] ||
                            (count[f] == count[top] && score[f] + 0 > score[top] + 0)) top = f
                    }
                    if (top == "") { print n + 0, "-", "-", "-", 0, "-", "-", "-"; exit }
                    for (f in count) if (f != top) list = list (list ? "," : "") name[f] ":" count[f]
                    print n + 0, name[top], top, clan[top], count[top], score[top], evalue[top], (list ? list : "-")
                }
            ' "$dir/${cluster}_rfam.tsv"
        )

        {
            printf 'cluster\t%s\n' "$cluster"
            printf 'round\t%s\n' "$step08_round"
            printf 'seed\t%s\n' "${seed[$cluster]}"
            printf 'seed_md5\t%s\n' "${seed_md5[$cluster]}"
            printf 'rfam_release\t%s\n' "$rfam_release"
            printf 'settings\t%s\n' "$settings"
            printf 'seqs\t%s\n' "$seqs"
            printf 'seqs_with_hit\t%s\n' "$seqs_hit"
            printf 'top_family\t%s\n' "$top_family"
            printf 'top_accession\t%s\n' "$top_accession"
            printf 'top_clan\t%s\n' "$top_clan"
            printf 'top_family_seqs\t%s\n' "$top_seqs"
            printf 'top_score\t%s\n' "$top_score"
            printf 'top_evalue\t%s\n' "$top_evalue"
            printf 'other_families\t%s\n' "$others"
            printf 'finished\t%s\n' "$finished_at"
        } > "$marker.tmp"

        mv "$marker.tmp" "$marker"
    done

    rm -f "$hits_tmp" "$scan_fasta"
elif ((todo_total > 0)); then
    for cluster in "${todo[@]}"; do
        printf '%s\t%s\n' "$cluster" "cmscan failed, see $scan_txt" >> "$failed.tmp"
    done
fi

LC_ALL=C sort -u "$failed.tmp" > "$failed"
rm -f "$failed.tmp"

failed_total=$(wc -l < "$failed")

# ---------------------------------------------------------------------------
# collated tables
# ---------------------------------------------------------------------------
#
# Over this run's seed list only, so a cluster dropped from the list is not
# reported, though its folder stays.

hits_table="$out_round/rfam_hits.tsv"
by_cluster="$out_round/rfam_by_cluster.tsv"
known_list="$out_round/known.txt"
novel_list="$out_round/novel.txt"
next_seeds="$out_round/seeds_round_${next_round}.tsv"

declare -A value=()

for cluster in "${clusters[@]}"; do
    marker="$out_round/$cluster/${cluster}_result_08.txt"
    [[ -s "$marker" ]] || continue
    # a marker from an older seed or Rfam is not this run's result
    is_current "$cluster" || continue

    while IFS=$'\t' read -r field content; do
        value[$cluster,$field]=$content
    done < "$marker"
done

done_total=0
known_total=0
partial_total=0
none_total=0

{
    printf 'cluster\tsequence\tfamily\taccession\tclan\tseq_from\tseq_to\tstrand\tscore\tevalue\n'
    for cluster in "${clusters[@]}"; do
        [[ -n "${value[$cluster,seqs]:-}" ]] || continue
        cat "$out_round/$cluster/${cluster}_rfam.tsv"
    done
} > "$hits_table.tmp"

mv "$hits_table.tmp" "$hits_table"

: > "$known_list.tmp"
: > "$novel_list.tmp"

{
    printf 'cluster\tseqs\tseqs_with_hit\ttop_family\ttop_accession\ttop_clan'
    printf '\ttop_family_seqs\ttop_family_fraction\ttop_score\ttop_evalue\tother_families\tstatus\n'

    for cluster in "${clusters[@]}"; do
        [[ -n "${value[$cluster,seqs]:-}" ]] || continue
        done_total=$((done_total + 1))

        fraction=$(
            awk -v a="${value[$cluster,top_family_seqs]}" -v b="${value[$cluster,seqs]}" \
                'BEGIN { printf "%.2f", (b > 0 ? a / b : 0) }'
        )

        if awk -v f="$fraction" -v m="$step08_min_fraction" 'BEGIN { exit !(f + 0 >= m + 0) }'; then
            status=known
            known_total=$((known_total + 1))
            printf '%s\n' "$cluster" >> "$known_list.tmp"
        elif ((${value[$cluster,seqs_with_hit]} > 0)); then
            status=partial
            partial_total=$((partial_total + 1))
            printf '%s\n' "$cluster" >> "$novel_list.tmp"
        else
            status=none
            none_total=$((none_total + 1))
            printf '%s\n' "$cluster" >> "$novel_list.tmp"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$cluster" \
            "${value[$cluster,seqs]}" "${value[$cluster,seqs_with_hit]}" \
            "${value[$cluster,top_family]}" "${value[$cluster,top_accession]}" \
            "${value[$cluster,top_clan]}" "${value[$cluster,top_family_seqs]}" \
            "$fraction" "${value[$cluster,top_score]}" "${value[$cluster,top_evalue]}" \
            "${value[$cluster,other_families]}" "$status"
    done
} > "$by_cluster.tmp"

mv "$by_cluster.tmp" "$by_cluster"
mv "$known_list.tmp" "$known_list"
mv "$novel_list.tmp" "$novel_list"

# The hand-off to step 6's next round: the seed lines as given, minus the
# known families and anything that failed here. Rebuilt every run: copy it
# before editing it by hand.
for cluster in "${clusters[@]}"; do
    [[ -n "${value[$cluster,seqs]:-}" ]] || continue
    grep -qxF "$cluster" "$known_list" && continue
    printf '%s\t%s\n' "$cluster" "${seed[$cluster]}"
done > "$next_seeds.tmp"

mv "$next_seeds.tmp" "$next_seeds"

next_total=$(wc -l < "$next_seeds")

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

printf '\n'
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "in the seed list"                    "$cluster_total"
printf '%-42s %10s\n' "already done, unchanged"             "$unchanged"
printf '%-42s %10s\n' "scanned this run"                    "$todo_total"
printf '%-42s %10s\n' "with a result"                       "$done_total"
printf '%-42s %10s\n' "  known family (>= $step08_min_fraction of seqs)" "$known_total"
printf '%-42s %10s\n' "  partial Rfam hits"                 "$partial_total"
printf '%-42s %10s\n' "  no Rfam hit"                       "$none_total"
printf '%-42s %10s\n' "failed"                              "$failed_total"
printf '\n'
printf '  By cluster:  %s\n' "$by_cluster"
printf '  All hits:    %s\n' "$hits_table"
printf '  Known:       %s\n' "$known_list"
printf '  Novel:       %s\n' "$novel_list"
printf '  Next round:  %s (%s clusters)\n' "$next_seeds" "$next_total"

if ((failed_total > 0)); then
    printf '  Failures:    %s (%s clusters)\n' "$failed" "$failed_total"
fi

if ((known_total > 0)); then
    printf '\n  Known families:\n'
    awk -F '\t' 'NR > 1 && $12 == "known" { printf "    %-40s %s (%s, %s of %s seqs)\n", $1, $4, $5, $7, $2 }' "$by_cluster"
fi

printf '\n'
date

((failed_total == 0))
