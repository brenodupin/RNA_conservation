#!/bin/bash

# Where a dataset stands from step 00 through step 03, and what step 4 has
# left to do: genomes, windows, clusters and their sizes, which clusters
# step 3 aligned and screened (and which it did not), how many passed, and how
# the passing clusters fall into the step 04a / 04b size ranges.
#
# Read-only, and fast enough for a login node: every count comes from one
# find or awk pass, never from a test per cluster. Run it before submitting
# step 4, or any time while step 3 or step 4 is going.
#
# A cluster counts as aligned when 03a left a non-empty <cluster>_aligned.aln,
# and as screened when 03b left a non-empty <cluster>_RNALalifold.out: the
# 03b worker creates its output file before RNALalifold runs, so an empty one
# is a failed screen, not a cluster without structure.
#
# Usage:
#   sharcnet/status_03.sh DATA_DIR [--list]
#
# --list also prints every cluster step 3 did not finish, not just the first
# few.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 DATA_DIR [--list]" >&2
    exit 1
fi

data_dir=$1
list_all=${2:-}

if [[ -n "$list_all" && "$list_all" != "--list" ]]; then
    echo "Unknown argument: $list_all (expected --list)" >&2
    exit 1
fi

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

shown=10
[[ -n "$list_all" ]] && shown=1000000000

row() {  # row <label> <value> [note]
    printf '  %-34s %10s  %s\n' "$1" "$2" "${3:-}"
}

count_fasta() {  # count_fasta <file> -> records, or "-" when missing
    if [[ -s "$1" ]]; then grep -c '^>' "$1"; else echo "-"; fi
}

printf 'Steps 00-03 status: %s\n' "$data_dir"

# ---------------------------------------------------------------------------
# step 00: genomes
# ---------------------------------------------------------------------------

printf '\n== step 00: genomes ==\n'

shopt -s nullglob
inputs=("$input_dir"/*.fasta)
onelines=("$step_00_dir"/*_oneLine.fasta)
shopt -u nullglob

row "input FASTA files" "${#inputs[@]}"
row "one-line FASTA files" "${#onelines[@]}"

if ((${#onelines[@]} > 0)); then
    # sequences, total length, shortest and longest
    read -r n_seqs total_nt min_nt max_nt < <(
        awk '
            /^>/ { next }
            {
                n++; len = length($0); total += len
                if (min == "" || len < min) min = len
                if (len > max) max = len
            }
            END { printf "%d %d %d %d\n", n, total, min, max }
        ' "${onelines[@]}"
    )
    row "sequences" "$n_seqs"
    row "total length (nt)" "$total_nt"
    row "shortest / longest (nt)" "$min_nt / $max_nt"

    duplicates=$(
        grep -h '^>' "${onelines[@]}" | awk '{ print $1 }' | LC_ALL=C sort | uniq -d | wc -l
    )
    row "duplicate accessions" "$duplicates" "$( ((duplicates > 0)) && echo '<- should be 0' )"
fi

# ---------------------------------------------------------------------------
# step 01: windows
# ---------------------------------------------------------------------------

printf '\n== step 01: windows (%s nt, overlap %s) ==\n' "$window_size" "$overlap"

stride=$((window_size - overlap))
raw_windows=$(count_fasta "$step_01_dir/${windows_prefix}.fasta")
no_ns=$(count_fasta "$step_01_dir/${windows_prefix}_processed_noNs.fasta")
no_polyn=$(count_fasta "$step_01_dir/${windows_prefix}_processed_noNs_polyN.fasta")
unique=$(count_fasta "$step_01_dir/${windows_prefix}_processed_noNs_polyN_uniq.fasta")

expected=""
if [[ -n "${total_nt:-}" && "$total_nt" -gt 0 ]]; then
    expected="(about $((total_nt / stride)) expected from the genome length)"
fi

row "windows" "$raw_windows" "$expected"
row "without Ns" "$no_ns"
row "without poly-N runs" "$no_polyn"
row "unique (step 2 input)" "$unique"

# ---------------------------------------------------------------------------
# step 02: clusters
# ---------------------------------------------------------------------------

cluster_count="$step_02_dir/cluster_count.tsv"

printf '\n== step 02: clusters (pass 1 %s id, pass 2 %s id, coverage %s) ==\n' \
    "$pid_pass_1" "$pid_pass_2" "$coverage"

if [[ ! -s "$cluster_count" ]]; then
    row "clusters with 2+ members" "-" "<- no $cluster_count yet"
    exit 0
fi

LC_ALL=C sort -t $'\t' -k2,2 "$cluster_count" > "$work/clusters.tsv"

clusters_total=$(wc -l < "$work/clusters.tsv")
splits=$(find "$step_02_dir/splits" -maxdepth 1 -name '*_cluster.fasta' 2> /dev/null | wc -l)

row "clusters with 2+ members" "$clusters_total"
row "cluster FASTA files (splits)" "$splits" \
    "$( ((splits != clusters_total)) && echo '<- should match the line above' )"

# size bins, in the step 4 ranges from info.sh
bin_of='
    function bin(s) {
        if (s < a_min) return "1:2-" (a_min - 1) " (below step 4)"
        if (s <= a_max) return "2:" a_min "-" a_max " (04a)"
        if (s >= b_min && s <= b_max) return "3:" b_min "-" b_max " (04b)"
        if (s > b_max) return "4:" (b_max + 1) "+ (no step 4 range)"
        return "5:" s " (between 04a and 04b)"
    }
'

size_table() {  # size_table <tsv with size in column 1>
    awk -F '\t' \
        -v a_min="$step04a_cluster_min" -v a_max="$step04a_cluster_max" \
        -v b_min="$step04b_cluster_min" -v b_max="$step04b_cluster_max" \
        "$bin_of"'
        { n[bin($1)]++; if ($1 > largest) largest = $1 }
        END {
            for (k in n) print k "\t" n[k]
            print "9:largest cluster\t" largest + 0
        }
    ' "$1" | LC_ALL=C sort | while IFS=$'\t' read -r label n; do
        row "${label#*:}" "$n"
    done
}

printf '\n  cluster sizes:\n'
size_table "$work/clusters.tsv"

# ---------------------------------------------------------------------------
# step 03: alignment, screening, passed
# ---------------------------------------------------------------------------

passed_list="$step_03_dir/RNALalifold_passedList.txt"

printf '\n== step 03: alignment and screening (RNALalifold at %s C) ==\n' "$fold_temperature"

# one find per file type: name and size
find "$step_03_dir" -mindepth 2 -maxdepth 2 -type f -name '*_aligned.aln' \
    -printf '%f\t%s\n' 2> /dev/null | sed 's/_aligned\.aln\t/\t/' > "$work/aln.tsv" || true
find "$step_03_dir" -mindepth 2 -maxdepth 2 -type f -name '*_RNALalifold.out' \
    -printf '%f\t%s\n' 2> /dev/null | sed 's/_RNALalifold\.out\t/\t/' > "$work/fold.tsv" || true

awk -F '\t' '
    FILENAME == ARGV[1] { if ($2 > 0) aln[$1] = 1; next }
    FILENAME == ARGV[2] { if ($2 > 0) fold[$1] = 1; else empty[$1] = 1; next }
    {
        c = $2
        if (!(c in aln))       print c "\tnot aligned (03a)"
        else if (c in empty)   print c "\tRNALalifold output empty (03b failed)"
        else if (!(c in fold)) print c "\tnot screened (03b)"
    }
' "$work/aln.tsv" "$work/fold.tsv" "$work/clusters.tsv" > "$work/unfinished.tsv"

aligned=$(awk -F '\t' '$2 > 0' "$work/aln.tsv" | wc -l)
screened=$(awk -F '\t' '$2 > 0' "$work/fold.tsv" | wc -l)
unfinished=$(wc -l < "$work/unfinished.tsv")

row "aligned (03a)" "$aligned" "of $clusters_total"
row "screened (03b)" "$screened" "of $clusters_total"
row "not finished" "$unfinished" "$( ((unfinished > 0)) && echo '<- rerun 03a/03b for these' )"

if ((unfinished > 0)); then
    awk -F '\t' '{ n[$2]++ } END { for (k in n) printf "    %-40s %6d\n", k, n[k] }' "$work/unfinished.tsv"
    printf '\n'
    head -n "$shown" "$work/unfinished.tsv" | sed 's/^/    /'
    ((unfinished > shown)) && printf '    ... %d more (--list shows all)\n' "$((unfinished - shown))"
fi

if [[ ! -s "$passed_list" ]]; then
    row "passed (03c)" "-" "<- no $passed_list yet: run 03c"
    exit 0
fi

LC_ALL=C sort -u "$passed_list" > "$work/passed.txt"
passed=$(wc -l < "$work/passed.txt")

awk -v p="$passed" -v s="$screened" 'BEGIN { if (s > 0) printf "%.1f\n", 100 * p / s; else print "-" }' \
    > "$work/rate"
row "passed (03c)" "$passed" "($(cat "$work/rate")% of screened)"

# a passed list older than the newest screen may be missing clusters
newest_fold=$(find "$step_03_dir" -mindepth 2 -maxdepth 2 -type f -name '*_RNALalifold.out' \
    -newer "$passed_list" -print -quit 2> /dev/null || true)
[[ -n "$newest_fold" ]] &&
    printf '  note: some RNALalifold output is newer than the passed list -- rerun 03c\n'

# ---------------------------------------------------------------------------
# step 04 workload
# ---------------------------------------------------------------------------

printf '\n== step 04 workload: passed clusters by size ==\n'

# passed clusters with their size, and whether step 4 already has a result
LC_ALL=C join -t $'\t' -1 1 -2 2 -o 2.1,0 "$work/passed.txt" "$work/clusters.tsv" > "$work/passed_sized.tsv"

find "$step_04_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result.stk' -size +0 \
    -printf '%h\n' 2> /dev/null | sed 's|.*/||' | LC_ALL=C sort -u > "$work/done04.txt" || true

size_table "$work/passed_sized.tsv"

printf '\n'
printf '  %-24s %10s %10s %10s %12s\n' range passed "done" "to do" submissions
for part in a b; do
    if [[ $part == a ]]; then lo=$step04a_cluster_min hi=$step04a_cluster_max
    else lo=$step04b_cluster_min hi=$step04b_cluster_max; fi

    read -r in_range done_here < <(
        awk -F '\t' -v lo="$lo" -v hi="$hi" '
            FILENAME == ARGV[1] { done[$1] = 1; next }
            $1 >= lo && $1 <= hi { n++; if ($2 in done) d++ }
            END { print n + 0, d + 0 }
        ' "$work/done04.txt" "$work/passed_sized.tsv"
    )
    todo=$((in_range - done_here))
    rounds=$(((todo + step04_max_submit - 1) / step04_max_submit))
    printf '  %-24s %10s %10s %10s %12s\n' "04$part ($lo-$hi members)" "$in_range" "$done_here" "$todo" "$rounds"
done

beyond=$(awk -F '\t' -v hi="$step04b_cluster_max" '$1 > hi' "$work/passed_sized.tsv" | wc -l)
if ((beyond > 0)); then
    printf '\n  %d passed cluster(s) have more than %d members: no step 4 range covers them.\n' \
        "$beyond" "$step04b_cluster_max"
    printf '  Raise step04b_cluster_max (DATA_DIR/info.sh) or add a range to fold them.\n'
fi

printf '\n  "submissions" = rounds of submit.sh 04a/04b at step04_max_submit=%s per round.\n' \
    "$step04_max_submit"
