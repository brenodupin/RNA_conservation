#!/bin/bash

# Check a finished step 04 array (step_04a.sh or step_04b.sh): what Slurm says
# happened to every task, whether each cluster really has a usable
# <cluster>_result.stk, how long the clusters took and how much memory they
# used against what was requested, and what is left to submit. Read-only.
#
# Usage:
#   sharcnet/check_04.sh DATA_DIR JOB_ID [a|b]
#
# e.g. sharcnet/check_04.sh data_protists_plt_July_31 22465833 a

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 DATA_DIR JOB_ID [a|b]" >&2
    exit 1
fi

data_dir=$1
job=$2
part=${3:-a}

sharcnet_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$sharcnet_dir/info.sh"

case $part in
    a) cluster_min=$step04a_cluster_min; cluster_max=$step04a_cluster_max
       req_mem=$step04a_mem; req_time=$step04a_time; req_cpus=$step04a_cpus ;;
    b) cluster_min=$step04b_cluster_min; cluster_max=$step04b_cluster_max
       req_mem=$step04b_mem; req_time=$step04b_time; req_cpus=$step04b_cpus ;;
    *) echo "Part must be a or b, got '$part'" >&2; exit 1 ;;
esac

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# ---------------------------------------------------------------------------
# which clusters this job ran: the task logs say
# ---------------------------------------------------------------------------

shopt -s nullglob
logs=("$logs_dir"/*/"${job}"_*.out)
shopt -u nullglob

if ((${#logs[@]} == 0)); then
    echo "No task logs for job $job under $logs_dir/*/" >&2
    exit 1
fi

for log in "${logs[@]}"; do
    task=${log##*_}
    task=${task%.out}
    awk -v task="$task" '
        /^Cluster: / { cluster = $2; members = $3; gsub(/[()]/, "", members) }
        /Already completed, skipping/ { skipped = 1 }
        END { printf "%s\t%s\t%s\t%s\n", task, cluster, members, skipped ? "skipped" : "ran" }
    ' "$log"
done | sort -n > "$work/tasks.tsv"

printf 'Step 04%s, job %s: %s task logs in %s\n\n' \
    "$part" "$job" "${#logs[@]}" "$(dirname "${logs[0]}")"

# ---------------------------------------------------------------------------
# 1. Slurm's view
# ---------------------------------------------------------------------------

printf '== 1. Slurm (sacct) ==\n'

if command -v sacct > /dev/null; then
    sacct -j "$job" -n -P \
        --format=JobID,State,ExitCode,ElapsedRaw,MaxRSS,TotalCPU,AllocCPUS \
        > "$work/sacct.txt" || true

    awk -F '|' -v req_mem="$req_mem" -v req_cpus="$req_cpus" '
        function seconds(t,    d, parts, n, s) {
            d = 0
            if (t ~ /-/) { split(t, parts, "-"); d = parts[1]; t = parts[2] }
            n = split(t, parts, ":")
            s = 0
            for (i = 1; i <= n; i++) s = s * 60 + parts[i]
            return d * 86400 + s
        }

        function megabytes(m,    v) {
            v = m + 0
            if (m ~ /K$/) return v / 1024
            if (m ~ /G$/) return v * 1024
            if (m ~ /T$/) return v * 1024 * 1024
            if (m ~ /M$/) return v
            return v / 1048576
        }

        $1 !~ /\./ {
            task = $1; sub(/.*_/, "", task)
            state[task] = $2
            split($2, first, " ")
            count[first[1]]++
            elapsed[task] = $4
        }

        $1 ~ /\.batch$/ {
            task = $1; sub(/\.batch$/, "", task); sub(/.*_/, "", task)
            rss[task] = megabytes($5)
            cpu[task] = seconds($6)
            cpus[task] = $7
        }

        END {
            for (s in count) printf "  %-14s %6d\n", s, count[s]

            n = 0
            for (t in rss) {
                n++
                if (rss[t] > max_rss) { max_rss = rss[t]; max_rss_task = t }
                if (elapsed[t] > 0 && cpus[t] > 0) {
                    eff_sum += cpu[t] / (elapsed[t] * cpus[t])
                    eff_n++
                }
            }

            req = megabytes(req_mem)
            if (n > 0) {
                printf "  peak memory:   %.0f MB (task %s) of %.0f MB requested\n", \
                    max_rss, max_rss_task, req
                over = 0
                for (t in rss) if (rss[t] > 0.9 * req) over++
                printf "  tasks above 90%% of the memory request: %d\n", over
            }
            if (eff_n > 0)
                printf "  mean CPU efficiency: %.0f%% of %s cpus\n", 100 * eff_sum / eff_n, req_cpus

            bad = ""
            for (t in state) if (state[t] !~ /^COMPLETED/) bad = bad " " t ":" state[t]
            if (bad != "") printf "  not COMPLETED:%s\n", bad
        }
    ' "$work/sacct.txt"
else
    printf '  sacct not available here\n'
fi

# ---------------------------------------------------------------------------
# 2. the results themselves
# ---------------------------------------------------------------------------
#
# A usable result.stk starts with "# STOCKHOLM", has a consensus structure,
# ends with "//", holds one sequence per cluster member, and has no LocARNA
# stand-in name (S000001...) left, i.e. worker_04.sh restored every identifier.

printf '\n== 2. Results ==\n'

: > "$work/problems.tsv"
ok=0

while IFS=$'\t' read -r task cluster members how; do
    result="$step_04_dir/$cluster/${cluster}_result.stk"

    if [[ -z "$cluster" ]]; then
        printf '%s\t-\tno cluster name in the task log\n' "$task" >> "$work/problems.tsv"
        continue
    fi

    if [[ ! -s "$result" ]]; then
        reason="no result.stk"
        timing="$step_04_dir/$cluster/${cluster}_timing.tsv"
        if [[ -s "$timing" ]]; then
            code=$(cut -f6 "$timing")
            [[ "$code" == 124 ]] && reason="no result.stk: timed out after $locarna_timeout" ||
                reason="no result.stk: mlocarna exit $code"
        fi
        printf '%s\t%s\t%s\n' "$task" "$cluster" "$reason" >> "$work/problems.tsv"
        continue
    fi

    verdict=$(
        awk -v members="$members" '
            NR == 1 && $0 !~ /^# STOCKHOLM/ { bad = bad "no STOCKHOLM header; " }
            /^#=GC SS_cons/ { ss = 1 }
            NF { last = $0 }
            $0 !~ /^#/ && $0 !~ /^\/\// && NF >= 2 {
                names[$1] = 1
                if ($1 ~ /^S[0-9]{6}$/) stand_in = 1
            }
            END {
                n = 0
                for (k in names) n++
                if (!ss) bad = bad "no SS_cons; "
                if (last != "//") bad = bad "no closing //; "
                if (stand_in) bad = bad "LocARNA stand-in names left; "
                if (members != "" && n != members) bad = bad n " sequences for " members " members; "
                print bad
            }
        ' "$result"
    )

    if [[ -n "$verdict" ]]; then
        printf '%s\t%s\t%s\n' "$task" "$cluster" "${verdict%; }" >> "$work/problems.tsv"
    else
        ok=$((ok + 1))
    fi
done < "$work/tasks.tsv"

tasks=$(wc -l < "$work/tasks.tsv")
skipped=$(awk -F '\t' '$4 == "skipped"' "$work/tasks.tsv" | wc -l)
problems=$(wc -l < "$work/problems.tsv")

printf '  tasks:                    %6d (%d skipped as already done)\n' "$tasks" "$skipped"
printf '  usable result.stk:        %6d\n' "$ok"
printf '  problems:                 %6d\n' "$problems"

if ((problems > 0)); then
    printf '\n  task\tcluster\tproblem\n'
    sort -n "$work/problems.tsv" | head -n 40 | sed 's/^/  /'
    ((problems > 40)) && printf '  ... %d more\n' "$((problems - 40))"
fi

# ---------------------------------------------------------------------------
# 3. how long clusters took, by size
# ---------------------------------------------------------------------------

printf '\n== 3. mlocarna wall time by cluster size (from *_timing.tsv) ==\n'

while IFS=$'\t' read -r _ cluster _ how; do
    [[ "$how" == ran ]] || continue
    timing="$step_04_dir/$cluster/${cluster}_timing.tsv"
    [[ -s "$timing" ]] && cat "$timing"
done < "$work/tasks.tsv" > "$work/timings.tsv"

if [[ -s "$work/timings.tsv" ]]; then
    sort -t $'\t' -k2,2n -k5,5n "$work/timings.tsv" |
        awk -F '\t' '
            function hms(s) { return sprintf("%d:%02d:%02d", s / 3600, (s % 3600) / 60, s % 60) }

            {
                size = $2
                t[size, ++n[size]] = $5
                total += $5
                if ($5 > longest) { longest = $5; longest_cluster = $1 " (" $2 " seqs)" }
                if ($6 != 0) failed++
            }

            END {
                printf "  %-8s %6s %10s %10s %10s\n", "members", "done", "median", "p90", "max"
                for (s = 1; s <= 200; s++) {
                    if (!(s in n)) continue
                    m = n[s]
                    printf "  %-8s %6d %10s %10s %10s\n", s, m, \
                        hms(t[s, int((m + 1) / 2)]), hms(t[s, int(0.9 * m + 0.999)]), hms(t[s, m])
                }
                printf "  longest: %s, %s\n", longest_cluster, hms(longest)
                printf "  total mlocarna time: %.1f core-hours\n", total / 3600
                if (failed) printf "  non-zero exit codes: %d\n", failed
            }
        '
    printf '  (requested per task: %s cpus, %s, %s; mlocarna timeout %s)\n' \
        "$req_cpus" "$req_mem" "$req_time" "$locarna_timeout"
else
    printf '  no timing files found\n'
fi

# ---------------------------------------------------------------------------
# 4. leftovers from tasks that did not finish
# ---------------------------------------------------------------------------

printf '\n== 4. Leftovers in %s ==\n' "$step_04_dir"

leftover_dirs=$(find "$step_04_dir" -mindepth 2 -maxdepth 2 -type d -name '*_locarnap' | wc -l)
leftover_tmp=$(find "$step_04_dir" -mindepth 2 -maxdepth 2 -type f -name '*.tmp' | wc -l)
printf '  *_locarnap working directories: %d (a finished cluster removes its own)\n' "$leftover_dirs"
printf '  *.tmp files:                    %d\n' "$leftover_tmp"

# ---------------------------------------------------------------------------
# 5. what is left in this size range
# ---------------------------------------------------------------------------

printf '\n== 5. Clusters with %s-%s members that passed step 3 ==\n' "$cluster_min" "$cluster_max"

passed_list="$step_03_dir/RNALalifold_passedList.txt"
cluster_count="$step_02_dir/cluster_count.tsv"

if [[ -s "$passed_list" && -s "$cluster_count" ]]; then
    awk -F '\t' -v min="$cluster_min" -v max="$cluster_max" '
        NR == FNR { passed[$1] = 1; next }
        $2 in passed && $1 >= min && $1 <= max { print $2 }
    ' "$passed_list" "$cluster_count" > "$work/in_range.txt"

    in_range=$(wc -l < "$work/in_range.txt")
    have=0

    while IFS= read -r cluster; do
        [[ -s "$step_04_dir/$cluster/${cluster}_result.stk" ]] && have=$((have + 1))
    done < "$work/in_range.txt"

    # A result.stk that exists but failed the checks above still counts as
    # done for step_04a.sh/04b.sh, which skip any non-empty result.stk.
    broken=$(awk -F '\t' '$3 !~ /^no / { n++ } END { print n + 0 }' "$work/problems.tsv")

    printf '  in range:        %6d\n' "$in_range"
    printf '  with a result:   %6d\n' "$have"
    printf '  still to do:     %6d\n' "$((in_range - have))"

    if ((in_range > have)); then
        printf '  -> sharcnet/submit.sh 04%s %s  submits exactly those\n' "$part" "$data_dir"
    fi

    if ((broken > 0)); then
        printf '  %d of the results above are broken (section 2) and will not be\n' "$broken"
        printf '  resubmitted as they are: delete their <cluster>_result.stk first.\n'
    fi
else
    printf '  cluster_count.tsv or RNALalifold_passedList.txt missing\n'
fi
