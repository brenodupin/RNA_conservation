#!/bin/bash
#SBATCH --job-name=rnac_06a_launch
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bdupin@uwo.ca

# Pick the seed alignments for this round of the homolog search, make sure the
# search database is in place, and submit worker_06a.sh as a Slurm array over
# them -- one task per cluster, each building, calibrating and searching one
# covariance model. Like step_04a.sh this job is only the launcher: it writes
# the manifest, submits the array, prints where to find it, and exits.
#
# Outputs go to DATA_DIR/06_homolog/round_<step06_round>/<cluster>.
#
# Seeds. In round 1 they come from step 5: every cluster whose R-scape result
# has at least covariation_min_percent of its base pairs covarying (the list
# step_05b.sh writes as passed_covariation.txt), seeded with its
# <cluster>_motif_cleaned.sto. --list FILE narrows or overrides that, one
# cluster per line (a listed cluster below the cutoff is still searched, as
# long as it has a cleaned alignment). Later rounds seed from step 7 and need
# --seeds FILE, one "<cluster>\t<alignment path>" per line.
#
# Database. step06_search_db in info.sh, or, when that is empty, every genome
# in 00_oneline concatenated into 06_homolog/search_db.fa, rebuilt only when
# its content would change. Do not let it change while an array from an
# earlier run is still searching it.
#
# A cluster is left out when its <cluster>_result_06a.txt marker already
# records the same seed md5, database md5 and search settings, so a rerun
# submits only what is new or changed. --force submits everything, and also
# rebuilds and recalibrates models that could have been reused. More clusters
# than step06_max_submit are submitted in batches: rerun once one drains.
#
# Usage:
#   sharcnet/submit.sh 06a DATA_DIR [--force] [--list FILE] [--seeds FILE]

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo \
        "Usage: sbatch $0 DATA_DIR SHARCNET_DIR [--force] [--list FILE] [--seeds FILE]" \
        >&2
    exit 1
fi

data_dir=$1
script_dir=$2
shift 2

force=no
list_file=
seeds_file=

while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            force=yes
            ;;
        --list | --seeds)
            [[ $# -ge 2 ]] || {
                echo "$1 needs a file" >&2
                exit 1
            }
            if [[ $1 == --list ]]; then list_file=$2; else seeds_file=$2; fi
            shift
            ;;
        *)
            echo "Unknown argument: $1 (expected --force, --list FILE or --seeds FILE)" >&2
            exit 1
            ;;
    esac
    shift
done

source "$script_dir/info.sh"

worker="$script_dir/worker_06a.sh"
round_dir="$step_06_dir/round_$step06_round"

[[ -x "$worker" ]] || {
    echo "Worker is not executable: $worker" >&2
    exit 1
}

[[ "$step06_round" =~ ^[1-9][0-9]*$ ]] || {
    echo "step06_round must be a positive whole number, got '$step06_round'" >&2
    exit 1
}

# ---------------------------------------------------------------------------
# seeds
# ---------------------------------------------------------------------------

declare -A seed=()

if [[ -n "$seeds_file" ]]; then
    [[ -s "$seeds_file" ]] || {
        echo "Seed list not found or empty: $seeds_file" >&2
        exit 1
    }

    while IFS=$'\t' read -r cluster path _; do
        [[ -z "$cluster" || "$cluster" == \#* ]] && continue
        [[ -n "${path:-}" ]] || {
            echo "No seed path for $cluster in $seeds_file" >&2
            exit 1
        }
        seed[$cluster]=$path
    done < "$seeds_file"
elif ((step06_round == 1)); then
    [[ -d "$step_05_dir" ]] || {
        echo "Missing step 5 output: $step_05_dir -- run step 05 first." >&2
        exit 1
    }

    if [[ -n "$list_file" ]]; then
        # A hand-picked list overrides the covariation cutoff.
        while IFS= read -r cluster; do
            seed[$cluster]="$step_05_dir/$cluster/${cluster}_motif_cleaned.sto"
        done < <(
            find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f \
                -name '*_motif_cleaned.sto' -printf '%P\n' |
                sed 's|/.*||'
        )
    else
        while IFS=$'\t' read -r cluster percent; do
            awk -v p="$percent" -v c="$covariation_min_percent" \
                'BEGIN { exit !(p != "" && p + 0 >= c + 0) }' || continue
            seed[$cluster]="$step_05_dir/$cluster/${cluster}_motif_cleaned.sto"
        done < <(
            find "$step_05_dir" -mindepth 2 -maxdepth 2 -type f \
                -name '*_result_05b.txt' -print0 |
                xargs -0 -r awk -F '\t' '
                    FNR == 1 {
                        split(FILENAME, path, "/")
                        cluster = path[length(path) - 1]
                    }

                    $1 == "percent_covarying" {
                        print cluster "\t" $2
                    }
                '
        )
    fi
else
    echo \
        "Round $step06_round seeds from step 7: pass them with --seeds FILE" \
        "(one '<cluster><TAB><alignment>' per line)." \
        >&2
    exit 1
fi

# --list keeps only the listed clusters.
if [[ -n "$list_file" ]]; then
    [[ -s "$list_file" ]] || {
        echo "Cluster list not found or empty: $list_file" >&2
        exit 1
    }

    declare -A listed=()

    while IFS= read -r cluster; do
        listed[$cluster]=yes
    done < <(grep -vE '^[[:space:]]*(#.*)?$' "$list_file" | awk '{ print $1 }')

    for cluster in "${!seed[@]}"; do
        [[ -n "${listed[$cluster]:-}" ]] || unset "seed[$cluster]"
    done

    for cluster in "${!listed[@]}"; do
        [[ -n "${seed[$cluster]:-}" ]] ||
            echo "Listed but no seed alignment, skipped: $cluster" >&2
    done
fi

mapfile -t clusters < <(printf '%s\n' "${!seed[@]}" | grep -v '^$' | LC_ALL=C sort)

cluster_total=${#clusters[@]}

if ((cluster_total == 0)); then
    echo \
        "No seed alignments for round $step06_round" \
        "(covariation cutoff $covariation_min_percent%) -- nothing to submit." \
        >&2
    exit 0
fi

missing=0

for cluster in "${clusters[@]}"; do
    if [[ ! -s "${seed[$cluster]}" ]]; then
        echo "Seed alignment not found, skipped: ${seed[$cluster]}" >&2
        unset "seed[$cluster]"
        missing=$((missing + 1))
    fi
done

if ((missing > 0)); then
    mapfile -t clusters < <(printf '%s\n' "${!seed[@]}" | grep -v '^$' | LC_ALL=C sort)
    cluster_total=${#clusters[@]}
    ((cluster_total > 0)) || exit 1
fi

# ---------------------------------------------------------------------------
# search database
# ---------------------------------------------------------------------------

mkdir -p "$round_dir"

if [[ -n "$step06_search_db" ]]; then
    search_db=$(realpath "$step06_search_db")

    [[ -s "$search_db" ]] || {
        echo "Search database not found: $step06_search_db" >&2
        exit 1
    }

    db_md5=$(md5sum < "$search_db" | cut -d ' ' -f 1)
else
    search_db="$step_06_dir/search_db.fa"
    db_note="$step_06_dir/search_db.md5"

    shopt -s nullglob
    genomes=("$step_00_dir"/*_oneLine.fasta)
    shopt -u nullglob

    ((${#genomes[@]} > 0)) || {
        echo "No genomes in $step_00_dir to build the search database from." >&2
        exit 1
    }

    # Built only when its content would change, so a rerun keeps the file
    # (and its md5, and therefore every search already done) as it is.
    db_md5=$(cat "${genomes[@]}" | md5sum | cut -d ' ' -f 1)

    if [[ ! -s "$search_db" || ! -s "$db_note" || "$(cat "$db_note")" != "$db_md5" ]]; then
        echo "Building the search database from ${#genomes[@]} genomes in $step_00_dir"
        cat "${genomes[@]}" > "$search_db.tmp"
        mv "$search_db.tmp" "$search_db"
        printf '%s\n' "$db_md5" > "$db_note"
    fi
fi

settings=$(
    printf 'evalue=%s cmsearch=%s cmcalibrate=%s' \
        "$step06_cmsearch_evalue" \
        "$step06_cmsearch_options" \
        "$step06_cmcalibrate_options"
)

# ---------------------------------------------------------------------------
# what still needs doing
# ---------------------------------------------------------------------------

seeds=()

for cluster in "${clusters[@]}"; do
    seeds+=("${seed[$cluster]}")
done

declare -A seed_md5=()

while read -r sum path; do
    seed_md5[$path]=$sum
done < <(printf '%s\0' "${seeds[@]}" | xargs -0 -r md5sum)

declare -A marker_seed=()
declare -A marker_db=()
declare -A marker_settings=()

while IFS=$'\t' read -r cluster field value; do
    case $field in
        seed_md5) marker_seed[$cluster]=$value ;;
        db_md5)   marker_db[$cluster]=$value ;;
        settings) marker_settings[$cluster]=$value ;;
    esac
done < <(
    find "$round_dir" -mindepth 2 -maxdepth 2 -type f -name '*_result_06a.txt' -print0 |
        xargs -0 -r awk -F '\t' '
            FNR == 1 {
                split(FILENAME, path, "/")
                cluster = path[length(path) - 1]
            }

            $1 == "seed_md5" || $1 == "db_md5" || $1 == "settings" {
                print cluster "\t" $1 "\t" $2
            }
        '
)

# One manifest per run, named down to the second: the retry command printed
# below points at it, so a later run the same day must not overwrite it.
run_tag="step_06a_r${step06_round}_$(date +%b_%d_%H%M%S)"
[[ "$force" == yes ]] && run_tag="${run_tag}_force"

manifest="$round_dir/${run_tag}.tsv"
unchanged=0

: > "$manifest.tmp"

for cluster in "${clusters[@]}"; do
    path=${seed[$cluster]}
    md5=${seed_md5[$path]:-}

    if [[ "$force" != yes &&
        "${marker_seed[$cluster]:-}" == "$md5" &&
        "${marker_db[$cluster]:-}" == "$db_md5" &&
        "${marker_settings[$cluster]:-}" == "$settings" ]]
    then
        unchanged=$((unchanged + 1))
        continue
    fi

    # The directory makes a queued cluster show up as pending in
    # step_06_status.tsv before its task has started.
    mkdir -p "$round_dir/$cluster"
    printf '%s\t%s\t%s\n' "$cluster" "$path" "$md5" >> "$manifest.tmp"
done

todo_total=$(wc -l < "$manifest.tmp")

printf 'Step 06a, round %s: %s seed(s), %s already searched, %s to submit.\n' \
    "$step06_round" "$cluster_total" "$unchanged" "$todo_total"
printf '  Database:  %s (md5 %s)\n' "$search_db" "$db_md5"
printf '  Settings:  %s\n' "$settings"

if ((todo_total == 0)); then
    rm -f "$manifest.tmp"
    printf '\nNothing to submit.\n'
    exit 0
fi

mv "$manifest.tmp" "$manifest"

submit_count=$todo_total
remaining_after=0

if ((todo_total > step06_max_submit)); then
    submit_count=$step06_max_submit
    remaining_after=$((todo_total - step06_max_submit))
fi

log_dir="$logs_dir/$run_tag"
mkdir -p "$log_dir"

worker_args=("$manifest" "$data_dir" "$script_dir" "$search_db" "$db_md5" "$settings")
[[ "$force" == yes ]] && worker_args+=(--force)

submission=$(
    sbatch \
        --parsable \
        --job-name=rnac_06a \
        --array="1-${submit_count}%${step06_jobs}" \
        --cpus-per-task="$step06a_cpus" \
        --mem="$step06a_mem" \
        --time="$step06a_time" \
        --mail-type=END,FAIL \
        --mail-user=bdupin@uwo.ca \
        --output="$log_dir/%A_%a.out" \
        "$worker" \
        "${worker_args[@]}"
)

job_id=${submission%%;*}

printf '  Manifest:  %s\n' "$manifest"
printf '  Array job: %s\n' "$job_id"
printf '  Array:     1-%s (%%%s concurrent)\n' "$submit_count" "$step06_jobs"
printf '  Resources: %s cpus, %s mem, %s time per task\n' \
    "$step06a_cpus" "$step06a_mem" "$step06a_time"
printf '  Logs:      %s/%s_*.out\n' "$log_dir" "$job_id"
printf '\n'
printf '  To retry failed tasks of this run (failed indices: sacct -j %s):\n' "$job_id"
printf '    sbatch --array=<indices> --cpus-per-task=%s --mem=%s --time=%s \\\n' \
    "$step06a_cpus" "$step06a_mem" "$step06a_time"
printf '        --output=%q %q' "$log_dir/%A_%a.out" "$worker"
printf ' %q' "${worker_args[@]}"
printf '\n'

if ((remaining_after > 0)); then
    # Without --force even if this run had it: the batch submitted now leaves
    # fresh markers, so the rerun picks up exactly the clusters left over.
    rerun_cmd="sharcnet/submit.sh 06a $data_dir"
    [[ -n "$list_file" ]] && rerun_cmd="$rerun_cmd --list $list_file"
    [[ -n "$seeds_file" ]] && rerun_cmd="$rerun_cmd --seeds $seeds_file"

    printf '\n'
    printf '  NOTE: %s clusters exceed the %s-submission cap (AssocMaxSubmitJobLimit=1000).\n' \
        "$todo_total" "$step06_max_submit"
    printf '        Submitted %s now; %s left for next time.\n' \
        "$submit_count" "$remaining_after"
    printf '        Once this batch drains, rerun: %s\n' "$rerun_cmd"
fi
