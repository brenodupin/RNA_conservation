#!/bin/bash

# Test one alignment for covariation with R-scape, in the rnatools container,
# and summarise what it found.
#
# Steps 6 and 7 run R-scape on their own alignments (the homolog hits, and the
# HMMER alignment with and without the LocARNA structure), so the alignment,
# the output names and the R-scape options are all arguments. OPTIONS is a
# single string of R-scape flags, e.g. "-s --cacofold".
#
# R-scape refuses to start when --outdir does not exist and writes into it
# without clearing it first, so the output directory is recreated here and a
# rerun never mixes old and new results.
#
# The summary is read from R-scape's power file, which reports
#   # BPAIRS 11
#   # BPAIRS expected to covary 1.9 +/- 1.2
#   # BPAIRS observed to covary 2
# Those are the numbers step 6 keeps motifs on (at least 5% of the base pairs
# covarying). It is written last, so its presence is what tells the step script
# that this alignment is done.
#
# Writes in the cluster directory:
#   <out name>.log                 R-scape's own output
#   <out dir>/                     everything R-scape produced
#   <out name>_covariation.tsv     alignment, base pairs, expected covarying,
#                                  its standard deviation, observed covarying,
#                                  percent covarying -- written last
#
# Arguments:
#   container, cluster directory, alignment, output directory, output name,
#   seed, timeout, R-scape options

set -euo pipefail

if [[ $# -ne 8 ]]; then
    echo \
        "Usage: $0 CONTAINER CLUSTER_DIR ALIGNMENT OUT_DIR OUT_NAME SEED TIMEOUT OPTIONS" \
        >&2
    exit 1
fi

rnatools=$1
cluster_dir=$2
alignment=$3
out_dir=$4
out_name=$5
seed=$6
run_timeout=$7
options=$8

[[ -n "$cluster_dir" && -n "$out_dir" && -n "$out_name" ]] || {
    echo "Cluster directory, output directory and output name cannot be empty" >&2
    exit 1
}

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

[[ -s "$cluster_dir/$alignment" ]] || {
    echo "Alignment not found: $cluster_dir/$alignment" >&2
    exit 1
}

log="$cluster_dir/${out_name}.log"
covariation="$cluster_dir/${out_name}_covariation.tsv"
power="$cluster_dir/$out_dir/${out_name}.power"

rm -rf "${cluster_dir:?}/${out_dir:?}"
rm -f "$log" "$covariation" "$covariation.tmp"
mkdir -p "$cluster_dir/$out_dir"

exit_code=0

# $options is deliberately unquoted: it carries the R-scape flags as a list.
# shellcheck disable=SC2086
timeout \
    --signal=TERM \
    --kill-after=1m \
    "$run_timeout" \
    apptainer exec \
        --bind "$cluster_dir:/work" \
        --pwd /work \
        "$rnatools" \
        R-scape \
        $options \
        --outdir "$out_dir" \
        --outname "$out_name" \
        --voutput \
        --seed "$seed" \
        "$alignment" \
        > "$log" 2>&1 || exit_code=$?

if ((exit_code == 124 || exit_code == 137)); then
    echo "R-scape timed out after $run_timeout for $cluster_dir/$alignment" >&2
    exit 1
fi

if ((exit_code != 0)); then
    echo "R-scape failed for $cluster_dir/$alignment: exit=$exit_code, see $log" >&2
    exit 1
fi

# Only the two-set test (-s) produces a power file, so an alignment folded
# without it is recorded with empty numbers rather than treated as a failure.
if [[ -s "$power" ]]; then
    awk -v alignment="$alignment" '
        /^# BPAIRS [0-9]+$/ {
            bpairs = $NF
            next
        }

        /^# BPAIRS expected to covary / {
            expected = $(NF - 2)
            deviation = $NF
            next
        }

        /^# BPAIRS observed to covary [0-9]+$/ {
            covarying = $NF
            next
        }

        END {
            percent = bpairs > 0 ? 100 * covarying / bpairs : 0

            printf "%s\t%d\t%s\t%s\t%d\t%.2f\n", \
                alignment, bpairs, expected, deviation, covarying, percent
        }
    ' "$power" > "$covariation.tmp"
else
    echo "No power file for $cluster_dir/$alignment: $power" >&2
    printf '%s\t\t\t\t\t\n' "$alignment" > "$covariation.tmp"
fi

mv "$covariation.tmp" "$covariation"
