#!/bin/bash

# Score one alignment with RNA-SCoRE: how many of its sequences can actually
# form the consensus structure, and how much confidence that gives the
# alignment (High, Mid or Low).
#
# Step 5 runs this on the trimmed LocARNA motif; step 6 runs it again on the
# homolog hits with its own thresholds, which is why everything comes in as an
# argument.
#
# Everything runs with the cluster directory as the working directory.
# RNA-SCoRE.pl takes the cluster name from the input file name with SUFFIX
# removed, and writes <cluster>_pair_matrix.txt and
# <cluster>_detected_hairpins.txt into the current directory under that bare
# name, so two clusters sharing a working directory would overwrite each
# other's side files.
#
# All four thresholds are passed on every call: RNA-SCoRE prints "threshold not
# provided, using default" on STDOUT when one is missing, and STDOUT is also
# where the rank line this worker keeps comes from.
#
# Perl randomises hash order per process, and RNA-SCoRE iterates a hash to pick
# which of a set of identical motif sequences to keep. PERL_HASH_SEED pins that,
# so the cleaned alignment is the same on a rerun; the rank is the same either
# way.
#
# Writes next to the alignment, where <base> is ALIGNMENT without .sto:
#   <base>_cleaned.sto     the sequences that passed (only when at least two did)
#   <base>_evaluated.tsv   per-sequence verdicts (only when the motif had >= 5 bp)
#   <base>_rnascore.err    STDERR, kept only when it is not empty
#   <base>_rank.tsv        the rank line without RNA-SCoRE's header, written
#                          last and only on success
#
# Arguments:
#   RNA-SCoRE.pl, cluster directory, alignment, suffix to strip (-e),
#   motif threshold (--mt), stem threshold (-t), GC threshold (--gc),
#   duplicates (-d), Perl hash seed

set -euo pipefail

if [[ $# -ne 9 ]]; then
    echo \
        "Usage: $0 SCORE_PL CLUSTER_DIR ALIGNMENT SUFFIX MT BP GC DUPL HASH_SEED" \
        >&2
    exit 1
fi

score_pl=$1
cluster_dir=$2
alignment=$3
suffix=$4
motif_threshold=$5
bp_threshold=$6
gc_threshold=$7
duplicates=$8
hash_seed=$9

cd "$cluster_dir"

[[ -s "$alignment" ]] || {
    echo "Alignment not found: $cluster_dir/$alignment" >&2
    exit 1
}

base=${alignment%.sto}
cleaned="${base}_cleaned.sto"
evaluated="${base}_evaluated.tsv"
errors="${base}_rnascore.err"
rank="${base}_rank.tsv"

rm -f "$cleaned" "$evaluated" "$errors" "$rank" "$rank.tmp" "$rank.line"

if ! PERL_HASH_SEED="$hash_seed" \
    perl "$score_pl" \
        -e "$suffix" \
        -d "$duplicates" \
        --mt "$motif_threshold" \
        -t "$bp_threshold" \
        --gc "$gc_threshold" \
        "$alignment" \
        "$cleaned" \
        "$evaluated" \
        > "$rank.tmp" \
        2> "$errors"
then
    echo "RNA-SCoRE failed for $cluster_dir/$alignment, see $errors" >&2
    rm -f "$rank.tmp"
    exit 1
fi

# STDOUT is a column header followed by the one line this alignment is about.
# Keeping only that line lets the step script concatenate every cluster under a
# single header.
grep -v '^clusterFile' "$rank.tmp" | grep -v '^[[:space:]]*$' > "$rank.line" || true
rm -f "$rank.tmp"

# A run that exits 0 without a rank line has not evaluated anything, which is a
# failure however it managed to happen.
if [[ ! -s "$rank.line" ]]; then
    echo "RNA-SCoRE produced no rank line for $cluster_dir/$alignment" >&2
    rm -f "$rank.line"
    exit 1
fi

mv "$rank.line" "$rank"

[[ -s "$errors" ]] || rm -f "$errors"
