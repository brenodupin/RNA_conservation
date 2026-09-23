#!/bin/bash

# Turn one cmsearch result into an alignment of its hits: extract the hit
# sequences with cmsearch_reformatv1_1.pl, then align them to the covariance
# model with cmalign. Both run in the rnatools container.
#
# cmalign writes Pfam-style Stockholm (one line per sequence) rather than the
# wrapped default: RNA-SCoRE reassembles wrapped alignments by matching names
# as regular expressions, and hit names like NC_000932.1/5906-6060_-1 are full
# of characters that mean something there. R-scape reads either form.
#
# cmalign is kept to one thread, because the step script already runs one of
# these per CPU.
#
# Writes in the cluster directory, <prefix> being e.g. <cluster>_hits:
#   <prefix>.fna          the hit sequences, from the reformat script
#   <prefix>_scores.tab   the hit scores, from the reformat script
#   <prefix>.sto          the hits aligned to the model -- written last
#   <prefix>_cmalign.log  what cmalign reported
#
# Arguments:
#   container, cluster directory, model (.cm), cmsearch output, output prefix,
#   directory of the riboswitch helper scripts inside the container

set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 CONTAINER CLUSTER_DIR MODEL HITS_TXT PREFIX SCRIPTS_DIR" >&2
    exit 1
fi

rnatools=$1
cluster_dir=$2
model=$3
hits_text=$4
prefix=$5
scripts_dir=$6

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

for file in "$model" "$hits_text"; do
    [[ -s "$cluster_dir/$file" ]] || {
        echo "Not found: $cluster_dir/$file" >&2
        exit 1
    }
done

fasta="${prefix}.fna"
scores="${prefix}_scores.tab"
alignment="${prefix}.sto"
log="$cluster_dir/${prefix}_cmalign.log"

rm -f \
    "$cluster_dir/$fasta" "$cluster_dir/$scores" \
    "$cluster_dir/$alignment" "$cluster_dir/$alignment.tmp" "$log"

in_container() {
    apptainer exec \
        --bind "$cluster_dir:/work" \
        --pwd /work \
        "$rnatools" \
        "$@"
}

if ! in_container perl "$scripts_dir/cmsearch_reformatv1_1.pl" \
    -s "$scores" "$hits_text" "$fasta" > "$log" 2>&1
then
    echo "cmsearch_reformatv1_1.pl failed for $cluster_dir/$hits_text, see $log" >&2
    exit 1
fi

if ! grep -q '^>' "$cluster_dir/$fasta" 2> /dev/null; then
    echo "No hit sequences extracted from $cluster_dir/$hits_text" >&2
    exit 1
fi

if ! in_container cmalign \
    --cpu 1 \
    --outformat pfam \
    -o "$alignment.tmp" \
    "$model" \
    "$fasta" \
    >> "$log" 2>&1
then
    echo "cmalign failed for $cluster_dir/$fasta, see $log" >&2
    rm -f "$cluster_dir/$alignment.tmp"
    exit 1
fi

mv "$cluster_dir/$alignment.tmp" "$cluster_dir/$alignment"
