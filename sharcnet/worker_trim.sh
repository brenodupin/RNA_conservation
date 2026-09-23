#!/bin/bash

# Trim one LocARNA alignment down to its motif region with trimAlignment.pl.
#
# The Stockholm file is written through a temporary name and only moved into
# place once it holds a consensus structure and at least one sequence.
# trimAlignment.pl cannot report that it produced nothing: when the LocARNA
# alignment fits in a single 120-column block it writes a header, a consensus
# line and no sequences at all (@aln_entires in its else branch), which
# RNA-SCoRE would then rank "Low" with zero sequences as if that were a real
# result. Such an alignment is kept as <cluster>_motif.sto.invalid for
# inspection and this worker exits non-zero instead.
#
# The FASTA (RNAz) and aligned-FASTA (SQUARNA) outputs are only kept with
# optional outputs on. trimAlignment.pl always wants three output names, so the
# two go to /dev/null otherwise. The clustal view (.aln) comes from
# esl-reformat, which lives in the container rather than on the host.
#
# Arguments:
#   trimAlignment.pl, LocARNA result.stk, cluster directory, cluster name,
#   optional outputs (yes|no), container

set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo \
        "Usage: $0 TRIM_PL RESULT_STK CLUSTER_DIR CLUSTER OPTIONAL_OUTPUTS CONTAINER" \
        >&2
    exit 1
fi

trim_pl=$1
result_stk=$2
cluster_dir=$3
cluster=$4
optional_outputs=$5
rnatools=$6

motif="$cluster_dir/${cluster}_motif.sto"
fasta="$cluster_dir/${cluster}_motif.fasta"
afa="$cluster_dir/${cluster}_motif.afa"
alignment="$cluster_dir/${cluster}_motif.aln"

[[ -s "$result_stk" ]] || {
    echo "LocARNA result not found: $result_stk" >&2
    exit 1
}

mkdir -p "$cluster_dir"

# Clear anything an interrupted run left behind, so a partial file can never be
# mistaken for finished work.
rm -f \
    "$motif" "$motif.tmp" "$motif.invalid" \
    "$fasta" "$fasta.tmp" \
    "$afa" "$afa.tmp" \
    "$alignment"

if [[ "$optional_outputs" == yes ]]; then
    fasta_out="$fasta.tmp"
    afa_out="$afa.tmp"
else
    fasta_out=/dev/null
    afa_out=/dev/null
fi

if ! error=$(perl "$trim_pl" "$result_stk" "$motif.tmp" "$fasta_out" "$afa_out" 2>&1)
then
    echo "trimAlignment.pl failed for $cluster${error:+: ${error//$'\n'/ }}" >&2
    rm -f "$motif.tmp" "$fasta.tmp" "$afa.tmp"
    exit 1
fi

# A usable motif has a consensus structure line and at least one sequence row,
# i.e. a line that starts with an identifier rather than # or / or a space.
if ! grep -q '#=GC SS_cons' "$motif.tmp" ||
    ! grep -qE '^[^#/[:space:]]' "$motif.tmp"
then
    mv "$motif.tmp" "$motif.invalid"
    rm -f "$fasta.tmp" "$afa.tmp"

    echo \
        "Empty motif for $cluster:" \
        "trimAlignment.pl produced no sequences or no SS_cons" \
        >&2

    exit 1
fi

mv "$motif.tmp" "$motif"

[[ "$optional_outputs" == yes ]] || exit 0

mv "$fasta.tmp" "$fasta"
mv "$afa.tmp" "$afa"

# The clustal view is for reading only; nothing downstream needs it, so a
# failure here is reported but does not fail the cluster.
apptainer exec \
    --bind "$cluster_dir:/work" \
    --pwd /work \
    "$rnatools" \
    esl-reformat \
    --informat stockholm \
    -o "${cluster}_motif.aln" \
    clustal \
    "${cluster}_motif.sto" \
    > /dev/null 2>&1 ||
    echo "esl-reformat failed for $cluster, .aln not written" >&2
