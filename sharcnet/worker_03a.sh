#!/bin/bash

# Align one cluster using Clustal Omega.
#
# Arguments:
#   source FASTA, copied FASTA, alignment output, distance-matrix output

set -euo pipefail

source_fasta=$1
cluster_fasta=$2
alignment=$3
distance_matrix=$4

mkdir -p "$(dirname "$cluster_fasta")"
cp "$source_fasta" "$cluster_fasta"

clustalo \
    -i "$cluster_fasta" \
    --threads=1 \
    --percent-id \
    --full \
    --distmat-out="$distance_matrix" \
    --outfmt=clu \
    --force \
    -o "$alignment" \
    > /dev/null