#!/bin/bash

# Screen one cluster alignment using RNALalifold.
#
# Arguments:
#   container, cluster directory, alignment name, output, temperature

set -euo pipefail

rnatools=$1
cluster_dir=$2
alignment=$3
output=$4
fold_temperature=$5

apptainer exec \
    --bind "$cluster_dir:/work" \
    --pwd /work \
    "$rnatools" \
    RNALalifold \
    -T "$fold_temperature" \
    --noLP \
    "$alignment" \
    > "$output" \
    2> /dev/null