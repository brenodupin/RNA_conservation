#!/bin/bash

# Render one Stockholm alignment as a coloured HTML page with
# stockholm_to_html.pl, in the rnatools container. For reading only: nothing
# in the pipeline uses the page.
#
# Arguments:
#   container, directory, alignment, HTML output (both relative to the
#   directory), directory of the riboswitch helper scripts inside the container

set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 CONTAINER DIR ALIGNMENT HTML SCRIPTS_DIR" >&2
    exit 1
fi

rnatools=$1
dir=$2
alignment=$3
html=$4
scripts_dir=$5

[[ -s "$dir/$alignment" ]] || {
    echo "Alignment not found: $dir/$alignment" >&2
    exit 1
}

rm -f "$dir/$html" "$dir/$html.tmp"

# BioPerl's exceptions run to a dozen lines; the first one says what failed.
if ! error=$(
    apptainer exec \
        --bind "$dir:/work" \
        --pwd /work \
        "$rnatools" \
        perl "$scripts_dir/stockholm_to_html.pl" "$alignment" "$html.tmp" \
        2>&1 > /dev/null
)
then
    echo "stockholm_to_html.pl failed for $dir/$alignment: $(head -n 1 <<< "$error")" >&2
    rm -f "$dir/$html.tmp"
    exit 1
fi

# It can also die without a useful exit status; no page means it failed.
if [[ ! -s "$dir/$html.tmp" ]]; then
    echo "stockholm_to_html.pl wrote nothing for $dir/$alignment: $(head -n 1 <<< "$error")" >&2
    rm -f "$dir/$html.tmp"
    exit 1
fi

mv "$dir/$html.tmp" "$dir/$html"
