#!/usr/bin/env bash
#
# 03_screen.sh — align each cluster with clustal omega, then screen for
#                consensus secondary structure with RNALalifold.
#
#   in : $RNAC_DATA/$RNAC_STEP_02/splits/<rep>_cluster.fasta
#        $RNAC_DATA/$RNAC_STEP_02/cluster_count.tsv
#   out: $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_cluster.fasta
#        $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_aligned.aln          clustalo, clustal format
#        $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_distMat.csv          percent identity matrix
#        $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_RNALalifold.out      dot-bracket structures
#        $RNAC_DATA/$RNAC_STEP_03/RNALalifold_passedList.txt       clusters worth folding
#
#   ./03_screen.sh                    normal run
#   RNAC_FORCE=1 ./03_screen.sh            realign and refold everything
#   RNAC_FOLD_TEMP=37 ./03_screen.sh       fold at 37 C instead of 21 C
#
# clustal omega runs on the host (it is NOT in the container); RNALalifold comes
# from the ViennaRNA package inside the container.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$RNAC_STEP_03"
INDIR="$RNAC_DATA/$RNAC_STEP_02"
SPLITDIR="$INDIR/splits"
COUNTS="$INDIR/cluster_count.tsv"
OUTDIR="$RNAC_DATA/$STEP"

start_log "$STEP"

require_dir "$SPLITDIR"
require_file "$COUNTS"
require_cmd clustalo
docker_preflight
mkdir -p "$OUTDIR"

log "fold temperature: ${RNAC_FOLD_TEMP} C"

# Upstream reads the cluster name from column 2 of cluster_count.tsv with
# `cut -f2`, which is the same tab dependency step 2 guards against.
mapfile -t clusters < <(cut -f2 "$COUNTS")
total=${#clusters[@]}
[ "$total" -gt 0 ] || die "no clusters listed in $COUNTS"
log "$total cluster(s) to screen"

# --- phase 1: clustal omega ------------------------------------------------
#
# Exact option set from clustalo_array_job.sh. --percent-id requires both
# --full and --distmat-out.
#
# clustalo refuses to overwrite an existing output file and exits non-zero, so
# on RNAC_FORCE=1 the previous outputs are removed rather than passing --force. That
# keeps the invocation identical to upstream's.

aligned=0 skipped_aln=0 missing=0
i=0

for name in "${clusters[@]}"; do
  i=$(( i + 1 ))
  src="$SPLITDIR/${name}_cluster.fasta"
  dir="$OUTDIR/$name"
  aln="$dir/${name}_aligned.aln"

  if [ ! -s "$src" ]; then
    warn "[$i/$total] $name: no cluster FASTA in splits/ — skipping"
    missing=$(( missing + 1 ))
    continue
  fi

  if [ -s "$aln" ] && [ "$RNAC_FORCE" != "1" ]; then
    skipped_aln=$(( skipped_aln + 1 ))
    continue
  fi

  mkdir -p "$dir"
  cp -f "$src" "$dir/"

  rm -f "$aln" "$dir/${name}_distMat.csv"
  log "[$i/$total] clustalo: $name"
  clustalo -i "$dir/${name}_cluster.fasta" \
    --percent-id --full \
    --distmat-out "$dir/${name}_distMat.csv" \
    -o "$aln" --outfmt clu

  require_file "$aln"
  aligned=$(( aligned + 1 ))
done

log "alignments: $aligned new, $skipped_aln reused, $missing missing input"

# --- phase 2: RNALalifold --------------------------------------------------
#
# Upstream binds the cluster directory into the container and passes the .aln
# through it. With docker the equivalent is mounting that directory at /work and
# setting it as the working directory, which also catches any PostScript side
# files RNALalifold may drop in the CWD.
#
# -u keeps outputs owned by you rather than root; -e HOME=/tmp gives the
# unmapped UID a writable home. The redirect happens on the host.

folded=0 skipped_fold=0
i=0

for name in "${clusters[@]}"; do
  i=$(( i + 1 ))
  dir="$OUTDIR/$name"
  aln="$dir/${name}_aligned.aln"
  out="$dir/${name}_RNALalifold.out"

  [ -s "$aln" ] || continue

  if [ -s "$out" ] && [ "$RNAC_FORCE" != "1" ]; then
    skipped_fold=$(( skipped_fold + 1 ))
    continue
  fi

  log "[$i/$total] RNALalifold: $name"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$dir":/work -w /work "$RNAC_IMAGE" \
    RNALalifold -T "$RNAC_FOLD_TEMP" --noLP "${name}_aligned.aln" \
    > "$out"

  folded=$(( folded + 1 ))
done

log "RNALalifold: $folded new, $skipped_fold reused"

# --- passed list -----------------------------------------------------------
#
# Upstream:
#   for i in `find -iname "*_RNALalifold.out"`; do grep -m1 -H "((" $i; done \
#     | cut -d':' -f1 | cut -d'/' -f3 | sed -e 's/_RNALalifold.out//g'
#
# `cut -d'/' -f3` hardcodes their directory depth; basename is equivalent and
# depth-independent.
#
# The test is the literal substring "((" -- two adjacent opening brackets, i.e.
# at least one stacked pair. The README describes this as "2 base-pairs", but
# note it is stricter than that: a structure like "(.(...).)" has two base pairs
# and no "((" . In practice --noLP already forbids lone pairs, so any structure
# RNALalifold reports here will contain a stack.

PASSED="$OUTDIR/$RNAC_PASSED_LIST"

shopt -s nullglob
outs=( "$OUTDIR"/*/*_RNALalifold.out )
shopt -u nullglob

: > "$PASSED"
if [ ${#outs[@]} -gt 0 ]; then
  grep -lF '((' "${outs[@]}" 2>/dev/null \
    | xargs -r -n1 basename \
    | sed 's/_RNALalifold\.out$//' \
    | LC_ALL=C sort > "$PASSED"
fi

n_screened=${#outs[@]}
n_passed=$(wc -l < "$PASSED")

# --- report ----------------------------------------------------------------

echo
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "listed in cluster_count.tsv"        "$total"
printf '%-42s %10s\n' "aligned by clustalo"                "$(( total - missing ))"
printf '%-42s %10s\n' "screened by RNALalifold"            "$n_screened"
printf '%-42s %10s\n' "with a predicted stem (passed)"     "$n_passed"
printf '%-42s %10s\n' "no consensus structure (dropped)"   "$(( n_screened - n_passed ))"

echo
if [ "$n_passed" -eq 0 ]; then
  warn "no cluster showed consensus structure potential"
  warn "this is a real result, not an error — step 4 would have nothing to fold"
else
  log "step 4 input: $PASSED ($n_passed clusters)"
fi
log "outputs in $OUTDIR"