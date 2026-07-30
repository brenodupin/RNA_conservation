#!/usr/bin/env bash
#
# worker_step_03.sh — process ONE cluster for step 3.
#
#   worker_step_03.sh align <clusterName>    clustal omega, on the host
#   worker_step_03.sh fold  <clusterName>    RNALalifold, in its own container
#
# One cluster per invocation, matching upstream: clustalo_array_job.sh and
# rnalalifold_array_job.sh each ran as a SLURM array task handling a single
# cluster. 03_screen.sh runs many of these concurrently via xargs -P.
#
# This is a separate script rather than a function exported into xargs.
# `export -f` ships code through the environment, requires every variable the
# function touches to be exported by hand, and does NOT carry `set -euo
# pipefail` into the child -- so one forgotten export silently becomes an empty
# string in a subshell with no nounset to catch it. Sourcing info.sh here gets
# the whole configuration with nothing exported.
#
# Exits non-zero on failure; 03_screen.sh derives the failure list by checking
# which clusters lack output afterwards, so nothing is written to a shared file.
#
# Debug a single cluster by hand:
#   ./worker_step_03.sh fold NC_008590.1_220_38500_38750

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

mode="${1:-}"
name="${2:-}"
[ -n "$mode" ] && [ -n "$name" ] \
  || die "usage: worker_step_03.sh <align|fold> <clusterName>"

OUTDIR="$RNAC_DATA/$RNAC_STEP_03"
SPLITDIR="$RNAC_DATA/$RNAC_STEP_02/splits"
dir="$OUTDIR/$name"

case "$mode" in

  align)
    src="$SPLITDIR/${name}_cluster.fasta"
    aln="$dir/${name}_aligned.aln"
    dm="$dir/${name}_distMat.csv"

    [ -s "$src" ] || { warn "no cluster FASTA: $name"; exit 1; }

    mkdir -p "$dir"
    cp -f "$src" "$dir/"

    # clustalo refuses to overwrite an existing output file and exits non-zero,
    # so clear any leftover target -- including a .tmp from an interrupted run.
    rm -f "$aln" "$aln.tmp" "$dm" "$dm.tmp"

    # clustalo uses every core by default; with RNAC_JOBS workers running at
    # once that oversubscribes badly, so pin each to one thread when parallel.
    # At RNAC_JOBS=1 no --threads is passed, identical to upstream.
    threads=()
    [ "$RNAC_JOBS" -gt 1 ] && threads=(--threads 1)

    # Capture stderr rather than discarding it: clustalo reports why it did
    # something, and throwing that away turns a clear diagnosis into a puzzle.
    # (2>&1 >/dev/null sends stderr to the capture and stdout to the bin.)
    if ! err=$(clustalo -i "$dir/${name}_cluster.fasta" \
           --percent-id --full "${threads[@]}" \
           --distmat-out "$dm.tmp" \
           -o "$aln.tmp" --outfmt clu 2>&1 >/dev/null)
    then
      rm -f "$aln.tmp" "$dm.tmp"
      warn "clustalo failed: $name${err:+ — ${err//$'\n'/ }}"
      exit 1
    fi

    # Rename only on success, so an interrupted run never leaves a truncated
    # file that a later run mistakes for finished work.
    mv -f "$aln.tmp" "$aln"

    # The distance matrix is OPTIONAL. clustalo emits
    #   "Have only two sequences: Will not calculate/print distance matrix"
    # and exits 0, so no file appears for any 2-member cluster -- typically the
    # majority of them. It is diagnostic output that nothing downstream reads,
    # so its absence is not a failure.
    if [ -f "$dm.tmp" ]; then
      mv -f "$dm.tmp" "$dm"
    fi
    ;;

  fold)
    aln="$dir/${name}_aligned.aln"
    out="$dir/${name}_RNALalifold.out"

    [ -s "$aln" ] || { warn "no alignment: $name"; exit 1; }

    # Upstream binds the cluster directory into the container and passes the
    # .aln through it. The docker equivalent is mounting that directory at /work
    # and making it the working directory, which also catches any PostScript
    # side files RNALalifold drops in the CWD.
    if ! docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
           -v "$dir":/work -w /work "$RNAC_IMAGE" \
           RNALalifold -T "$RNAC_FOLD_TEMP" --noLP "${name}_aligned.aln" \
           > "$out.tmp" 2>/dev/null
    then
      rm -f "$out.tmp"
      warn "RNALalifold failed: $name"
      exit 1
    fi

    mv -f "$out.tmp" "$out"
    ;;

  *)
    die "unknown mode '$mode' (expected: align | fold)"
    ;;
esac