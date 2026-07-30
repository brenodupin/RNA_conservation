#!/usr/bin/env bash
#
# 03_screen.sh — align each cluster with clustal omega, then screen for
#                consensus secondary structure with RNALalifold.
#
#   in : $RNAC_DATA/$RNAC_STEP_02/splits/<rep>_cluster.fasta
#        $RNAC_DATA/$RNAC_STEP_02/cluster_count.tsv
#   out: $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_aligned.aln      clustalo, clustal format
#        $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_distMat.csv      percent identity matrix
#        $RNAC_DATA/$RNAC_STEP_03/<rep>/<rep>_RNALalifold.out  dot-bracket structures
#        $RNAC_DATA/$RNAC_STEP_03/RNALalifold_passedList.txt   clusters worth folding
#        $RNAC_DATA/$RNAC_STEP_03/failed_clusters.txt          anything that errored
#
#   ./03_screen.sh                     normal run
#   RNAC_JOBS=8 ./03_screen.sh         8 clusters at a time
#   RNAC_FORCE=1 ./03_screen.sh        realign and refold everything
#
# This script only orchestrates: it decides which clusters still need work, runs
# worker_step_03.sh once per cluster via xargs -P, and verifies the results. The
# per-cluster work lives in worker_step_03.sh.
#
# One worker invocation per cluster, matching upstream's SLURM array tasks --
# deliberately not batched, so each RNALalifold runs in its own container
# exactly as the published method did.
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
WORKER="$HERE/worker_step_03.sh"

start_log "$STEP"

require_dir "$SPLITDIR"
require_file "$COUNTS"
[ -x "$WORKER" ] || die "worker script is not executable: $WORKER  (fix with: chmod +x \"$WORKER\")"
require_cmd clustalo xargs
docker_preflight
mkdir -p "$OUTDIR"

[[ "$RNAC_JOBS" =~ ^[0-9]+$ ]] && [ "$RNAC_JOBS" -ge 1 ] \
  || die "RNAC_JOBS must be a positive integer (got '$RNAC_JOBS')"

log "fold temperature: ${RNAC_FOLD_TEMP} C   parallel jobs: $RNAC_JOBS"

# Upstream reads the cluster name from column 2 of cluster_count.tsv with
# `cut -f2`, the same tab dependency step 2 guards against.
mapfile -t clusters < <(cut -f2 "$COUNTS")
total=${#clusters[@]}
[ "$total" -gt 0 ] || die "no clusters listed in $COUNTS"
log "$total cluster(s) listed"

# --- helpers ---------------------------------------------------------------

count_out() {  # count_out <glob>
  find "$OUTDIR" -mindepth 2 -maxdepth 2 -type f -name "$1" | wc -l
}

# run_phase <mode> <label> <target-count> <result-glob> <names...>
# Runs 03_worker.sh once per cluster, RNAC_JOBS at a time.
run_phase() {
  local mode=$1 label=$2 target=$3 glob=$4; shift 4

  log "$label: $# cluster(s), $RNAC_JOBS at a time"

  progress_start "$label" "$target" \
    "$OUTDIR" -mindepth 2 -maxdepth 2 -name "$glob"

  local rc=0
  printf '%s\n' "$@" \
    | xargs -r -P "$RNAC_JOBS" -n1 "$WORKER" "$mode" || rc=$?

  progress_stop
  [ "$rc" -eq 0 ] || warn "$label: $rc worker invocation(s) reported failure"
}

# failures_among <result-file-suffix> <names...>
# A cluster that produced no output failed. Deriving the list this way avoids
# workers appending to a shared file from RNAC_JOBS processes at once.
failures_among() {
  local suffix=$1; shift
  local n
  for n in "$@"; do
    [ -s "$OUTDIR/$n/${n}${suffix}" ] || printf '%s\n' "$n"
  done
}

# --- phase 1: clustal omega ------------------------------------------------

todo_aln=()
missing=0
for name in "${clusters[@]}"; do
  if [ ! -s "$SPLITDIR/${name}_cluster.fasta" ]; then
    missing=$(( missing + 1 ))
    continue
  fi
  [ -s "$OUTDIR/$name/${name}_aligned.aln" ] && [ "$RNAC_FORCE" != "1" ] && continue
  todo_aln+=("$name")
done

[ "$missing" -eq 0 ] || warn "$missing cluster(s) in cluster_count.tsv have no FASTA in splits/"
expected_aln=$(( total - missing ))

if [ ${#todo_aln[@]} -eq 0 ]; then
  log "all $expected_aln alignment(s) already present"
else
  run_phase align "aligned" "$expected_aln" '*_aligned.aln' "${todo_aln[@]}"
fi

n_aln=$(count_out '*_aligned.aln')
[ "$n_aln" -eq "$expected_aln" ] \
  || die "only $n_aln of $expected_aln alignments present — re-run to finish before screening"
log "alignments complete: $n_aln"

# --- phase 2: RNALalifold --------------------------------------------------

todo_fold=()
for name in "${clusters[@]}"; do
  [ -s "$OUTDIR/$name/${name}_aligned.aln" ] || continue
  [ -s "$OUTDIR/$name/${name}_RNALalifold.out" ] && [ "$RNAC_FORCE" != "1" ] && continue
  todo_fold+=("$name")
done

if [ ${#todo_fold[@]} -eq 0 ]; then
  log "all $n_aln structure prediction(s) already present"
else
  run_phase fold "folded" "$n_aln" '*_RNALalifold.out' "${todo_fold[@]}"
fi

# --- collect failures ------------------------------------------------------

FAILED="$OUTDIR/failed_clusters.txt"
{
  [ ${#todo_aln[@]}  -gt 0 ] && failures_among '_aligned.aln'      "${todo_aln[@]}"
  [ ${#todo_fold[@]} -gt 0 ] && failures_among '_RNALalifold.out'  "${todo_fold[@]}"
  true
} | LC_ALL=C sort -u > "$FAILED.tmp"
mv -f "$FAILED.tmp" "$FAILED"
n_failed=$(wc -l < "$FAILED")
[ "$n_failed" -eq 0 ] || warn "$n_failed cluster(s) failed — see $FAILED"

# --- completeness gate -----------------------------------------------------
#
# The passed list is derived from whichever *_RNALalifold.out files exist, so it
# must not be built from a partial screen: a cluster with no .out file is
# indistinguishable from one that legitimately found no structure, and would be
# silently dropped from step 4.

n_out=$(count_out '*_RNALalifold.out')
[ "$n_out" -eq "$n_aln" ] \
  || die "only $n_out of $n_aln clusters screened — re-run to finish before the passed list is built"
log "structure predictions complete: $n_out"

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
# at least one stacked pair. The README calls this "2 base-pairs", but it is
# stricter: "(.(...).)" has two base pairs and no "((" . In practice --noLP
# already forbids lone pairs, so any structure reported here contains a stack.

PASSED="$OUTDIR/$RNAC_PASSED_LIST"

find "$OUTDIR" -mindepth 2 -maxdepth 2 -type f -name '*_RNALalifold.out' -print0 \
  | xargs -0 -r grep -lF '((' \
  | xargs -r -n1 basename \
  | sed 's/_RNALalifold\.out$//' \
  | LC_ALL=C sort > "$PASSED.tmp"
mv -f "$PASSED.tmp" "$PASSED"

n_passed=$(wc -l < "$PASSED")

# --- report ----------------------------------------------------------------

echo
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "listed in cluster_count.tsv"        "$total"
printf '%-42s %10s\n' "missing from splits/"               "$missing"
printf '%-42s %10s\n' "aligned by clustalo"                "$n_aln"
printf '%-42s %10s\n' "screened by RNALalifold"            "$n_out"
printf '%-42s %10s\n' "failed"                             "$n_failed"
printf '%-42s %10s\n' "with a predicted stem (passed)"     "$n_passed"
printf '%-42s %10s\n' "no consensus structure (dropped)"   "$(( n_out - n_passed ))"

echo
if [ "$n_passed" -eq 0 ]; then
  warn "no cluster showed consensus structure potential"
  warn "this is a real result, not an error — step 4 would have nothing to fold"
else
  log "step 4 input: $PASSED ($n_passed clusters)"
fi
log "outputs in $OUTDIR"