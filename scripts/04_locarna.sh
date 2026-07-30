#!/usr/bin/env bash
#
# 04_locarna.sh — predict a common secondary structure for each cluster that
#                 passed step 3, using mlocarna (LocARNA-P).
#
#   in : $RNAC_DATA/$RNAC_STEP_03/RNALalifold_passedList.txt   clusters worth folding
#        $RNAC_DATA/$RNAC_STEP_02/splits/<rep>_cluster.fasta   unaligned sequences
#        $RNAC_DATA/$RNAC_STEP_02/cluster_count.tsv            members per cluster (scheduling only)
#   out: $RNAC_DATA/$RNAC_STEP_04/<rep>/<rep>_locarnap/        full mlocarna tree
#        $RNAC_DATA/$RNAC_STEP_04/<rep>/<rep>_result.stk       the conserved structure
#        $RNAC_DATA/$RNAC_STEP_04/<rep>/<rep>_mlocarna.log     --moreverbose output
#        $RNAC_DATA/$RNAC_STEP_04/timings.tsv                  per-cluster wall time
#        $RNAC_DATA/$RNAC_STEP_04/failed_clusters.txt          anything that errored
#
#   ./04_locarna.sh                          normal run (and resumes an interrupted one)
#   RNAC_JOBS=10 ./04_locarna.sh             10 clusters at a time
#   RNAC_LOCARNA_THREADS=1 ./04_locarna.sh   1 thread each, so more clusters fit
#   RNAC_FORCE=1 ./04_locarna.sh             refold everything from scratch
#
# Like step 3, this script only orchestrates: it decides which clusters still
# need work, runs worker_step_04.sh once per cluster via xargs -P, and verifies
# the results afterwards. The per-cluster work lives in worker_step_04.sh.
#
# This is by far the longest step in the pipeline. mlocarna cost grows roughly
# with (members^2 x length), so a handful of large clusters can outlast the
# other few thousand put together -- which is why work is scheduled largest
# first (see RNAC_LOCARNA_ORDER) and why every cluster is timed.
#
# Everything runs in the container; nothing here is needed on the host.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$RNAC_STEP_04"
SPLITDIR="$RNAC_DATA/$RNAC_STEP_02/splits"
COUNTS="$RNAC_DATA/$RNAC_STEP_02/cluster_count.tsv"
PASSED="$RNAC_DATA/$RNAC_STEP_03/$RNAC_PASSED_LIST"
OUTDIR="$RNAC_DATA/$STEP"
WORKER="$HERE/worker_step_04.sh"

start_log "$STEP"

require_dir "$SPLITDIR"
[ -f "$PASSED" ] \
  || die "no passed list at $PASSED — run 03_screen.sh first"
[ -x "$WORKER" ] || die "worker script is not executable: $WORKER  (fix with: chmod +x \"$WORKER\")"
require_cmd xargs
docker_preflight
mkdir -p "$OUTDIR"

[[ "$RNAC_JOBS" =~ ^[0-9]+$ ]] && [ "$RNAC_JOBS" -ge 1 ] \
  || die "RNAC_JOBS must be a positive integer (got '$RNAC_JOBS')"
[[ "$RNAC_LOCARNA_THREADS" =~ ^[0-9]+$ ]] && [ "$RNAC_LOCARNA_THREADS" -ge 1 ] \
  || die "RNAC_LOCARNA_THREADS must be a positive integer (got '$RNAC_LOCARNA_THREADS')"

# An empty passed list is a legitimate result from step 3, not a failure.
mapfile -t clusters < <(grep -v '^[[:space:]]*$' "$PASSED" || true)
total=${#clusters[@]}
if [ "$total" -eq 0 ]; then
  warn "$PASSED is empty — step 3 found no cluster with consensus structure potential"
  warn "nothing to fold; this is a real result, not an error"
  exit 0
fi

cores=$(( RNAC_JOBS * RNAC_LOCARNA_THREADS ))
log "$total cluster(s) passed screening"
log "fold temperature: ${RNAC_FOLD_TEMP} C   parallel jobs: $RNAC_JOBS x ${RNAC_LOCARNA_THREADS} thread(s) = $cores core(s)"
[ -n "$RNAC_LOCARNA_MEM" ] && log "per-container memory cap: $RNAC_LOCARNA_MEM"
[ "$cores" -le "$(nproc)" ] \
  || warn "requesting $cores cores on a $(nproc)-core host — expect contention"

# --- work list -------------------------------------------------------------
#
# A cluster is finished when <name>_result.stk exists. The worker writes that
# file last, by copying it out of the mlocarna tgtdir only after mlocarna has
# exited 0, so its presence means the whole run completed -- unlike the tgtdir
# itself, which appears immediately and is populated throughout.

todo=()
missing=0
for name in "${clusters[@]}"; do
  if [ ! -s "$SPLITDIR/${name}_cluster.fasta" ]; then
    missing=$(( missing + 1 ))
    continue
  fi
  [ -s "$OUTDIR/$name/${name}_result.stk" ] && [ "$RNAC_FORCE" != "1" ] && continue
  todo+=("$name")
done

[ "$missing" -eq 0 ] || warn "$missing passed cluster(s) have no FASTA in splits/"
expected=$(( total - missing ))
[ "$expected" -gt 0 ] || die "none of the passed clusters have a FASTA in $SPLITDIR"

# --- scheduling ------------------------------------------------------------
#
# With runtimes spanning orders of magnitude, order matters: start the largest
# cluster last and every other worker sits idle waiting for it. Largest first
# keeps the tail busy, and costs nothing, since clusters are independent and the
# output is per-cluster.
#
# Set RNAC_LOCARNA_ORDER=list to process the passed list in its own order
# instead (upstream's behaviour -- useful only for lockstep comparison).

if [ "$RNAC_LOCARNA_ORDER" = "size" ] && [ -s "$COUNTS" ] && [ ${#todo[@]} -gt 1 ]; then
  # cluster_count.tsv is <count>\t<name>; join on name, sort numerically
  # descending, keep unlisted clusters (count 0) at the end rather than dropping
  # them -- the count is a scheduling hint, not a gate.
  mapfile -t todo < <(
    awk -F'\t' '
      NR==FNR { n[$2] = $1; next }
      { printf "%d\t%s\n", ($1 in n ? n[$1] : 0), $1 }
    ' "$COUNTS" <(printf '%s\n' "${todo[@]}") \
    | LC_ALL=C sort -k1,1nr -k2,2 \
    | cut -f2
  )
  log "scheduling largest cluster first"
fi

# --- run -------------------------------------------------------------------

if [ ${#todo[@]} -eq 0 ]; then
  log "all $expected cluster(s) already folded"
else
  log "folding: ${#todo[@]} cluster(s), $RNAC_JOBS at a time"
  log "this is the slow step — follow along with: tail -f $LOGFILE"

  progress_start "folded" "$expected" \
    "$OUTDIR" -mindepth 2 -maxdepth 2 -name '*_result.stk'

  rc=0
  printf '%s\n' "${todo[@]}" \
    | xargs -r -P "$RNAC_JOBS" -n1 "$WORKER" || rc=$?

  progress_stop

  # xargs reports a status, not a count: 123 means "at least one invocation
  # exited 1-125", which is the ordinary case when some clusters fail. The real
  # tally comes from the failure list below. 126/127 mean the worker itself
  # could not be run, which is a setup problem rather than a data problem and
  # would otherwise be misread as "every cluster failed".
  case "$rc" in
    0)   ;;
    123) log "some clusters failed — counted below" ;;
    125) warn "xargs itself failed" ;;
    126) die "worker is not executable: $WORKER" ;;
    127) die "worker not found: $WORKER" ;;
    *)   warn "xargs exited $rc" ;;
  esac
fi

# --- failures --------------------------------------------------------------
#
# Derived from which clusters lack output, so workers never append to a shared
# file from RNAC_JOBS processes at once.

FAILED="$OUTDIR/failed_clusters.txt"
{
  for name in "${clusters[@]}"; do
    [ -s "$SPLITDIR/${name}_cluster.fasta" ] || continue
    [ -s "$OUTDIR/$name/${name}_result.stk" ] || printf '%s\n' "$name"
  done
} | LC_ALL=C sort -u > "$FAILED.tmp"
mv -f "$FAILED.tmp" "$FAILED"
n_failed=$(wc -l < "$FAILED")

# --- timings ---------------------------------------------------------------
#
# Collated from the per-cluster files the workers wrote. Rebuilt from scratch on
# every run so it always describes what is currently on disk, and sorted by
# duration so the clusters worth worrying about are at the top.

TIMINGS="$OUTDIR/timings.tsv"
{
  printf 'cluster\tn_seqs\tstart\tend\tduration_sec\texit_code\n'
  find "$OUTDIR" -mindepth 2 -maxdepth 2 -type f -name '*_timing.tsv' -print0 \
    | sort -z \
    | xargs -0 -r cat \
    | LC_ALL=C sort -t$'\t' -k5,5nr
} > "$TIMINGS.tmp"
mv -f "$TIMINGS.tmp" "$TIMINGS"

n_done=$(find "$OUTDIR" -mindepth 2 -maxdepth 2 -type f -name '*_result.stk' | wc -l)

# --- report ----------------------------------------------------------------

echo
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "passed step 3 screening"            "$total"
printf '%-42s %10s\n' "missing from splits/"               "$missing"
printf '%-42s %10s\n' "attempted this run"                 "${#todo[@]}"
printf '%-42s %10s\n' "with a predicted structure"         "$n_done"
printf '%-42s %10s\n' "failed"                             "$n_failed"

# Slowest few, since that is what determines whether the next run needs
# different RNAC_JOBS / RNAC_LOCARNA_THREADS settings.
if [ "$(wc -l < "$TIMINGS")" -gt 1 ]; then
  echo
  printf 'slowest clusters this run:\n'
  awk -F'\t' 'NR>1 && NR<=6 { printf "  %-44s %5s seqs %8s s\n", $1, $2, $5 }' "$TIMINGS"
fi

echo
[ "$n_failed" -eq 0 ] || warn "$n_failed cluster(s) failed — see $FAILED, then re-run to retry just those"
if [ "$n_done" -eq 0 ]; then
  die "no cluster produced a structure — check $OUTDIR/*/*_mlocarna.log"
fi
log "step 5 input: $OUTDIR/<cluster>/<cluster>_result.stk ($n_done clusters)"
log "timings: $TIMINGS"
log "outputs in $OUTDIR"