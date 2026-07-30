#!/usr/bin/env bash
#
# 05_eval.sh — evaluate each folded cluster: trim the LocARNA alignment to the
#              motif, score it with RNA-SCoRE, then test High/Mid clusters for
#              covariation with R-scape.
#
#   in : $RNAC_DATA/$RNAC_STEP_04/<rep>/<rep>_result.stk        LocARNA structure
#        $RNAC_REPO/step5_evaluationOfRNAstructures/*.pl        trim + RNA-SCoRE
#   out: $RNAC_DATA/$RNAC_STEP_05/<rep>/<rep>_motif.sto|.fasta|.afa   trimmed motif
#        $RNAC_DATA/$RNAC_STEP_05/<rep>/<rep>_motif_cleaned.sto       passed seqs (High/Mid)
#        $RNAC_DATA/$RNAC_STEP_05/<rep>/<rep>_motif_evaluated.tsv     per-sequence verdicts
#        $RNAC_DATA/$RNAC_STEP_05/<rep>/<rep>_rank.tsv                this cluster's rank line
#        $RNAC_DATA/$RNAC_STEP_05/<rep>/<rep>_<label>_rscape_out/     R-scape output (High/Mid)
#        $RNAC_DATA/$RNAC_STEP_05/allClusters_evaluation.tsv          all rank lines, collated
#        $RNAC_DATA/$RNAC_STEP_05/passed_High.txt                     ranked High
#        $RNAC_DATA/$RNAC_STEP_05/passed_Mid.txt                      ranked Mid
#        $RNAC_DATA/$RNAC_STEP_05/failed_clusters.txt                 errored in any phase
#
#   ./05_eval.sh                          normal run (and resumes an interrupted one)
#   RNAC_JOBS=10 ./05_eval.sh             10 clusters at a time
#   RNAC_SCORE_MIN_RANK=High ./05_eval.sh only run R-scape on High clusters
#   RNAC_FORCE=1 ./05_eval.sh             redo trim, eval and R-scape from scratch
#
# Structured like 03_screen.sh: this script only orchestrates. It runs
# worker_step_05.sh once per cluster via xargs -P, in three phases, and gates
# between them --
#
#   trim   : every folded cluster
#   eval   : every cluster that trimmed to a motif
#   rscape : only clusters RNA-SCoRE ranked >= RNAC_SCORE_MIN_RANK
#
# The gate is the reason the phases cannot be collapsed into one worker call:
# the rscape work list is not known until eval has ranked everything.
#
# trim and eval run on the host (Perl); esl-reformat and R-scape run in the
# container. clustalo is NOT needed here.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$RNAC_STEP_05"
LOCARNA_DIR="$RNAC_DATA/$RNAC_STEP_04"
OUTDIR="$RNAC_DATA/$STEP"
WORKER="$HERE/worker_step_05.sh"
STEPSRC="$RNAC_REPO/step5_evaluationOfRNAstructures"

start_log "$STEP"

require_dir "$LOCARNA_DIR"
[ -x "$WORKER" ] || die "worker script is not executable: $WORKER  (fix with: chmod +x \"$WORKER\")"
require_cmd perl xargs
require_file "$STEPSRC/trimAlignment.pl"
require_file "$STEPSRC/RNA-SCoRE.pl" \
  || die "RNA-SCoRE.pl not found in $STEPSRC — clone it from RodrigoReisLab/RNA-SCoRE"
docker_preflight                       # esl-reformat + R-scape run in the container
mkdir -p "$OUTDIR"

[[ "$RNAC_JOBS" =~ ^[0-9]+$ ]] && [ "$RNAC_JOBS" -ge 1 ] \
  || die "RNAC_JOBS must be a positive integer (got '$RNAC_JOBS')"
case "$RNAC_SCORE_MIN_RANK" in
  High|Mid) ;;
  *) die "RNAC_SCORE_MIN_RANK must be High or Mid (got '$RNAC_SCORE_MIN_RANK')" ;;
esac

# Clusters to consider: those step 4 actually folded (a real _result.stk).
mapfile -t clusters < <(
  find "$LOCARNA_DIR" -mindepth 2 -maxdepth 2 -type f -name '*_result.stk' \
    | sed -E 's#.*/([^/]+)/[^/]+_result\.stk$#\1#' \
    | LC_ALL=C sort
)
total=${#clusters[@]}
[ "$total" -gt 0 ] \
  || die "no *_result.stk under $LOCARNA_DIR — run 04_locarna.sh first"
log "$total folded cluster(s) to evaluate   parallel jobs: $RNAC_JOBS"
log "RNA-SCoRE thresholds: mt=$RNAC_SCORE_MT t=$RNAC_SCORE_BP gc=$RNAC_SCORE_GC dupl=$RNAC_SCORE_DUPL"
log "R-scape on clusters ranked >= $RNAC_SCORE_MIN_RANK"

# --- generic phase runner (mirrors 03_screen.sh) ---------------------------
#
# run_phase <mode> <label> <target> <marker-glob> <names...>
# Runs the worker once per name, RNAC_JOBS at a time, with a filesystem
# progress meter counting <marker-glob> under the output tree.
run_phase() {
  local mode=$1 label=$2 target=$3 glob=$4; shift 4
  [ "$#" -gt 0 ] || { log "$label: nothing to do"; return 0; }

  log "$label: $# cluster(s), $RNAC_JOBS at a time"
  progress_start "$label" "$target" \
    "$OUTDIR" -mindepth 2 -maxdepth 2 -name "$glob"

  local rc=0
  printf '%s\n' "$@" \
    | xargs -r -P "$RNAC_JOBS" -n1 "$WORKER" "$mode" || rc=$?

  progress_stop
  case "$rc" in
    0|123) ;;                                   # 123 = some clusters failed; counted later
    126) die "worker not executable: $WORKER" ;;
    127) die "worker not found: $WORKER" ;;
    *)   warn "$label: xargs exited $rc" ;;
  esac
}

# has_marker <name> <relative-file>   — is this cluster's phase output present?
has_marker() { [ -s "$OUTDIR/$1/$2" ]; }

# ===========================================================================
# phase 1 — trim
# ===========================================================================

todo=()
for name in "${clusters[@]}"; do
  { has_marker "$name" "${name}_motif.sto" && [ "$RNAC_FORCE" != "1" ]; } && continue
  todo+=("$name")
done
run_phase trim "trimmed" "$total" '*_motif.sto' "${todo[@]}"

mapfile -t trimmed < <(
  for name in "${clusters[@]}"; do
    has_marker "$name" "${name}_motif.sto" && printf '%s\n' "$name"
  done
)
n_trim=${#trimmed[@]}
log "trimmed to a motif: $n_trim / $total"
[ "$n_trim" -gt 0 ] || die "no cluster produced a trimmed motif — check $OUTDIR/*/result.stk headers"

# ===========================================================================
# phase 2 — RNA-SCoRE  (marker: <name>_rank.tsv, always written on success)
# ===========================================================================

todo=()
for name in "${trimmed[@]}"; do
  { has_marker "$name" "${name}_rank.tsv" && [ "$RNAC_FORCE" != "1" ]; } && continue
  todo+=("$name")
done
run_phase eval "scored" "$n_trim" '*_rank.tsv' "${todo[@]}"

mapfile -t scored < <(
  for name in "${trimmed[@]}"; do
    has_marker "$name" "${name}_rank.tsv" && printf '%s\n' "$name"
  done
)
n_scored=${#scored[@]}
log "scored by RNA-SCoRE: $n_scored / $n_trim"

# --- collate ranks ---------------------------------------------------------
#
# One allClusters_evaluation.tsv assembled from the per-cluster rank files,
# rebuilt every run so it always matches what is on disk. Same reasoning as
# step 4's timings.tsv: no concurrent appends, nothing to interleave.
#
# The rank is the last tab-separated field of each line: High, Mid, or Low
# (or "Low - No evaluation ..."). Split by rank into passed_High/Mid lists that
# drive the R-scape phase and feed step 6.

EVALALL="$OUTDIR/allClusters_evaluation.tsv"
find "$OUTDIR" -mindepth 2 -maxdepth 2 -type f -name '*_rank.tsv' -print0 \
  | sort -z | xargs -0 -r cat \
  | LC_ALL=C sort > "$EVALALL.tmp"
mv -f "$EVALALL.tmp" "$EVALALL"

# rank = first whitespace-token of the final tab-separated field
rank_of() { awk -F'\t' '{split($NF,a," "); print a[1]}'; }

HIGH="$OUTDIR/passed_High.txt"
MID="$OUTDIR/passed_Mid.txt"
awk -F'\t' '{split($NF,a," "); if(a[1]=="High") print $1}' "$EVALALL" | LC_ALL=C sort > "$HIGH"
awk -F'\t' '{split($NF,a," "); if(a[1]=="Mid")  print $1}' "$EVALALL" | LC_ALL=C sort > "$MID"
n_high=$(wc -l < "$HIGH"); n_mid=$(wc -l < "$MID")
n_low=$(( n_scored - n_high - n_mid ))

# R-scape input list = High, plus Mid unless MIN_RANK=High.
mapfile -t rscape_pool < <(
  cat "$HIGH"
  [ "$RNAC_SCORE_MIN_RANK" = "Mid" ] && cat "$MID"
)
# drop possible blank line if a list was empty
mapfile -t rscape_pool < <(printf '%s\n' "${rscape_pool[@]}" | grep -v '^[[:space:]]*$' || true)

# ===========================================================================
# phase 3 — R-scape  (marker: <name>_rscape.done)
# ===========================================================================

todo=()
for name in "${rscape_pool[@]}"; do
  [ -s "$OUTDIR/$name/${name}_motif_cleaned.sto" ] || continue   # defensive
  { has_marker "$name" "${name}_rscape.done" && [ "$RNAC_FORCE" != "1" ]; } && continue
  todo+=("$name")
done
run_phase rscape "covariation" "${#rscape_pool[@]}" '*_rscape.done' "${todo[@]}"

n_rscape=$(find "$OUTDIR" -mindepth 2 -maxdepth 2 -type f -name '*_rscape.done' | wc -l)

# ===========================================================================
# failures + report
# ===========================================================================
#
# A cluster failed if it trimmed but never got a rank, or was promoted to
# R-scape but has no .done. Derived from missing markers so workers never share
# a file.

FAILED="$OUTDIR/failed_clusters.txt"
{
  for name in "${trimmed[@]}"; do
    has_marker "$name" "${name}_rank.tsv" || printf '%s\ttrim_or_eval\n' "$name"
  done
  for name in "${rscape_pool[@]}"; do
    has_marker "$name" "${name}_rscape.done" || printf '%s\trscape\n' "$name"
  done
} | LC_ALL=C sort -u > "$FAILED.tmp"
mv -f "$FAILED.tmp" "$FAILED"
n_failed=$(wc -l < "$FAILED")

echo
printf '%-42s %10s\n' STAGE CLUSTERS
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "folded (step 4 input)"              "$total"
printf '%-42s %10s\n' "trimmed to a motif"                 "$n_trim"
printf '%-42s %10s\n' "scored by RNA-SCoRE"                "$n_scored"
printf '%-42s %10s\n' "  ranked High"                      "$n_high"
printf '%-42s %10s\n' "  ranked Mid"                       "$n_mid"
printf '%-42s %10s\n' "  ranked Low (dropped)"             "$n_low"
printf '%-42s %10s\n' "R-scape covariation tested"         "$n_rscape"
printf '%-42s %10s\n' "failed"                             "$n_failed"

echo
[ "$n_failed" -eq 0 ] || warn "$n_failed cluster(s) failed — see $FAILED, then re-run to retry just those"
if [ "$(( n_high + n_mid ))" -eq 0 ]; then
  warn "no cluster ranked High or Mid — nothing reliable to carry into step 6"
else
  log "step 6 input: $HIGH ($n_high) + $MID ($n_mid) — motif_cleaned.sto per cluster"
fi
log "evaluation table: $EVALALL"
log "outputs in $OUTDIR"