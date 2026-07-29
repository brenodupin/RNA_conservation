#!/usr/bin/env bash
#
# 01_windows.sh — slice each genome into overlapping windows, then filter.
#
#   in : $DATA/$STEP_00/*_oneLine.fasta
#   out: $DATA/$STEP_01/<sample>_windows.fa            per-sample windows
#        $DATA/$STEP_01/all_windows.fasta              combined
#        $DATA/$STEP_01/all_windows_processed_noNs.fasta
#        $DATA/$STEP_01/all_windows_processed_noNs_polyN.fasta
#        $DATA/$STEP_01/all_windows_processed_noNs_polyN_uniq.fasta   <-- step 2 input
#        $DATA/$STEP_01/all_windows_nuclComposition.tsv
#        $DATA/$STEP_01/<sample>_nuclInfo.csv           (if RUN_COUNTNS=1)
#
#   ./01_windows.sh                        normal run
#   FORCE=1 ./01_windows.sh                rebuild everything
#   WINDOW=500 OVERLAP=100 ./01_windows.sh different window geometry
#   SORT_WINDOWS=0 ./01_windows.sh         random output ordering

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$STEP_01"
INDIR="$DATA/$STEP_00"
OUTDIR="$DATA/$STEP"
SRC="$REPO/step1_createWindows"

start_log "$STEP"
log "window=$WINDOW overlap=$OVERLAP stride=$(( WINDOW - OVERLAP ))  force=$FORCE"

require_cmd perl
require_dir "$INDIR"
require_file "$SRC/createWindows.pl"
require_file "$SRC/removeNs_polyN_windows.pl"
require_file "$SRC/removeDuplicates.pl"
mkdir -p "$OUTDIR"

# --- validate geometry -----------------------------------------------------
#
# createWindows.pl advances with  repos += WINDOW - OVERLAP  while testing
# (repos + WINDOW) < size. If OVERLAP >= WINDOW the stride is zero or negative,
# repos never moves, and the loop writes the same window forever until the disk
# fills. It has no guard of its own, so check here.
#
# Its usage() also calls say() without `use feature`, so any bad-argument path
# dies with a confusing perl error rather than printing usage. All the more
# reason to validate before invoking it.

[[ "$WINDOW"  =~ ^[0-9]+$ ]] || die "WINDOW must be a positive integer (got '$WINDOW')"
[[ "$OVERLAP" =~ ^[0-9]+$ ]] || die "OVERLAP must be a non-negative integer (got '$OVERLAP')"
[ "$WINDOW" -gt 0 ] || die "WINDOW must be > 0"
[ "$OVERLAP" -lt "$WINDOW" ] \
  || die "OVERLAP ($OVERLAP) must be < WINDOW ($WINDOW): stride would be <= 0 and createWindows.pl would loop forever"

# --- collect inputs --------------------------------------------------------

shopt -s nullglob
mapfile -t inputs < <(printf '%s\n' "$INDIR"/*_oneLine.fasta | LC_ALL=C sort)
shopt -u nullglob

[ ${#inputs[@]} -gt 0 ] \
  || die "no *_oneLine.fasta in $INDIR — run 00_oneline.sh first"
log "found ${#inputs[@]} one-line genome(s)"

# Step 0 guarantees strict header/sequence pairs, but re-check: createWindows.pl
# reads the file as line pairs and produces silent garbage on anything else,
# which is exactly the failure that is hardest to notice later.
for f in "${inputs[@]}"; do
  assert_oneline_fasta "$f"
done

# --- window each genome ----------------------------------------------------

for f in "${inputs[@]}"; do
  b=$(basename "$f" _oneLine.fasta)
  out="$OUTDIR/${b}_windows.fa"

  if [ -s "$out" ] && [ "$FORCE" != "1" ]; then
    log "skip $b (windows exist; FORCE=1 to rebuild)"
    continue
  fi

  log "windowing $b (w=$WINDOW p=$OVERLAP)"
  rm -f "$out"
  perl "$SRC/createWindows.pl" -f "$f" -w "$WINDOW" -p "$OVERLAP" -O "$out"

  [ -s "$out" ] || die "createWindows.pl produced no output for $b (is the genome shorter than WINDOW=$WINDOW?)"
done

# --- combine ---------------------------------------------------------------
#
# Concatenated in the sorted sample order above rather than by re-globbing, so
# the combined file is reproducible across runs and machines.

COMBINED="$OUTDIR/${WINDOWS_PREFIX}.fasta"

if [ -s "$COMBINED" ] && [ "$FORCE" != "1" ]; then
  log "skip combine (${WINDOWS_PREFIX}.fasta exists)"
else
  log "combining ${#inputs[@]} window file(s) -> $(basename "$COMBINED")"
  : > "$COMBINED"
  for f in "${inputs[@]}"; do
    b=$(basename "$f" _oneLine.fasta)
    cat "$OUTDIR/${b}_windows.fa" >> "$COMBINED"
  done
  ensure_trailing_newline "$COMBINED"
fi

assert_oneline_fasta "$COMBINED"

# --- filter ----------------------------------------------------------------
#
# removeNs_polyN_windows.pl takes  <input.fasta> <output PREFIX>  and writes
# <prefix>_noNs.fasta and <prefix>_noNs_polyN.fasta. Its stdout is a per-window
# nucleotide composition table, which must be redirected or it floods the log.
#
# removeDuplicates.pl takes a PREFIX, not a filename — it appends .fasta itself
# and writes <prefix>_uniq.fasta. Passing a filename makes it look for
# "....fasta.fasta", open nothing, and silently emit an empty result.

FILT_PREFIX="$OUTDIR/${WINDOWS_PREFIX}_processed"
NONS="${FILT_PREFIX}_noNs.fasta"
POLYN="${FILT_PREFIX}_noNs_polyN.fasta"
UNIQ="${FILT_PREFIX}_noNs_polyN_uniq.fasta"
COMPOSITION="$OUTDIR/${WINDOWS_PREFIX}_nuclComposition.tsv"

if [ -s "$POLYN" ] && [ "$FORCE" != "1" ]; then
  log "skip N/poly-N filter (output exists)"
else
  log "removing windows containing N and poly-A/T/G/C windows"
  perl "$SRC/removeNs_polyN_windows.pl" "$COMBINED" "$FILT_PREFIX" > "$COMPOSITION"
  require_file "$NONS"
  require_file "$POLYN"
fi

if [ -s "$UNIQ" ] && [ "$FORCE" != "1" ]; then
  log "skip deduplication (output exists)"
else
  log "removing duplicate window sequences"
  # NOTE: prefix, not filename.
  perl "$SRC/removeDuplicates.pl" "${FILT_PREFIX}_noNs_polyN"
  require_file "$UNIQ"
fi

# --- optional deterministic ordering ---------------------------------------
#
# Both perl filters iterate `keys %hash`, and perl randomises hash order per
# process, so reruns produce the same records in a different order. That is
# harmless for correctness but makes byte-level diffing useless, and mmseqs2 in
# step 2 can pick different cluster representatives from a different input
# order. Set SORT_WINDOWS=1 to sort by header and make runs reproducible.

if [ "$SORT_WINDOWS" = "1" ]; then
  log "sorting filtered outputs by header for reproducibility"
  for target in "$NONS" "$POLYN" "$UNIQ"; do
    [ -s "$target" ] || continue
    paste - - < "$target" | LC_ALL=C sort -k1,1 | tr '\t' '\n' > "${target}.sorted"
    mv "${target}.sorted" "$target"
  done
fi

for target in "$NONS" "$POLYN" "$UNIQ"; do
  assert_oneline_fasta "$target"
done

# --- optional per-genome N content -----------------------------------------
#
# countNs.pl emits: header <tab> length <tab> N_count <tab> percent_N
# It counts anything outside AUTGCautgc, so ambiguity codes (R, Y, W...) are
# included in the N count, not just literal N.

if [ "$RUN_COUNTNS" = "1" ]; then
  for f in "${inputs[@]}"; do
    b=$(basename "$f" _oneLine.fasta)
    csv="$OUTDIR/${b}_nuclInfo.csv"
    if [ -s "$csv" ] && [ "$FORCE" != "1" ]; then continue; fi
    perl "$SRC/countNs.pl" "$f" > "$csv"
  done
  log "per-genome N content written to <sample>_nuclInfo.csv"
fi

# --- report ----------------------------------------------------------------

echo
printf '%-16s %10s  %s\n' SAMPLE WINDOWS OUTPUT
printf '%-16s %10s  %s\n' ---------------- ---------- ------
for f in "${inputs[@]}"; do
  b=$(basename "$f" _oneLine.fasta)
  printf '%-16s %10s  %s\n' "$b" "$(grep -c '^>' "$OUTDIR/${b}_windows.fa")" "${b}_windows.fa"
done

n_all=$(grep -c '^>' "$COMBINED")
n_nons=$(grep -c '^>' "$NONS")
n_poly=$(grep -c '^>' "$POLYN")
n_uniq=$(grep -c '^>' "$UNIQ")

echo
printf '%-34s %10s %12s\n' STAGE WINDOWS REMOVED
printf '%-34s %10s %12s\n' ---------------------------------- ---------- ------------
printf '%-34s %10s %12s\n' "created"                "$n_all"  "-"
printf '%-34s %10s %12s\n' "after removing N-containing" "$n_nons" "$(( n_all  - n_nons ))"
printf '%-34s %10s %12s\n' "after removing poly-A/T/G/C" "$n_poly" "$(( n_nons - n_poly ))"
printf '%-34s %10s %12s\n' "after deduplication"    "$n_uniq" "$(( n_poly - n_uniq ))"

echo
[ "$n_uniq" -gt 0 ] || die "no windows survived filtering — nothing for step 2 to cluster"

if [ "$n_uniq" -lt 100 ]; then
  warn "only $n_uniq unique windows; clustering will produce very few multi-member clusters"
fi

log "step 2 input: $(basename "$UNIQ") ($n_uniq windows)"
log "outputs in $OUTDIR"