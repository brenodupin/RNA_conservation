#!/usr/bin/env bash
#
# 01_windows.sh — slice each genome into overlapping windows, then filter.
#
#   in : $RNAC_DATA/$RNAC_STEP_00/*_oneLine.fasta
#   out: $RNAC_DATA/$RNAC_STEP_01/<sample>_windows.fa            per-sample windows
#        $RNAC_DATA/$RNAC_STEP_01/all_windows.fasta              combined
#        $RNAC_DATA/$RNAC_STEP_01/all_windows_processed_noNs.fasta
#        $RNAC_DATA/$RNAC_STEP_01/all_windows_processed_noNs_polyN.fasta
#        $RNAC_DATA/$RNAC_STEP_01/all_windows_processed_noNs_polyN_uniq.fasta   <-- step 2 input
#        $RNAC_DATA/$RNAC_STEP_01/all_windows_nuclComposition.tsv
#        $RNAC_DATA/$RNAC_STEP_01/<sample>_nuclInfo.csv           (if RNAC_RUN_COUNTNS=1)
#
#   ./01_windows.sh                        normal run
#   RNAC_FORCE=1 ./01_windows.sh                rebuild everything
#   RNAC_WINDOW=500 RNAC_OVERLAP=100 ./01_windows.sh different window geometry
#   RNAC_SORT_WINDOWS=0 ./01_windows.sh         random output ordering

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$RNAC_STEP_01"
INDIR="$RNAC_DATA/$RNAC_STEP_00"
OUTDIR="$RNAC_DATA/$STEP"
SRC="$RNAC_REPO/step1_createWindows"

start_log "$STEP"
log "window=$RNAC_WINDOW overlap=$RNAC_OVERLAP stride=$(( RNAC_WINDOW - RNAC_OVERLAP ))  force=$RNAC_FORCE"

require_cmd perl
require_dir "$INDIR"
require_file "$SRC/createWindows.pl"
require_file "$SRC/removeNs_polyN_windows.pl"
require_file "$SRC/removeDuplicates.pl"
mkdir -p "$OUTDIR"

# --- validate geometry -----------------------------------------------------
#
# createWindows.pl advances with  repos += RNAC_WINDOW - RNAC_OVERLAP  while testing
# (repos + RNAC_WINDOW) < size. If RNAC_OVERLAP >= RNAC_WINDOW the stride is zero or negative,
# repos never moves, and the loop writes the same window forever until the disk
# fills. It has no guard of its own, so check here.
#
# Its usage() also calls say() without `use feature`, so any bad-argument path
# dies with a confusing perl error rather than printing usage. All the more
# reason to validate before invoking it.

[[ "$RNAC_WINDOW"  =~ ^[0-9]+$ ]] || die "RNAC_WINDOW must be a positive integer (got '$RNAC_WINDOW')"
[[ "$RNAC_OVERLAP" =~ ^[0-9]+$ ]] || die "RNAC_OVERLAP must be a non-negative integer (got '$RNAC_OVERLAP')"
[ "$RNAC_WINDOW" -gt 0 ] || die "RNAC_WINDOW must be > 0"
[ "$RNAC_OVERLAP" -lt "$RNAC_WINDOW" ] \
  || die "RNAC_OVERLAP ($RNAC_OVERLAP) must be < RNAC_WINDOW ($RNAC_WINDOW): stride would be <= 0 and createWindows.pl would loop forever"

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

  if [ -s "$out" ] && [ "$RNAC_FORCE" != "1" ]; then
    log "skip $b (windows exist; RNAC_FORCE=1 to rebuild)"
    continue
  fi

  log "windowing $b (w=$RNAC_WINDOW p=$RNAC_OVERLAP)"
  rm -f "$out"
  perl "$SRC/createWindows.pl" -f "$f" -w "$RNAC_WINDOW" -p "$RNAC_OVERLAP" -O "$out"

  [ -s "$out" ] || die "createWindows.pl produced no output for $b (is the genome shorter than RNAC_WINDOW=$RNAC_WINDOW?)"
done

# --- combine ---------------------------------------------------------------
#
# Concatenated in the sorted sample order above rather than by re-globbing, so
# the combined file is reproducible across runs and machines.

COMBINED="$OUTDIR/${RNAC_WINDOWS_PREFIX}.fasta"

if [ -s "$COMBINED" ] && [ "$RNAC_FORCE" != "1" ]; then
  log "skip combine (${RNAC_WINDOWS_PREFIX}.fasta exists)"
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

FILT_PREFIX="$OUTDIR/${RNAC_WINDOWS_PREFIX}_processed"
NONS="${FILT_PREFIX}_noNs.fasta"
POLYN="${FILT_PREFIX}_noNs_polyN.fasta"
UNIQ="${FILT_PREFIX}_noNs_polyN_uniq.fasta"
COMPOSITION="$OUTDIR/${RNAC_WINDOWS_PREFIX}_nuclComposition.tsv"

if [ -s "$POLYN" ] && [ "$RNAC_FORCE" != "1" ]; then
  log "skip N/poly-N filter (output exists)"
else
  log "removing windows containing N and poly-A/T/G/C windows"
  perl "$SRC/removeNs_polyN_windows.pl" "$COMBINED" "$FILT_PREFIX" > "$COMPOSITION"
  require_file "$NONS"
  require_file "$POLYN"
fi

if [ -s "$UNIQ" ] && [ "$RNAC_FORCE" != "1" ]; then
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
# order. Set RNAC_SORT_WINDOWS=1 to sort by header and make runs reproducible.

if [ "$RNAC_SORT_WINDOWS" = "1" ]; then
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

if [ "$RNAC_RUN_COUNTNS" = "1" ]; then
  for f in "${inputs[@]}"; do
    b=$(basename "$f" _oneLine.fasta)
    csv="$OUTDIR/${b}_nuclInfo.csv"
    if [ -s "$csv" ] && [ "$RNAC_FORCE" != "1" ]; then continue; fi
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