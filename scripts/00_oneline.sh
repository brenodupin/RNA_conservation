#!/usr/bin/env bash
#
# 00_oneline.sh — convert every input FASTA to strict one-line format.
#
#   in :  $RNAC_INPUT/*.{fasta,fa,fna}
#   out:  $RNAC_DATA/$RNAC_STEP_00/<sample>_oneLine.fasta
#
# Every downstream perl script in this pipeline reads FASTA as header/sequence
# LINE PAIRS. Multi-line input does not error, it just produces silent garbage,
# so this step is mandatory and validates its own output.
#
#   ./00_oneline.sh                      normal run
#   RNAC_FORCE=1 ./00_oneline.sh              rebuild existing outputs
#   RNAC_ONELINE_METHOD=repo ./00_oneline.sh  use upstream fasta_oneLiner.sh instead

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$RNAC_STEP_00"
OUTDIR="$RNAC_DATA/$STEP"

start_log "$STEP"
log "RNAC_REPO=$RNAC_REPO"
log "RNAC_INPUT=$RNAC_INPUT"
log "method=$RNAC_ONELINE_METHOD  force=$RNAC_FORCE"

require_dir "$RNAC_INPUT"
mkdir -p "$OUTDIR"

# --- collect inputs --------------------------------------------------------

shopt -s nullglob
inputs=( "$RNAC_INPUT"/*.fasta "$RNAC_INPUT"/*.fa "$RNAC_INPUT"/*.fna )
shopt -u nullglob

[ ${#inputs[@]} -gt 0 ] || die "no .fasta/.fa/.fna files found in $RNAC_INPUT"
log "found ${#inputs[@]} input file(s)"

# Upstream derives the sample name with  basename | cut -d. -f1  so
# AP019006.1.fasta becomes AP019006. We replicate that exactly to stay
# consistent with the repo, which means two inputs sharing a prefix before the
# first dot would overwrite each other. Catch that before it happens.
declare -A seen=()
for f in "${inputs[@]}"; do
  b=$(basename "$f" | cut -d. -f1)
  if [ -n "${seen[$b]:-}" ]; then
    die "name collision: '$(basename "$f")' and '${seen[$b]}' both reduce to '$b'"
  fi
  seen["$b"]=$(basename "$f")
done

# --- convert ---------------------------------------------------------------

converted=0 skipped=0

for f in "${inputs[@]}"; do
  b=$(basename "$f" | cut -d. -f1)
  out="$OUTDIR/${b}_oneLine.fasta"

  if [ -s "$out" ] && [ "$RNAC_FORCE" != "1" ]; then
    log "skip $b (output exists; RNAC_FORCE=1 to rebuild)"
    skipped=$((skipped + 1))
  else
    log "converting $(basename "$f") -> $(basename "$out")"
    rm -f "$out"

    case "$RNAC_ONELINE_METHOD" in
      awk)
        # Strips all whitespace inside sequences and emits a final record even
        # when the input lacks a trailing newline (upstream's `while read` loop
        # drops that last line).
        awk '
          /^>/ { if (seq != "") print seq; print; seq = ""; next }
                { gsub(/[[:space:]]/, ""); seq = seq $0 }
          END   { if (seq != "") print seq }
        ' "$f" > "$out"
        ;;
      repo)
        # Upstream writes to the current directory, so run it from OUTDIR.
        ( cd "$OUTDIR" \
          && bash "$RNAC_REPO/step0_convertInputSequencesToReqFASTAformat/fasta_oneLiner.sh" "$f" )
        ;;
      *)
        die "unknown RNAC_ONELINE_METHOD '$RNAC_ONELINE_METHOD' (expected: awk | repo)"
        ;;
    esac

    # Upstream leaves no trailing newline on the final sequence line; step 1
    # concatenates these files, so normalise before anything downstream sees it.
    ensure_trailing_newline "$out"
    converted=$((converted + 1))
  fi

  assert_oneline_fasta "$out"

  # Residue-count guard. The upstream script drives a `while read` loop, which
  # drops the final line of any file that does not end in a newline -- silently,
  # and only for the last record. Compare in against out so no method can lose
  # sequence without failing loudly.
  in_bp=$(grep -v '^>' "$f" | tr -d '[:space:]' | wc -c)
  out_bp=$(seq_length "$out")
  if [ "$in_bp" -ne "$out_bp" ]; then
    die "$b lost sequence during conversion: input $in_bp bp, output $out_bp bp"$'\n'"       (RNAC_ONELINE_METHOD=repo truncates files with no trailing newline; use awk)"
  fi
done

# --- report ----------------------------------------------------------------

echo
printf '%-14s %8s %12s  %s\n' SAMPLE RECORDS LENGTH_BP OUTPUT
printf '%-14s %8s %12s  %s\n' -------------- -------- ------------ ------

for f in "${inputs[@]}"; do
  b=$(basename "$f" | cut -d. -f1)
  out="$OUTDIR/${b}_oneLine.fasta"

  heads=$(grep -c '^>' "$out")
  len=$(seq_length "$out")

  printf '%-14s %8s %12s  %s\n' "$b" "$heads" "$len" "$(basename "$out")"

  # Non-ACGT characters: Ns are fine here but every window containing one is
  # discarded in step 1, and a U means the input is RNA rather than DNA.
  odd=$(grep -v '^>' "$out" | tr -d '\n' | tr -d 'ACGTacgt' | wc -c)
  if [ "$odd" -gt 0 ]; then
    warn "$b has $odd non-ACGT character(s); windows containing N are dropped in step 1"
  fi
  if grep -v '^>' "$out" | grep -qi 'u'; then
    warn "$b appears to contain U — this pipeline expects DNA, not RNA"
  fi
done

echo
log "converted=$converted skipped=$skipped"

log "outputs in $OUTDIR"