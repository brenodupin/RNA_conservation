#!/usr/bin/env bash
#
# common.sh — helpers shared by every step script.
# Source this AFTER info.sh.

# --- messaging -------------------------------------------------------------

log()  { printf '[%s] %s\n'          "$(date +'%H:%M:%S')" "$*" >&2; }
warn() { printf '[%s] WARN  %s\n'    "$(date +'%H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR %s\n'    "$(date +'%H:%M:%S')" "$*" >&2; exit 1; }

# --- preflight checks ------------------------------------------------------

# require_cmd clustalo mmseqs ...
require_cmd() {
  local c
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 \
      || die "required command not on PATH: $c"
  done
}

require_dir() {
  [ -d "$1" ] || die "directory not found: $1"
}

require_file() {
  [ -s "$1" ] || die "file missing or empty: $1"
}

# --- logging ---------------------------------------------------------------

# start_log 00_oneline
# Tees everything from here on into $LOGDIR/<RUN_ID>_<step>.log as well as the
# terminal, and points $LOGDIR/latest_<step>.log at it.
#
# An EXIT trap restores the original file descriptors, which makes tee see EOF
# and flush. Without that, a script dying mid-run can lose its last few lines --
# exactly the lines you need when something failed.
start_log() {
  LOGSTEP=$1
  mkdir -p "$LOGDIR"
  LOGFILE="$LOGDIR/${RUN_ID}_${LOGSTEP}.log"

  exec 3>&1 4>&2                        # save real stdout/stderr
  exec > >(tee "$LOGFILE") 2>&1

  # Stable name so you never have to look up a timestamp:
  #   tail -f data/logs/latest_00_oneline.log
  ln -sfn "$(basename "$LOGFILE")" "$LOGDIR/latest_${LOGSTEP}.log"

  trap _close_log EXIT
  log "===== $LOGSTEP started (run $RUN_ID) ====="
}

_close_log() {
  local rc=$?
  if [ "$rc" -eq 0 ]; then
    log "===== $LOGSTEP finished OK ====="
  else
    log "===== $LOGSTEP FAILED (exit $rc) ====="
  fi
  log "log: $LOGFILE"
  exec 1>&3 2>&4 3>&- 4>&-              # restore; tee flushes and exits
}

# --- FASTA validation ------------------------------------------------------

# ensure_trailing_newline <file>
# The upstream fasta_oneLiner.sh leaves the final sequence line with no
# terminating newline. That is invisible until step 1 does
# `cat *_windows.fa > all_windows.fasta`, at which point the last sequence of
# one file is silently glued onto the first header of the next. Normalise here.
ensure_trailing_newline() {
  local f=$1
  [ -s "$f" ] || return 0
  if [ -n "$(tail -c1 "$f")" ]; then
    printf '\n' >> "$f"
  fi
}

# assert_oneline_fasta <file>
# createWindows.pl and removeNs_polyN_windows.pl both read the file as strict
# header/sequence LINE PAIRS (for i=0; i<$#lines; i+=2). Anything else silently
# produces garbage rather than erroring, so verify the shape explicitly.
assert_oneline_fasta() {
  local f=$1 lines heads
  require_file "$f"
  # awk counts a final line lacking a newline; wc -l does not.
  lines=$(awk 'END{print NR}' "$f")
  heads=$(grep -c '^>' "$f" || true)

  [ "$heads" -gt 0 ] || die "$f contains no FASTA headers"
  (( lines % 2 == 0 )) \
    || die "$f has an odd line count ($lines); expected strict header/sequence pairs"
  (( lines == heads * 2 )) \
    || die "$f has $heads headers across $lines lines; expected exactly one sequence line per header"
}

# seq_length <file>  -> total residues, gaps and newlines excluded
seq_length() {
  grep -v '^>' "$1" | tr -d '\n' | wc -c
}