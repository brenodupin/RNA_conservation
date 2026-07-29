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
# Tees everything from here on into $RNAC_LOGDIR/<RNAC_RUN_ID>_<step>.log as well as the
# terminal, and points $RNAC_LOGDIR/latest_<step>.log at it.
#
# An EXIT trap restores the original file descriptors, which makes tee see EOF
# and flush. Without that, a script dying mid-run can lose its last few lines --
# exactly the lines you need when something failed.
start_log() {
  LOGSTEP=$1
  mkdir -p "$RNAC_LOGDIR"
  LOGFILE="$RNAC_LOGDIR/${RNAC_RUN_ID}_${LOGSTEP}.log"

  exec 3>&1 4>&2                        # save real stdout/stderr
  exec > >(tee "$LOGFILE") 2>&1

  # Stable name so you never have to look up a timestamp:
  #   tail -f data/logs/latest_00_oneline.log
  ln -sfn "$(basename "$LOGFILE")" "$RNAC_LOGDIR/latest_${LOGSTEP}.log"

  trap _close_log EXIT
  log "===== $LOGSTEP started (run $RNAC_RUN_ID) ====="

  # Surface any RNAC_* overrides in effect. A per-command prefix and a stale
  # `export` are indistinguishable from inside the script, so the only defence
  # against the latter is to print whatever is actually set.
  if [ "${#RNAC_OVERRIDES[@]}" -gt 0 ]; then
    local v
    for v in "${RNAC_OVERRIDES[@]}"; do
      log "override: $v=${!v}"
    done
  fi
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

# --- conda ------------------------------------------------------------------
#
# conda_activate <env> [required_cmd ...]
#
# `conda activate` is a shell function, not an executable, and it is undefined
# in a non-interactive script until conda's hook is evaluated -- without it you
# get "CommandNotFoundError: Your shell has not been properly configured".
# Conda's own scripts also reference unset variables, so `set -u` is relaxed
# around them.
#
# The resulting PATH change persists for the remainder of the script and is
# inherited by every command it runs. It does NOT persist if this is called
# from inside ( ) or $( ), so call it at the top level of a step.
conda_activate() {
  local env=$1; shift

  command -v conda >/dev/null 2>&1 \
    || die "conda not on PATH — needed to activate '$env'"

  set +u
  if [ -z "${_CONDA_HOOKED:-}" ]; then
    eval "$(conda shell.bash hook)"
    _CONDA_HOOKED=1
  fi
  set -u

  conda env list | awk '{print $1}' | grep -qx -- "$env" \
    || die "conda env '$env' does not exist — create it with: conda env create -f $RNAC_REPO/step2_clustering/${env}*.yml"

  set +u
  if ! conda activate "$env"; then set -u; die "failed to activate conda env '$env'"; fi
  set -u
  log "conda env active: $env"

  local c
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 \
      || die "'$c' not found in conda env '$env'"
  done
}

conda_deactivate() {
  set +u
  conda deactivate 2>/dev/null || true
  set -u
}

# --- docker -----------------------------------------------------------------
#
# Checked once per step, before any long loop, so a missing image or a
# permission problem fails immediately rather than after the first cluster.
docker_preflight() {
  command -v docker >/dev/null 2>&1 \
    || die "docker not on PATH"
  docker info >/dev/null 2>&1 \
    || die "cannot reach the docker daemon. If this is a permission error: sudo usermod -aG docker \$USER, then log out and back in"
  docker image inspect "$RNAC_IMAGE" >/dev/null 2>&1 \
    || die "image '$RNAC_IMAGE' not present locally — run: docker pull $RNAC_IMAGE"
  log "docker image: $RNAC_IMAGE"
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