#!/usr/bin/env bash
#
# worker_step_04.sh — run mlocarna on ONE cluster, in its own container.
#
#   worker_step_04.sh <clusterName>
#
# One cluster per invocation, matching upstream: locarnap_array_job.sh was a
# SLURM array task per cluster. 04_locarna.sh runs many of these concurrently
# via xargs -P, exactly as 03_screen.sh does with worker_step_03.sh.
#
# A separate script rather than a function exported into xargs, for the reasons
# spelled out at the top of worker_step_03.sh: `export -f` loses `set -euo
# pipefail` and every unexported variable silently becomes an empty string.
#
# Exits non-zero on failure; 04_locarna.sh derives the failure list by checking
# which clusters lack output afterwards, so nothing is written to a shared file.
#
# Debug a single cluster by hand:
#   ./worker_step_04.sh NC_008590.1_220_38500_38750

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

name="${1:-}"
[ -n "$name" ] || die "usage: worker_step_04.sh <clusterName>"

SPLITDIR="$RNAC_DATA/$RNAC_STEP_02/splits"
OUTDIR="$RNAC_DATA/$RNAC_STEP_04"

src="$SPLITDIR/${name}_cluster.fasta"
dir="$OUTDIR/$name"
tgt="$dir/${name}_locarnap"          # mlocarna's --tgtdir
result="$dir/${name}_result.stk"     # completion marker, see below
timing="$dir/${name}_timing.tsv"
mlog="$dir/${name}_mlocarna.log"

[ -s "$src" ] || { warn "no cluster FASTA: $name"; exit 1; }

mkdir -p "$dir"
cp -f "$src" "$dir/"

# mlocarna will happily start inside a half-finished --tgtdir from an interrupted
# run and produce a mixture of old and new intermediates. There is no resume
# worth trusting here, so start clean every time we decide to run at all.
# (04_locarna.sh has already decided this cluster needs work.)
rm -rf "$tgt" "$result" "$timing"

nseqs=$(grep -c '^>' "$dir/${name}_cluster.fasta" || true)

start_ts="$(date '+%Y-%m-%d %H:%M:%S')"
start_epoch=$(date +%s)

# Upstream binds ../${clusterName} to /input and runs with that directory as the
# CWD, so --tgtdir lands beside the FASTA. The docker equivalent is mounting the
# cluster directory at /work and making it the working directory -- one mount
# instead of two, and it also catches anything mlocarna drops in the CWD.
#
# --cpus matches --threads so that RNAC_JOBS workers use RNAC_JOBS x
# RNAC_LOCARNA_THREADS cores in total and no more.
cpus=( --cpus="$RNAC_LOCARNA_THREADS" )
mem=()
[ -n "$RNAC_LOCARNA_MEM" ] && mem=( --memory="$RNAC_LOCARNA_MEM" )

# `set -e` would abort the script the moment docker returned non-zero, before
# the exit code could be recorded, so capture it explicitly instead.
rc=0
docker run --rm "${cpus[@]}" "${mem[@]}" -u "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$dir":/work -w /work "$RNAC_IMAGE" \
  mlocarna --probabilistic --tgtdir "${name}_locarnap" --moreverbose \
    --stockholm --local-progressive \
    --threads="$RNAC_LOCARNA_THREADS" \
    --rnafold-temperature="$RNAC_FOLD_TEMP" \
    "${name}_cluster.fasta" \
  > "$mlog" 2>&1 || rc=$?

end_ts="$(date '+%Y-%m-%d %H:%M:%S')"
duration=$(( $(date +%s) - start_epoch ))

# One timing file per cluster rather than RNAC_JOBS processes appending to a
# shared log. 04_locarna.sh collates them at the end. Same reasoning as the
# failure list in step 3: no locking, and nothing to interleave.
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$name" "$nseqs" "$start_ts" "$end_ts" "$duration" "$rc" > "$timing"

if [ "$rc" -ne 0 ]; then
  # 137 = SIGKILL, which for a memory-capped container is almost always the
  # kernel OOM killer rather than anything mlocarna did wrong. Worth naming,
  # because the mlocarna log just stops mid-sentence and looks like a crash.
  if [ "$rc" -eq 137 ] && [ -n "$RNAC_LOCARNA_MEM" ]; then
    warn "mlocarna killed (OOM?): $name — ${nseqs} seqs hit RNAC_LOCARNA_MEM=$RNAC_LOCARNA_MEM"
  else
    warn "mlocarna failed: $name (exit $rc, ${nseqs} seqs) — see $mlog"
  fi
  exit 1
fi

# --- locate result.stk -----------------------------------------------------
#
# mlocarna --stockholm writes the alignment to <tgtdir>/results/result.stk.
# Copying it to <dir>/<name>_result.stk gives step 4 the same
# "<cluster dir>/<cluster name>_<thing>" shape as every other step, so the
# driver can count finished clusters at a fixed depth and step 5 has one
# predictable path per cluster instead of a nested tgtdir layout.
#
# The fallback search covers mlocarna versions that place it elsewhere under
# the tgtdir -- cheap, and better than failing a run that actually succeeded.
stk="$tgt/results/result.stk"
if [ ! -s "$stk" ]; then
  stk=$(find "$tgt" -type f -name 'result.stk' -print -quit 2>/dev/null || true)
fi

if [ -z "$stk" ] || [ ! -s "$stk" ]; then
  warn "mlocarna exited 0 but produced no result.stk: $name — see $mlog"
  exit 1
fi

cp -f "$stk" "$result"