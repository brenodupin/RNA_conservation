#!/usr/bin/env bash
#
# worker_step_05.sh — process ONE cluster for step 5.
#
#   worker_step_05.sh trim   <clusterName>   trimAlignment.pl (+ esl-reformat), host Perl / container esl
#   worker_step_05.sh eval   <clusterName>   RNA-SCoRE.pl, host Perl
#   worker_step_05.sh rscape <clusterName>   R-scape, in its own container
#
# One cluster per invocation, matching upstream: submit_trimAlignment_job.sh,
# submit_evalStruct_v3_job.sh and rscape_eval_twoTest_gc15.sh each looped one
# cluster at a time. 05_eval.sh runs many of these concurrently via xargs -P and
# gates the phases: only clusters RNA-SCoRE ranks High/Mid reach `rscape`.
#
# A separate sourced script rather than a function exported into xargs, for the
# reasons in worker_step_03.sh: `export -f` drops `set -euo pipefail` in the
# child and turns every unexported variable into a silent empty string.
#
# EVERY tool runs with the cluster directory as CWD. This is not cosmetic:
# RNA-SCoRE.pl writes side files (<cluster>_pair_matrix.txt,
# <cluster>_detected_hairpins.txt) to the CWD under a bare name, so two clusters
# sharing a CWD would clobber each other's side files. cd-ing per worker keeps
# each cluster's output in its own directory.
#
# Exits non-zero on failure; 05_eval.sh derives the per-phase failure lists by
# checking which clusters lack the expected marker afterwards, so nothing is
# written to a shared file from many processes at once.
#
# Debug a single cluster by hand:
#   ./worker_step_05.sh trim NC_008590.1_220_38500_38750
#   ./worker_step_05.sh eval NC_008590.1_220_38500_38750

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

mode="${1:-}"
name="${2:-}"
[ -n "$mode" ] && [ -n "$name" ] \
  || die "usage: worker_step_05.sh <trim|eval|rscape> <clusterName>"

# Perl programs vendored in the step directory (trimAlignment.pl ships with the
# repo; RNA-SCoRE.pl is added from RodrigoReisLab/RNA-SCoRE).
STEPSRC="$RNAC_REPO/step5_evaluationOfRNAstructures"
TRIM_PL="$STEPSRC/trimAlignment.pl"
SCORE_PL="$STEPSRC/RNA-SCoRE.pl"

LOCARNA_DIR="$RNAC_DATA/$RNAC_STEP_04"
OUTDIR="$RNAC_DATA/$RNAC_STEP_05"
dir="$OUTDIR/$name"

# --- helper: the trimmed alignment must be well-formed before anything reads it
#
# trimAlignment.pl and RNA-SCoRE.pl are strict line-structure parsers with a
# hardcoded header size ($n=3 and $n=2 respectively). If mlocarna's result.stk
# does not have exactly that header shape they emit a malformed .sto rather than
# erroring, which would then silently poison eval and R-scape. Catch it here.
assert_motif_sto() {
  local f=$1
  [ -s "$f" ] || { warn "empty/absent motif sto: $name ($f)"; return 1; }
  grep -q '#=GC SS_cons' "$f" || { warn "no SS_cons in motif sto: $name"; return 1; }
  # at least one sequence line (not blank, not a # comment, not the // terminator)
  grep -qE '^[^#/[:space:]]' "$f" || { warn "no sequences in motif sto: $name"; return 1; }
  return 0
}

case "$mode" in

  # --- trim ----------------------------------------------------------------
  #
  # in : $STEP_04/<name>/<name>_result.stk
  # out: <name>_motif.sto  <name>_motif.fasta  <name>_motif.afa  [<name>_motif.aln]
  #
  # marker for the phase: <name>_motif.sto
  trim)
    src="$LOCARNA_DIR/$name/${name}_result.stk"
    [ -s "$src" ] || { warn "no locarna result.stk: $name"; exit 1; }

    mkdir -p "$dir"
    cp -f "$src" "$dir/result.stk"     # upstream trims a file literally named result.stk
    cd "$dir"

    sto="${name}_motif.sto"
    fa="${name}_motif.fasta"
    afa="${name}_motif.afa"

    # Clear any partial output from an interrupted run so a leftover file can
    # never be mistaken for finished work.
    rm -f "$sto" "$fa" "$afa" "${name}_motif.aln"

    if ! err=$(perl "$TRIM_PL" result.stk "$sto" "$fa" "$afa" 2>&1); then
      warn "trimAlignment.pl failed: $name${err:+ — ${err//$'\n'/ }}"
      exit 1
    fi

    assert_motif_sto "$sto" || { rm -f "$sto"; exit 1; }

    # Optional clustal view via esl-reformat, in the container (it lives at
    # /rscape_v2.0.0.q/bin inside the image, not on the host).
    if [ "$RNAC_TRIM_MAKE_ALN" = "1" ]; then
      if ! docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
             -v "$dir":/work -w /work "$RNAC_IMAGE" \
             esl-reformat --informat stockholm -o "${name}_motif.aln" clustal "$sto" \
             >/dev/null 2>&1
      then
        warn "esl-reformat failed (non-fatal): $name"   # .aln is optional; keep going
      fi
    fi
    ;;

  # --- eval ----------------------------------------------------------------
  #
  # in : <name>_motif.sto
  # out: <name>_motif_cleaned.sto    (only if any sequence passed)
  #      <name>_motif_evaluated.tsv  (only if the cluster was evaluated)
  #      <name>_rank.tsv             (ALWAYS — the STDOUT rank line)
  #
  # marker for the phase: <name>_rank.tsv
  #
  # RNA-SCoRE prints one tab-separated line to STDOUT ending in High/Mid/Low.
  # That line is the ONLY output produced for every cluster: cleaned.sto is
  # written only when a sequence passes, and evaluated.tsv only when the cluster
  # reaches evaluation at all, so a Low cluster legitimately produces neither.
  # Capturing STDOUT per-cluster is therefore both the completion marker and the
  # gate signal for the rscape phase, and it avoids many processes appending to
  # one shared allClusters_evaluation.tsv (05_eval.sh collates the ranks at the
  # end instead).
  eval)
    cd "$dir" 2>/dev/null || { warn "no trim output dir: $name"; exit 1; }
    sto="${name}_motif.sto"
    assert_motif_sto "$sto" || { warn "trim output missing/bad, run trim first: $name"; exit 1; }

    cleaned="${name}_motif_cleaned.sto"
    evaluated="${name}_motif_evaluated.tsv"
    rank="${name}_rank.tsv"
    rm -f "$cleaned" "$evaluated" "$rank" "$rank.tmp"

    # STDERR kept separate so a warning never lands in the rank line. All four
    # thresholds passed explicitly (see info.sh: an omitted one prints a notice
    # to STDOUT and would corrupt the capture).
    if ! perl "$SCORE_PL" \
           -e _motif.sto -d "$RNAC_SCORE_DUPL" \
           --mt "$RNAC_SCORE_MT" -t "$RNAC_SCORE_BP" --gc "$RNAC_SCORE_GC" \
           "$sto" "$cleaned" "$evaluated" > "$rank.tmp" 2>"${name}_rnascore.err"
    then
      warn "RNA-SCoRE.pl failed: $name — see $dir/${name}_rnascore.err"
      rm -f "$rank.tmp"
      exit 1
    fi

    # A successful run must emit exactly the rank line. No line = something is
    # wrong even though perl exited 0 (e.g. an empty alignment), so fail rather
    # than record an empty rank.
    [ -s "$rank.tmp" ] || { warn "RNA-SCoRE produced no rank line: $name"; rm -f "$rank.tmp"; exit 1; }
    mv -f "$rank.tmp" "$rank"
    ;;

  # --- rscape --------------------------------------------------------------
  #
  # in : <name>_motif_cleaned.sto        (only High/Mid clusters have this)
  # out: <name>_<label>_rscape_out/      R-scape's output directory
  #      <name>_rscape.log               R-scape STDOUT/STDERR
  #      <name>_rscape.done              written last, on success — phase marker
  #
  # Why a sentinel rather than a real R-scape output file: R-scape's exact
  # output filenames vary with version and with whether covariation was found,
  # so keying completion on e.g. ".cacofold.sto" risks marking a genuinely
  # finished run as failed. The worker writes <name>_rscape.done only after
  # R-scape exits 0, which is unambiguous and version-independent.
  rscape)
    cd "$dir" 2>/dev/null || { warn "no cluster dir: $name"; exit 1; }
    cleaned="${name}_motif_cleaned.sto"
    [ -s "$cleaned" ] || { warn "no cleaned alignment (not High/Mid?): $name"; exit 1; }

    outdir="${name}_${RNAC_RSCAPE_LABEL}_rscape_out"
    outname="${name}_${RNAC_RSCAPE_LABEL}_rscape"
    done_marker="${name}_rscape.done"

    # R-scape refuses to reuse a populated --outdir; clear any partial one.
    rm -rf "$outdir" "$done_marker"
    mkdir -p "$outdir"                  # README: outdir must exist beforehand

    if ! docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
           -v "$dir":/work -w /work "$RNAC_IMAGE" \
           R-scape -s --cacofold --outdir "$outdir" --outname "$outname" \
             --voutput --seed "$RNAC_RSCAPE_SEED" "$cleaned" \
           > "${name}_rscape.log" 2>&1
    then
      warn "R-scape failed: $name — see $dir/${name}_rscape.log"
      exit 1
    fi

    printf 'ok\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" > "$done_marker"
    ;;

  *)
    die "unknown mode '$mode' (expected: trim | eval | rscape)"
    ;;
esac