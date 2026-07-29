#!/usr/bin/env bash
#
# 02_clusters.sh — add reverse-complement windows, cluster in two mmseqs2
#                  passes, and split each cluster into its own FASTA.
#
#   in : $DATA/$STEP_01/<prefix>_processed_noNs_polyN_uniq.fasta
#   out: $DATA/$STEP_02/<prefix>_revComp.fasta
#        $DATA/$STEP_02/<prefix>_bothStrands.fasta        clustering input
#        $DATA/$STEP_02/mmseq2_covmode0_PID95_cov80_*     pass 1 (redundancy)
#        $DATA/$STEP_02/mmseq2_repSeqs_covmode0_PID50_cov80_*   pass 2 (clusters)
#        $DATA/$STEP_02/cluster_count.tsv                 members per cluster
#        $DATA/$STEP_02/splits/<rep>_cluster.fasta        one file per cluster
#
#   ./02_clusters.sh                     normal run
#   FORCE=1 ./02_clusters.sh             rebuild everything
#   RUN_REVCOMP=0 ./02_clusters.sh       strand-specific input, skip minus strand
#   SPLIT_METHOD=awk ./02_clusters.sh    fast splitter (identical output)
#   PID_PASS2=0.40 ./02_clusters.sh      looser cluster threshold
#
# Upstream marks this whole step optional: skip it if your input sequences are
# already a defined, related set (e.g. UTRs of one gene family) rather than
# genome-wide windows of unknown similarity.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=info.sh
source "$HERE/info.sh"
# shellcheck source=common.sh
source "$HERE/common.sh"

STEP="$STEP_02"
INDIR="$DATA/$STEP_01"
OUTDIR="$DATA/$STEP"
SPLITDIR="$OUTDIR/splits"
SRC="$REPO/step2_clustering"

start_log "$STEP"

UNIQ="$INDIR/${WINDOWS_PREFIX}_processed_noNs_polyN_uniq.fasta"
require_file "$UNIQ"
mkdir -p "$OUTDIR"

# --- validate thresholds ---------------------------------------------------

frac_ok() { awk -v v="$1" 'BEGIN{ exit !(v > 0 && v <= 1) }'; }
pct()     { awk -v v="$1" 'BEGIN{ printf "%d", v*100 + 0.5 }'; }

for v in PID_PASS1 PID_PASS2 COVERAGE; do
  frac_ok "${!v}" || die "$v must be a fraction in (0,1], got '${!v}'"
done
awk -v a="$PID_PASS2" -v b="$PID_PASS1" 'BEGIN{ exit !(a <= b) }' \
  || die "PID_PASS2 ($PID_PASS2) should not exceed PID_PASS1 ($PID_PASS1): pass 2 clusters the survivors of pass 1 at a looser threshold"

log "pass1 min-seq-id=$PID_PASS1  pass2 min-seq-id=$PID_PASS2  coverage=$COVERAGE  threads=$THREADS"

# --- reverse complement ----------------------------------------------------
#
# seqtk emits one line per sequence by default (-l 0), which the rest of the
# pipeline depends on. Passed explicitly rather than relied upon.
#
# The trailing 'r' appended to each header is what marks a window as minus
# strand; it survives all the way to the final motif coordinates, so do not
# change it without checking step 6's coordinate handling.

REVCOMP="$OUTDIR/${WINDOWS_PREFIX}_revComp.fasta"
BOTH="$OUTDIR/${WINDOWS_PREFIX}_bothStrands.fasta"

if [ "$RUN_REVCOMP" = "1" ]; then
  if [ -s "$BOTH" ] && [ "$FORCE" != "1" ]; then
    log "skip reverse complement (output exists)"
  else
    conda_activate "$CONDA_ENV_SEQTK" seqtk
    log "generating reverse-complement windows"
    seqtk seq -r -l 0 "$UNIQ" > "$REVCOMP"
    sed -i '/^>/s/$/r/' "$REVCOMP"
    ensure_trailing_newline "$REVCOMP"
    assert_oneline_fasta "$REVCOMP"

    cat "$UNIQ" "$REVCOMP" > "$BOTH"
    ensure_trailing_newline "$BOTH"
  fi

  assert_oneline_fasta "$BOTH"
  n_uniq=$(grep -c '^>' "$UNIQ")
  n_both=$(grep -c '^>' "$BOTH")
  [ "$n_both" -eq $(( n_uniq * 2 )) ] \
    || die "expected $(( n_uniq * 2 )) sequences after adding reverse complements, got $n_both"
  CLUSTER_INPUT="$BOTH"
else
  log "RUN_REVCOMP=0 — clustering the plus strand only"
  CLUSTER_INPUT="$UNIQ"
fi

# --- mmseqs2 ---------------------------------------------------------------
#
# Two passes, matching upstream:
#   1. --min-seq-id 0.95 collapses near-identical windows; its _rep_seq.fasta
#      keeps one representative of each redundant group.
#   2. --min-seq-id 0.50 groups those representatives into the clusters that
#      the rest of the pipeline treats as candidate RNA families.
# Only pass 1 gets --kmer-per-seq; upstream omits it from pass 2.

P1="$OUTDIR/mmseq2_covmode0_PID$(pct "$PID_PASS1")_cov$(pct "$COVERAGE")"
P2="$OUTDIR/mmseq2_repSeqs_covmode0_PID$(pct "$PID_PASS2")_cov$(pct "$COVERAGE")"
TMP1="$OUTDIR/tmp_covmode0_PID$(pct "$PID_PASS1")_cov$(pct "$COVERAGE")"
TMP2="$OUTDIR/tmp_repSeqs_covmode0_PID$(pct "$PID_PASS2")_cov$(pct "$COVERAGE")"

if [ -s "${P1}_rep_seq.fasta" ] && [ "$FORCE" != "1" ]; then
  log "skip mmseqs pass 1 (output exists)"
else
  conda_activate "$CONDA_ENV_MMSEQS" mmseqs
  log "mmseqs pass 1: collapsing windows >= ${PID_PASS1} identical"
  mmseqs easy-cluster "$CLUSTER_INPUT" "$P1" "$TMP1" \
    -c "$COVERAGE" --threads "$THREADS" --kmer-per-seq "$KMER_PER_SEQ" \
    --min-seq-id "$PID_PASS1" --cov-mode 0 --filter-hits 1
  require_file "${P1}_rep_seq.fasta"
fi

if [ -s "${P2}_cluster.tsv" ] && [ "$FORCE" != "1" ]; then
  log "skip mmseqs pass 2 (output exists)"
else
  conda_activate "$CONDA_ENV_MMSEQS" mmseqs
  log "mmseqs pass 2: clustering representatives at >= ${PID_PASS2} identity"
  mmseqs easy-cluster "${P1}_rep_seq.fasta" "$P2" "$TMP2" \
    -c "$COVERAGE" --threads "$THREADS" \
    --min-seq-id "$PID_PASS2" --cov-mode 0 --filter-hits 1
  require_file "${P2}_cluster.tsv"
  require_file "${P2}_all_seqs.fasta"
fi

# --- cluster_count.tsv -----------------------------------------------------
#
# Upstream: cut -f1 ..._cluster.tsv | sort | uniq -dc > cluster_count.tsv
#
# The README calls the result "tab-separated", but `uniq -c` emits right-aligned
# space padding, e.g. "     38 NC_008590.1_220_38500_38750". getClusterSequences.sh
# then reads it with IFS=$'\t', so every field lands in $i, $j stays empty, and
# it writes a directory full of empty _cluster.fasta files without erroring.
# The sed below is what makes the two halves agree.
#
# `uniq -d` drops clusters seen only once: a single-member cluster has nothing
# to align against, so it is deliberately excluded from everything downstream.

COUNTS="$OUTDIR/cluster_count.tsv"

if [ -s "$COUNTS" ] && [ "$FORCE" != "1" ]; then
  log "skip cluster_count.tsv (exists)"
else
  cut -f1 "${P2}_cluster.tsv" | sort | uniq -dc \
    | sed -e 's/^ *//' -e 's/ \+/\t/' > "$COUNTS"
fi

# Fail loudly if the tab conversion ever regresses.
if [ -s "$COUNTS" ] && ! head -1 "$COUNTS" | grep -q $'\t'; then
  die "cluster_count.tsv is not tab-separated — getClusterSequences.sh would produce empty files"
fi

n_multi=$(wc -l < "$COUNTS")
log "clusters with >= 2 members: $n_multi"

if [ "$n_multi" -eq 0 ]; then
  warn "no multi-member clusters — every window is unique at ${PID_PASS2} identity"
  warn "with few input genomes this is expected; nothing for steps 3+ to align"
fi

# --- split into per-cluster FASTAs -----------------------------------------

if [ -d "$SPLITDIR" ] && [ -n "$(ls -A "$SPLITDIR" 2>/dev/null)" ] && [ "$FORCE" != "1" ]; then
  log "skip cluster splitting (splits/ already populated)"
else
  rm -rf "$SPLITDIR"
  mkdir -p "$SPLITDIR"

  case "$SPLIT_METHOD" in
    repo)
      # getClusterSequences.sh hardcodes both the input filename and the output
      # location (its working directory), so the all_seqs file is copied in
      # under the exact name it expects -- which matters if PID_PASS2 or
      # COVERAGE were changed, since our own filename would no longer match.
      log "splitting $n_multi cluster(s) with upstream getClusterSequences.sh"
      require_file "$SRC/getClusterSequences.sh"
      cp "${P2}_all_seqs.fasta" \
         "$SPLITDIR/mmseq2_repSeqs_covmode0_PID50_cov80_all_seqs.fasta"
      cp "$COUNTS" "$SPLITDIR/"
      ( cd "$SPLITDIR" && bash "$SRC/getClusterSequences.sh" cluster_count.tsv >/dev/null )
      rm -f "$SPLITDIR/mmseq2_repSeqs_covmode0_PID50_cov80_all_seqs.fasta" \
            "$SPLITDIR/cluster_count.tsv"
      ;;
    awk)
      # Single pass over _all_seqs.fasta. A cluster marker is a header line
      # immediately followed by another header line; everything until the next
      # marker belongs to that cluster.
      log "splitting $n_multi cluster(s) with the awk splitter"
      awk -v dir="$SPLITDIR" '
        NR==FNR { want[$2] = 1; next }
        {
          if (prev ~ /^>/ && $0 ~ /^>/) {
            if (out != "") { close(out); out = "" }
            c = substr(prev, 2); sub(/[ \t].*/, "", c)
            out = (c in want) ? dir "/" c "_cluster.fasta" : ""
          } else if (prev ~ /^>/ && out != "") {
            print prev > out
            print $0  > out
          }
          prev = $0
        }
      ' "$COUNTS" "${P2}_all_seqs.fasta"
      ;;
    *)
      die "unknown SPLIT_METHOD '$SPLIT_METHOD' (expected: repo | awk)"
      ;;
  esac
fi

# --- verify the splits -----------------------------------------------------

shopt -s nullglob
splits=( "$SPLITDIR"/*_cluster.fasta )
shopt -u nullglob

[ ${#splits[@]} -eq "$n_multi" ] \
  || warn "wrote ${#splits[@]} cluster FASTAs but expected $n_multi"

# The classic symptom of the tab bug is files that exist but are empty, so
# check contents rather than just presence.
empty=0
mismatch=0
while IFS=$'\t' read -r want name; do
  f="$SPLITDIR/${name}_cluster.fasta"
  if [ ! -s "$f" ]; then
    empty=$(( empty + 1 ))
    continue
  fi
  got=$(grep -c '^>' "$f")
  [ "$got" -eq "$want" ] || mismatch=$(( mismatch + 1 ))
done < "$COUNTS"

[ "$empty" -eq 0 ] || die "$empty cluster FASTA(s) are empty or missing"
[ "$mismatch" -eq 0 ] || warn "$mismatch cluster(s) have a different member count than cluster_count.tsv"

# --- cleanup ---------------------------------------------------------------

if [ "$CLEAN_MMSEQS_TMP" = "1" ]; then
  rm -rf "$TMP1" "$TMP2"
  log "removed mmseqs tmp directories"
fi

conda_deactivate

# --- report ----------------------------------------------------------------

n_in=$(grep -c '^>' "$UNIQ")
n_clust_in=$(grep -c '^>' "$CLUSTER_INPUT")
n_rep=$(grep -c '^>' "${P1}_rep_seq.fasta")
n_total=$(cut -f1 "${P2}_cluster.tsv" | sort -u | wc -l)

echo
printf '%-42s %10s\n' STAGE SEQUENCES
printf '%-42s %10s\n' ------------------------------------------ ----------
printf '%-42s %10s\n' "unique windows from step 1"          "$n_in"
printf '%-42s %10s\n' "clustering input (both strands)"     "$n_clust_in"
printf '%-42s %10s\n' "after ${PID_PASS1} redundancy removal"   "$n_rep"

echo
printf '%-42s %10s\n' "clusters at ${PID_PASS2} identity"       "$n_total"
printf '%-42s %10s\n' "  with >= 2 members (kept)"          "$n_multi"
printf '%-42s %10s\n' "  singletons (dropped)"              "$(( n_total - n_multi ))"
printf '%-42s %10s\n' "cluster FASTAs written"              "${#splits[@]}"

echo
if [ "${#splits[@]}" -gt 0 ]; then
  biggest=$(sort -k1,1nr "$COUNTS" | head -1)
  log "largest cluster: $(echo "$biggest" | cut -f2) ($(echo "$biggest" | cut -f1) members)"
fi
log "step 3 input: $SPLITDIR/*_cluster.fasta"
log "outputs in $OUTDIR"