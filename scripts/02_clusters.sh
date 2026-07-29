#!/usr/bin/env bash
#
# 02_clusters.sh — add reverse-complement windows, cluster in two mmseqs2
#                  passes, and split each cluster into its own FASTA.
#
#   in : $RNAC_DATA/$RNAC_STEP_01/<prefix>_processed_noNs_polyN_uniq.fasta
#   out: $RNAC_DATA/$RNAC_STEP_02/<prefix>_revComp.fasta
#        $RNAC_DATA/$RNAC_STEP_02/<prefix>_bothStrands.fasta        clustering input
#        $RNAC_DATA/$RNAC_STEP_02/mmseq2_covmode0_PID95_cov80_*     pass 1 (redundancy)
#        $RNAC_DATA/$RNAC_STEP_02/mmseq2_repSeqs_covmode0_PID50_cov80_*   pass 2 (clusters)
#        $RNAC_DATA/$RNAC_STEP_02/cluster_count.tsv                 members per cluster
#        $RNAC_DATA/$RNAC_STEP_02/splits/<rep>_cluster.fasta        one file per cluster
#
#   ./02_clusters.sh                     normal run
#   RNAC_FORCE=1 ./02_clusters.sh             rebuild everything
#   RNAC_RUN_REVCOMP=0 ./02_clusters.sh       strand-specific input, skip minus strand
#   RNAC_SPLIT_METHOD=awk ./02_clusters.sh    fast splitter (identical output)
#   RNAC_PID_PASS2=0.40 ./02_clusters.sh      looser cluster threshold
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

STEP="$RNAC_STEP_02"
INDIR="$RNAC_DATA/$RNAC_STEP_01"
OUTDIR="$RNAC_DATA/$STEP"
SPLITDIR="$OUTDIR/splits"
SRC="$RNAC_REPO/step2_clustering"

start_log "$STEP"

UNIQ="$INDIR/${RNAC_WINDOWS_PREFIX}_processed_noNs_polyN_uniq.fasta"
require_file "$UNIQ"
mkdir -p "$OUTDIR"

# --- validate thresholds ---------------------------------------------------

frac_ok() { awk -v v="$1" 'BEGIN{ exit !(v > 0 && v <= 1) }'; }
pct()     { awk -v v="$1" 'BEGIN{ printf "%d", v*100 + 0.5 }'; }

for v in RNAC_PID_PASS1 RNAC_PID_PASS2 RNAC_COVERAGE; do
  frac_ok "${!v}" || die "$v must be a fraction in (0,1], got '${!v}'"
done
awk -v a="$RNAC_PID_PASS2" -v b="$RNAC_PID_PASS1" 'BEGIN{ exit !(a <= b) }' \
  || die "RNAC_PID_PASS2 ($RNAC_PID_PASS2) should not exceed RNAC_PID_PASS1 ($RNAC_PID_PASS1): pass 2 clusters the survivors of pass 1 at a looser threshold"

log "pass1 min-seq-id=$RNAC_PID_PASS1  pass2 min-seq-id=$RNAC_PID_PASS2  coverage=$RNAC_COVERAGE  threads=$RNAC_THREADS"

# --- reverse complement ----------------------------------------------------
#
# seqtk emits one line per sequence by default (-l 0), which the rest of the
# pipeline depends on. Passed explicitly rather than relied upon.
#
# The trailing 'r' appended to each header is what marks a window as minus
# strand; it survives all the way to the final motif coordinates, so do not
# change it without checking step 6's coordinate handling.

REVCOMP="$OUTDIR/${RNAC_WINDOWS_PREFIX}_revComp.fasta"
BOTH="$OUTDIR/${RNAC_WINDOWS_PREFIX}_bothStrands.fasta"

if [ "$RNAC_RUN_REVCOMP" = "1" ]; then
  if [ -s "$BOTH" ] && [ "$RNAC_FORCE" != "1" ]; then
    log "skip reverse complement (output exists)"
  else
    conda_activate "$RNAC_CONDA_ENV_SEQTK" seqtk
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
  log "RNAC_RUN_REVCOMP=0 — clustering the plus strand only"
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

P1="$OUTDIR/mmseq2_covmode0_PID$(pct "$RNAC_PID_PASS1")_cov$(pct "$RNAC_COVERAGE")"
P2="$OUTDIR/mmseq2_repSeqs_covmode0_PID$(pct "$RNAC_PID_PASS2")_cov$(pct "$RNAC_COVERAGE")"
TMP1="$OUTDIR/tmp_covmode0_PID$(pct "$RNAC_PID_PASS1")_cov$(pct "$RNAC_COVERAGE")"
TMP2="$OUTDIR/tmp_repSeqs_covmode0_PID$(pct "$RNAC_PID_PASS2")_cov$(pct "$RNAC_COVERAGE")"

if [ -s "${P1}_rep_seq.fasta" ] && [ "$RNAC_FORCE" != "1" ]; then
  log "skip mmseqs pass 1 (output exists)"
else
  conda_activate "$RNAC_CONDA_ENV_MMSEQS" mmseqs
  log "mmseqs pass 1: collapsing windows >= ${RNAC_PID_PASS1} identical"
  mmseqs easy-cluster "$CLUSTER_INPUT" "$P1" "$TMP1" \
    -c "$RNAC_COVERAGE" --threads "$RNAC_THREADS" --kmer-per-seq "$RNAC_KMER_PER_SEQ" \
    --min-seq-id "$RNAC_PID_PASS1" --cov-mode 0 --filter-hits 1
  require_file "${P1}_rep_seq.fasta"
fi

if [ -s "${P2}_cluster.tsv" ] && [ "$RNAC_FORCE" != "1" ]; then
  log "skip mmseqs pass 2 (output exists)"
else
  conda_activate "$RNAC_CONDA_ENV_MMSEQS" mmseqs
  log "mmseqs pass 2: clustering representatives at >= ${RNAC_PID_PASS2} identity"
  mmseqs easy-cluster "${P1}_rep_seq.fasta" "$P2" "$TMP2" \
    -c "$RNAC_COVERAGE" --threads "$RNAC_THREADS" \
    --min-seq-id "$RNAC_PID_PASS2" --cov-mode 0 --filter-hits 1
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

if [ -s "$COUNTS" ] && [ "$RNAC_FORCE" != "1" ]; then
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
  warn "no multi-member clusters — every window is unique at ${RNAC_PID_PASS2} identity"
  warn "with few input genomes this is expected; nothing for steps 3+ to align"
fi

# --- split into per-cluster FASTAs -----------------------------------------

count_splits() {
  find "$SPLITDIR" -maxdepth 1 -type f -name '*_cluster.fasta' 2>/dev/null | wc -l
}

n_have=0
[ -d "$SPLITDIR" ] && n_have=$(count_splits)

# Completeness is a count comparison, not "is the directory non-empty". An
# interrupted split leaves a populated directory that looks finished, and
# skipping it then fails verification thousands of clusters later.
if [ "$n_multi" -gt 0 ] && [ "$n_have" -eq "$n_multi" ] && [ "$RNAC_FORCE" != "1" ]; then
  log "skip cluster splitting ($n_have/$n_multi already present)"
else
  if [ "$n_have" -gt 0 ]; then
    log "splits/ holds $n_have of $n_multi expected files — rebuilding from scratch"
  fi
  rm -rf "$SPLITDIR"
  mkdir -p "$SPLITDIR"

  case "$RNAC_SPLIT_METHOD" in
    repo)
      # getClusterSequences.sh hardcodes both the input filename and the output
      # location (its working directory), so the all_seqs file is copied in
      # under the exact name it expects -- which matters if RNAC_PID_PASS2 or
      # RNAC_COVERAGE were changed, since our own filename would no longer match.
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
      die "unknown RNAC_SPLIT_METHOD '$RNAC_SPLIT_METHOD' (expected: repo | awk)"
      ;;
  esac
fi

# --- verify the splits -----------------------------------------------------
#
# All three checks are single-pass. The previous version ran one grep per
# cluster, which cost ~1 minute at 44k clusters and ran on every invocation,
# including ones that skipped the split entirely.

n_have=$(count_splits)
[ "$n_have" -eq "$n_multi" ] \
  || die "wrote $n_have cluster FASTAs but expected $n_multi"

# Empty-but-present files are the classic symptom of cluster_count.tsv losing
# its tabs, since getClusterSequences.sh then reads an empty cluster name.
n_empty=$(find "$SPLITDIR" -maxdepth 1 -type f -name '*_cluster.fasta' -empty | wc -l)
[ "$n_empty" -eq 0 ] \
  || die "$n_empty cluster FASTA(s) are empty — check that cluster_count.tsv is tab-separated"

# Member counts, in one grep pass over every file rather than one grep each.
# -H forces the filename prefix even when xargs hands grep a single file.
if [ "$n_have" -gt 0 ]; then
  mismatch=$(
    find "$SPLITDIR" -maxdepth 1 -type f -name '*_cluster.fasta' -print0 \
      | xargs -0 grep -cH '^>' \
      | sed 's/_cluster\.fasta:/\t/' \
      | awk -F'\t' '
          NR==FNR { want[$2] = $1; next }
          { name = $1; sub(/^.*\//, "", name); if (want[name] != $2) n++ }
          END { print n+0 }
        ' "$COUNTS" -
  )
  [ "$mismatch" -eq 0 ] \
    || warn "$mismatch cluster(s) have a member count differing from cluster_count.tsv"
fi

# --- cleanup ---------------------------------------------------------------

if [ "$RNAC_CLEAN_MMSEQS_TMP" = "1" ]; then
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
printf '%-42s %10s\n' "after ${RNAC_PID_PASS1} redundancy removal"   "$n_rep"

echo
printf '%-42s %10s\n' "clusters at ${RNAC_PID_PASS2} identity"       "$n_total"
printf '%-42s %10s\n' "  with >= 2 members (kept)"          "$n_multi"
printf '%-42s %10s\n' "  singletons (dropped)"              "$(( n_total - n_multi ))"
printf '%-42s %10s\n' "cluster FASTAs written"              "$n_have"

echo
if [ "$n_have" -gt 0 ]; then
  # A `sort | head -1` here dies with SIGPIPE (exit 141) once cluster_count.tsv
  # is large enough that sort is still writing when head closes the pipe --
  # which pipefail then propagates and set -e turns into a failed step. One awk
  # pass has no pipe to break, and is O(n) rather than O(n log n).
  #
  # The size buckets matter for planning: RNA-SCoRE ranks an alignment High at
  # >=10 sequences and Mid at 7-9, so clusters below 7 cannot produce a usable
  # result no matter how much time step 3 and step 4 spend on them.
  awk -F'\t' '
    { n = $1
      if (n > max) { max = n; maxname = $2 }
      total++
      if (n == 2)      b2++
      else if (n <= 4) b34++
      else if (n <= 6) b56++
      else if (n <= 9) b79++
      else             b10++
    }
    END {
      printf "%-42s %10s\n", "largest cluster", max " (" maxname ")"
      printf "\n%-42s %10s %9s\n", "CLUSTER SIZE", "COUNT", "SHARE"
      printf "%-42s %10s %9s\n", "------------------------------------------", "----------", "---------"
      printf "%-42s %10d %8.1f%%\n", "2 members",              b2+0,  100*(b2+0)/total
      printf "%-42s %10d %8.1f%%\n", "3-4 members",            b34+0, 100*(b34+0)/total
      printf "%-42s %10d %8.1f%%\n", "5-6 members",            b56+0, 100*(b56+0)/total
      printf "%-42s %10d %8.1f%%\n", "7-9 members (RNA-SCoRE Mid)",  b79+0, 100*(b79+0)/total
      printf "%-42s %10d %8.1f%%\n", ">=10 members (RNA-SCoRE High)", b10+0, 100*(b10+0)/total
      printf "\n%-42s %10d\n", "usable for a Mid/High rank (>=7)", (b79+0)+(b10+0)
    }
  ' "$COUNTS"
fi
echo
log "step 3 input: $SPLITDIR/*_cluster.fasta"
log "outputs in $OUTDIR"