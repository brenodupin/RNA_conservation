#!/usr/bin/env bash
#
# info.sh — central configuration for the RNA_conservation pipeline.
#
# Sourced by every step script; not meant to be executed on its own.
#
# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# Repo root. Physical path (`pwd -P`), not the route through the ~/zfs-data
# symlink, so nothing downstream inherits the symlink -- docker resolves bind
# mounts on the host and symlink/physical mismatches there fail confusingly.
#
# Hardcoded with no override at all: this is a fixed property of the machine,
# and it is the one value that must never be wrong.
_repo_from_env="${RNAC_REPO:-}"
RNAC_REPO=$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd -P)
if [ -n "$_repo_from_env" ] && [ "$_repo_from_env" != "$RNAC_REPO" ]; then
  echo "info.sh: ignoring RNAC_REPO=$_repo_from_env from the environment (RNAC_REPO is hardcoded here)" >&2
fi
unset _repo_from_env

[ -d "$RNAC_REPO" ] || {
  echo "info.sh: RNAC_REPO does not exist: $RNAC_REPO" >&2
  return 1 2>/dev/null || exit 1
}

# Root for all pipeline data. Every step writes into $RNAC_DATA/<stepname>/.
# Overridable so you can point a test run at a scratch tree:
#     RNAC_DATA=$HOME/sandbox/rna_test/data ./00_oneline.sh
: "${RNAC_DATA:=$RNAC_REPO/data}"

# Raw, unmodified input FASTAs. Nothing ever writes here.
: "${RNAC_INPUT:=$RNAC_DATA/input}"

# Per-step run logs.
: "${RNAC_LOGDIR:=$RNAC_DATA/logs}"

# Timestamp identifying this run; becomes the log filename prefix, e.g.
#   20260729T143518_00_oneline.log
# Each script gets its own unless you export one to group a whole invocation:
#   export RNAC_RUN_ID=$(date +%Y%m%dT%H%M%S); ./00_oneline.sh && ./01_windows.sh
: "${RNAC_RUN_ID:=$(date +%Y%m%dT%H%M%S)}"

# ---------------------------------------------------------------------------
# Step directories
# ---------------------------------------------------------------------------
#
# Each step reads the previous step's output directory, so these names are
# shared between scripts -- 01_windows.sh needs to know what 00_oneline.sh
# called its output. Defining them once here means renaming a step is a one-line
# change instead of a hunt through every downstream script.
#
# Keyed by number rather than by description: the number is the step's fixed
# position in the pipeline, whereas a descriptive key would go stale the moment
# you renamed the directory it points at.
#
# The same name is used for that step's log file, so a rename carries through
# to $RNAC_LOGDIR automatically.
#
# Renaming after a run leaves the old directory in place; the next step will
# fail with "no input in <new dir>" rather than silently reusing stale data.

: "${RNAC_STEP_00:=00_oneline}"    # normalise input FASTA to one line per record
: "${RNAC_STEP_01:=01_windows}"    # slice into overlapping windows, filter, dedupe
: "${RNAC_STEP_02:=02_clusters}"   # reverse complement + mmseqs2 clustering
: "${RNAC_STEP_03:=03_screen}"     # clustal omega + RNALalifold screening
: "${RNAC_STEP_04:=04_locarna}"    # mlocarna structure prediction
: "${RNAC_STEP_05:=05_eval}"       # trim, RNA-SCoRE, R-scape
: "${RNAC_STEP_06:=06_homolog}"    # Infernal homolog search
: "${RNAC_STEP_07:=07_compat}"     # HMMER structure/sequence compatibility

# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------

: "${RNAC_IMAGE:=dollycm/rnatools:v2.1}"

# ---------------------------------------------------------------------------
# Pipeline parameters
# ---------------------------------------------------------------------------

# Step 1 windowing. 250/75 are the values from the paper's windowing_job.sh,
# not the 100/10 shown in the README example. Stride is WINDOW - OVERLAP.
# OVERLAP must stay below WINDOW or createWindows.pl loops forever; 01_windows.sh
# checks this before invoking it.
: "${RNAC_WINDOW:=250}"
: "${RNAC_OVERLAP:=75}"

# Steps 3-4 RNA folding temperature in degrees C. 21 = plant growth conditions
# in the original study; change this for your organism.
: "${RNAC_FOLD_TEMP:=21}"

# Name of the list of clusters that passed RNALalifold screening. Written by
# step 3, read by step 4.
: "${RNAC_PASSED_LIST:=RNALalifold_passedList.txt}"

# Step 2 clustering identity thresholds (two-pass).
: "${RNAC_PID_PASS1:=0.95}"
: "${RNAC_PID_PASS2:=0.50}"
: "${RNAC_COVERAGE:=0.80}"

# --kmer-per-seq is passed to the first mmseqs pass only. The upstream README
# omits it from the second command, so we do too.
: "${RNAC_KMER_PER_SEQ:=200}"

: "${RNAC_THREADS:=8}"

# Clusters processed concurrently in step 3 (and later per-cluster steps).
# The work is independent per cluster and dominated by container startup rather
# than compute, so this scales close to linearly until the docker daemon or the
# filesystem saturates. With ~44k clusters the difference is hours.
#
# Distinct from RNAC_THREADS, which is threads *within* one mmseqs invocation.
: "${RNAC_JOBS:=4}"


# ---------------------------------------------------------------------------
# Behaviour
# ---------------------------------------------------------------------------

# Step 0 conversion backend:
#   awk  - fast, handles files with no trailing newline    (default)
#   repo - call the upstream fasta_oneLiner.sh verbatim    (for cross-checking)
# Upstream drops the final line of any file that does not end in a newline, so
# awk is the safer default; 00_oneline.sh compares residue counts either way.
: "${RNAC_ONELINE_METHOD:=awk}"

# 1 = rebuild outputs that already exist, 0 = skip them.
: "${RNAC_FORCE:=0}"

# Basename for step 1's combined/filtered window files.
: "${RNAC_WINDOWS_PREFIX:=all_windows}"

# Sort step 1's filtered FASTAs by header. The upstream perl filters iterate
# `keys %hash`, and perl randomises hash order per process, so reruns emit the
# same records in a different order -- which also lets mmseqs2 in step 2 pick
# different cluster representatives. Default 1 for byte-reproducible runs; set
# to 0 to match upstream's unordered behaviour exactly.
: "${RNAC_SORT_WINDOWS:=1}"

# Run countNs.pl to record per-genome N content alongside the windows.
: "${RNAC_RUN_COUNTNS:=1}"

# Step 2: generate reverse-complement windows before clustering. RNA folding
# only considers the strand it is given. Set to 0 when the input is already
# strand-specific (e.g. UTRs of coding genes) -- upstream marks this optional.
: "${RNAC_RUN_REVCOMP:=1}"

# Step 2 cluster splitting backend:
#   repo - upstream getClusterSequences.sh          (default)
#   awk  - single-pass equivalent, far faster at scale
# The upstream script re-scans the whole _all_seqs.fasta once per cluster, so
# its cost is O(clusters x filesize); with thousands of clusters that dominates
# the step. Both produce byte-identical output.
: "${RNAC_SPLIT_METHOD:=repo}"

# Remove mmseqs2 tmp directories when the step completes (they are large).
: "${RNAC_CLEAN_MMSEQS_TMP:=1}"

# Conda environments, assumed to already exist (created from the .yml files in
# step2_clustering/). Each step activates what it needs and verifies the tool is
# present before running anything.
: "${RNAC_CONDA_ENV_SEQTK:=seqtk_env}"
: "${RNAC_CONDA_ENV_MMSEQS:=mmseq2env}"

# ---------------------------------------------------------------------------
# Override visibility
# ---------------------------------------------------------------------------
#
# Collect any RNAC_* variables present in the environment so start_log can
# report them. Bash cannot distinguish a per-command prefix from a stale export,
# so the best we can do is make whatever is in effect visible in every log.

mapfile -t RNAC_OVERRIDES < <(compgen -e 2>/dev/null | grep '^RNAC_' | grep -v '^RNAC_REPO$' | LC_ALL=C sort || true)