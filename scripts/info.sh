#!/usr/bin/env bash
#
# info.sh — central configuration for the RNA_conservation pipeline.
#
# Sourced by every step script; not meant to be executed on its own.
#
# Every value uses the  : "${VAR:=default}"  form, which means "set VAR to this
# default ONLY if it isn't already set". So anything here can be overridden from
# the environment without editing this file:
#
#     FORCE=1 ./00_oneline.sh
#     WINDOW=500 OVERLAP=100 ./01_windows.sh
#
# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

# Repo root. Set to the physical path (`pwd -P`) rather than a path through the
# ~/zfs-data symlink, so nothing downstream inherits the symlink.
#
# Kept in the ${VAR:=default} form so it can still be overridden per-invocation
# without editing this file -- useful when someone else runs the pipeline from
# their own clone:   REPO=/home/alice/RNA_conservation ./00_oneline.sh
: "${REPO:=$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd -P)}"

# A hardcoded path is only as good as its accuracy, so verify it up front
# rather than letting every downstream path fail one at a time.
[ -d "$REPO" ] || { echo "info.sh: REPO does not exist: $REPO" >&2; return 1 2>/dev/null || exit 1; }

# Root for all pipeline data. Every step writes into $DATA/<stepname>/.
: "${DATA:=$REPO/data}"

# Raw, unmodified input FASTAs. Nothing ever writes here.
: "${INPUT:=$DATA/input}"

# Per-step run logs.
: "${LOGDIR:=$DATA/logs}"

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
# to $LOGDIR automatically.
#
# Renaming after a run leaves the old directory in place; the next step will
# fail with "no input in <new dir>" rather than silently reusing stale data.

: "${STEP_00:=00_oneline}"    # normalise input FASTA to one line per record
: "${STEP_01:=01_windows}"    # slice into overlapping windows, filter, dedupe
: "${STEP_02:=02_clusters}"   # reverse complement + mmseqs2 clustering
: "${STEP_03:=03_screen}"     # clustal omega + RNALalifold screening
: "${STEP_04:=04_locarna}"    # mlocarna structure prediction
: "${STEP_05:=05_eval}"       # trim, RNA-SCoRE, R-scape
: "${STEP_06:=06_homolog}"    # Infernal homolog search
: "${STEP_07:=07_compat}"     # HMMER structure/sequence compatibility

# Timestamp identifying this run; becomes the log filename prefix, e.g.
#   20260729T143518_00_oneline.log
# Each script that sources info.sh independently gets its own RUN_ID. Export it
# beforehand to make every step of one pipeline invocation share a prefix:
#   export RUN_ID=$(date +%Y%m%dT%H%M%S); ./00_oneline.sh && ./01_windows.sh
: "${RUN_ID:=$(date +%Y%m%dT%H%M%S)}"

# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------

: "${IMAGE:=dollycm/rnatools:v2.1}"

# ---------------------------------------------------------------------------
# Pipeline parameters
# ---------------------------------------------------------------------------

# Step 1 windowing. 250/75 are the values from the paper's windowing_job.sh,
# not the 100/10 shown in the README example. Step size is WINDOW - OVERLAP.
: "${WINDOW:=250}"
: "${OVERLAP:=75}"

# Steps 3-4 RNA folding temperature in degrees C. 21 = plant growth conditions
# in the original study; change this for your organism.
: "${FOLD_TEMP:=21}"

# Name of the list of clusters that passed RNALalifold screening. Written by
# step 3, read by step 4.
: "${PASSED_LIST:=RNALalifold_passedList.txt}"

# Step 2 clustering identity thresholds (two-pass).
: "${PID_PASS1:=0.95}"
: "${PID_PASS2:=0.50}"
: "${COVERAGE:=0.80}"

# --kmer-per-seq is passed to the first mmseqs pass only. The upstream README
# omits it from the second command, so we do too.
: "${KMER_PER_SEQ:=200}"

# Step 2: generate reverse-complement windows before clustering. RNA folding
# only considers the strand it is given. Set to 0 when the input is already
# strand-specific (e.g. UTRs of coding genes) -- upstream marks this optional.
: "${RUN_REVCOMP:=1}"

# Step 2 cluster splitting backend:
#   repo - upstream getClusterSequences.sh          (default)
#   awk  - single-pass equivalent, far faster at scale
# The upstream script re-scans the whole _all_seqs.fasta once per cluster, so
# its cost is O(clusters x filesize); with thousands of clusters that dominates
# the step. Both produce byte-identical output.
: "${SPLIT_METHOD:=repo}"

# Remove mmseqs2 tmp directories when the step completes (they are large).
: "${CLEAN_MMSEQS_TMP:=1}"

# Conda environments, assumed to already exist (created from the .yml files in
# step2_clustering/). Each step activates what it needs and verifies the tool is
# present before running anything.
: "${CONDA_ENV_SEQTK:=seqtk_env}"
: "${CONDA_ENV_MMSEQS:=mmseq2env}"

: "${THREADS:=8}"

# ---------------------------------------------------------------------------
# Behaviour
# ---------------------------------------------------------------------------

# Step 0 conversion backend:
#   awk  - fast, handles files with no trailing newline    (default)
#   repo - call the upstream fasta_oneLiner.sh verbatim    (for cross-checking)
: "${ONELINE_METHOD:=awk}"

# 1 = rebuild outputs that already exist, 0 = skip them.
: "${FORCE:=0}"

# Basename for step 1's combined/filtered window files.
: "${WINDOWS_PREFIX:=all_windows}"

# Sort step 1's filtered FASTAs by header. The upstream perl filters iterate
# `keys %hash`, and perl randomises hash order per process, so reruns emit the
# same records in a different order -- which also lets mmseqs2 in step 2 pick
# different cluster representatives. Default 1 for byte-reproducible runs; set
# to 0 to match upstream's unordered behaviour exactly.
: "${SORT_WINDOWS:=1}"

# Run countNs.pl to record per-genome N content alongside the windows.
: "${RUN_COUNTNS:=1}"