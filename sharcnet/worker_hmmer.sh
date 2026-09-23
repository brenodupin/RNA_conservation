#!/bin/bash

# Realign one motif around a single representative sequence with HMMER, and
# carry the representative's structure over to the new alignment -- the
# sequence/structure compatibility check of step 7 (nhmmer_local_job.sh and
# addStructToHMMERaln_job.sh upstream).
#
#   1. esl-reformat      the alignment's sequences, ungapped    <prefix>_seqs.fa
#   2. esl-alimanip +    the representative alone, its gap
#      esl-alimask       columns removed, with its structure    <prefix>_one.sto
#   3. esl-sfetch        the representative's sequence          <prefix>_one.fa
#   4. hmmbuild          a profile HMM of that one sequence     <prefix>_one.hmm
#   5. hmmsearch -A      every sequence searched with it; the
#                        ones it includes, aligned              <prefix>.sto
#   6. addStructToHMMERaln.pl
#                        the representative's structure placed
#                        on the HMM's match columns             <prefix>_withss.sto
#
# Sequences the one-sequence HMM does not include are left out of the new
# alignment: upstream's "remove sequences which do not fit".
#
# HMMER and easel run in the rnatools container, where they are not on the
# PATH, hence the two directories. addStructToHMMERaln.pl runs on the host.
# It reads the structure line of <prefix>_one.sto in one piece, so that file
# is written one line per sequence (Pfam format).
#
# <prefix>_withss.sto is written last, after checking the structure it
# carries has one character per HMM match column.
#
# Arguments:
#   container, directory, alignment, representative name, output prefix,
#   HMMER directory, easel directory, addStructToHMMERaln.pl,
#   extra hmmsearch options

set -euo pipefail

if [[ $# -ne 9 ]]; then
    echo \
        "Usage: $0 CONTAINER DIR ALIGNMENT REPRESENTATIVE PREFIX HMMER_BIN EASEL_BIN ADDSTRUCT_PL OPTIONS" \
        >&2
    exit 1
fi

rnatools=$1
dir=$2
alignment=$3
representative=$4
prefix=$5
hmmer_bin=$6
easel_bin=$7
addstruct_pl=$8
options=$9

[[ -s "$rnatools" ]] || {
    echo "Container not found: $rnatools" >&2
    exit 1
}

[[ -s "$dir/$alignment" ]] || {
    echo "Alignment not found: $dir/$alignment" >&2
    exit 1
}

[[ -s "$addstruct_pl" ]] || {
    echo "addStructToHMMERaln.pl not found: $addstruct_pl" >&2
    exit 1
}

log="$dir/${prefix}.log"

rm -f \
    "$dir/${prefix}_seqs.fa" "$dir/${prefix}_rep.txt" \
    "$dir/${prefix}_one.v0.sto" "$dir/${prefix}_one.v1.sto" "$dir/${prefix}_one.sto" \
    "$dir/${prefix}_one.fa" "$dir/${prefix}_one.hmm" \
    "$dir/${prefix}.sto" "$dir/${prefix}_hmmsearch.txt" \
    "$dir/${prefix}_withss.tmp" "$dir/${prefix}_withss.sto" \
    "$log"

in_container() {
    apptainer exec \
        --bind "$dir:/work" \
        --pwd /work \
        "$rnatools" \
        "$@" \
        >> "$log" 2>&1
}

step() {  # step <what> <command...>
    local what=$1
    shift

    if ! in_container "$@"; then
        echo "$what failed for $dir/$alignment, see $log" >&2
        exit 1
    fi
}

if ! awk -v name="$representative" '
        $0 !~ /^#/ && $0 !~ /^\/\// && NF >= 2 && $1 == name { found = 1; exit }
        END { exit !found }
    ' "$dir/$alignment"
then
    echo "Representative $representative is not a sequence of $dir/$alignment" >&2
    exit 1
fi

printf '%s\n' "$representative" > "$dir/${prefix}_rep.txt"

step esl-reformat \
    "$easel_bin/esl-reformat" -o "${prefix}_seqs.fa" fasta "$alignment"

step esl-alimanip \
    "$easel_bin/esl-alimanip" -o "${prefix}_one.v0.sto" \
    --seq-k "${prefix}_rep.txt" "$alignment"

step esl-alimask \
    "$easel_bin/esl-alimask" -o "${prefix}_one.v1.sto" \
    --gapthresh 0.9 -g "${prefix}_one.v0.sto"

step esl-reformat \
    "$easel_bin/esl-reformat" -o "${prefix}_one.sto" pfam "${prefix}_one.v1.sto"

# esl-sfetch writes the sequence to STDOUT, which the log would swallow.
if ! apptainer exec \
    --bind "$dir:/work" \
    --pwd /work \
    "$rnatools" \
    "$easel_bin/esl-sfetch" "$alignment" "$representative" \
    > "$dir/${prefix}_one.fa" 2>> "$log"
then
    echo "esl-sfetch failed for $representative in $dir/$alignment, see $log" >&2
    exit 1
fi

step hmmbuild \
    "$hmmer_bin/hmmbuild" "${prefix}_one.hmm" "${prefix}_one.fa"

# shellcheck disable=SC2086 # a list of hmmsearch options, possibly empty
step hmmsearch \
    "$hmmer_bin/hmmsearch" $options -o "${prefix}_hmmsearch.txt" \
    -A "${prefix}.sto" "${prefix}_one.hmm" "${prefix}_seqs.fa"

# hmmsearch writes no alignment at all when nothing reaches its inclusion
# threshold -- not even the representative, when its HMM is too short to call
# significant.
[[ -s "$dir/${prefix}.sto" ]] || {
    echo "hmmsearch included no sequence for $dir/$alignment" >&2
    exit 1
}

if ! perl "$addstruct_pl" "$dir/${prefix}_one.sto" "$dir/${prefix}.sto" \
    > "$dir/${prefix}_withss.tmp" 2>> "$log"
then
    echo "addStructToHMMERaln.pl failed for $dir/$alignment, see $log" >&2
    exit 1
fi

# The transfer is only right when every HMM match column (x in the RF line)
# received exactly one structure character: the HMM has one match state per
# residue of the representative, and _one.sto has one column per residue.
read -r match_columns structure_length transferred < <(
    awk '
        FNR == NR {
            if ($1 == "#=GC" && $2 == "SS_cons") one = one $3
            next
        }

        $1 == "#=GC" && $2 == "RF" {
            rf = rf $3
        }

        $1 == "#=GC" && $2 == "SS_cons" {
            ss = ss $3
        }

        END {
            x = rf
            gsub(/[^x]/, "", x)
            fits = (length(ss) == length(rf) && length(ss) > 0)
            printf "%d %d %d\n", length(x), length(one), fits
        }
    ' "$dir/${prefix}_one.sto" "$dir/${prefix}_withss.tmp"
)

if ((match_columns != structure_length || transferred != 1)); then
    echo \
        "Structure transfer does not fit for $dir/$alignment:" \
        "$match_columns match columns, representative structure $structure_length long" \
        >&2
    exit 1
fi

step esl-reformat \
    "$easel_bin/esl-reformat" -o "${prefix}_withss.sto" pfam "${prefix}_withss.tmp"

rm -f \
    "$dir/${prefix}_withss.tmp" "$dir/${prefix}_one.v0.sto" \
    "$dir/${prefix}_one.v1.sto" "$dir/${prefix}_rep.txt"
