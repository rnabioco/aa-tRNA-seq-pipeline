#!/usr/bin/env bash
#
# Run `dorado basecaller` so a job killed on wall clock loses a tail of the
# basecall instead of all of it.
#
# dorado already knows how to resume; the whole job of this wrapper is to keep
# the partial alive long enough for `--resume-from` to be handed it, and to
# refuse a partial that would be wrong to resume from.
#
# WHY THE PARTIAL CANNOT BE THE RULE'S OUTPUT
#
# Snakemake deletes a failed job's outputs ("Removing output files of failed
# job ... since they might be corrupted"), and the basecall is additionally
# maybe_temp(tier="basecall"). So the one artifact worth keeping is destroyed by
# the failure path that creates it. The partial is therefore written to a
# sibling dot-file Snakemake knows nothing about, and renamed onto the real
# output only once dorado exits 0. `clean` still reclaims it: that rule removes
# the whole bam/rebasecall{,_run} directory.
#
# WHAT --resume-from ACTUALLY DOES  (dorado ResumeLoader.cpp, basecaller.cpp)
#
#   * `copy_completed_reads()` pushes every record of the resume file straight
#     into the writer, so the NEW output is self-contained -- the two files are
#     never concatenated, and no read is counted twice.
#   * The read ids it saw are then excluded from the basecall, using the `pi`
#     parent tag for split reads so a split child resumes as its parent.
#   * The record loop is wrapped in `try { ... } catch { /* end of properly
#     formatted records */ }`, i.e. a half-written trailing record is expected
#     and ignored. A BAM cut off mid-write by a kill is the designed input.
#   * It refuses to run with `--output-dir`, so the stdout redirect stays.
#   * A missing resume file is a hard error, hence the branch below.
#   * dorado re-parses the checkpoint's `@PG ID:basecaller CL:` and ERRORS if
#     the models differ. That is a good backstop, but leaning on it would wedge
#     the rule: after a model change every retry would fail on the same stale
#     file. The argv is compared here first, so a superseded checkpoint is
#     dropped and the basecall simply starts over.
#
# usage: dorado_basecall_resume.sh <final_bam> <args to `dorado basecaller`...>

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: $(basename "$0") <final_bam> <dorado basecaller args...>" >&2
    exit 2
fi

final=$1
shift

dir=$(dirname "$final")
base=$(basename "$final")
part="$dir/.$base.partial" # what this attempt writes
prev="$dir/.$base.resume"  # what dorado reads; never also written
args="$dir/.$base.args"    # the argv that produced them

mkdir -p "$dir"

# The argv is the identity of a basecall: model, modified-base models, the
# --read-ids path, the POD5 root. A checkpoint written under a different one
# describes different work and must not be continued.
this_args=$(printf '%s\n' "$@")

resumable() {
    # dorado needs a header carrying `@PG ID:basecaller ... CL:` -- it re-parses
    # that command line to check the model, and dies if it is absent. Checking
    # here means a checkpoint truncated before its header ever landed is dropped
    # quietly rather than failing every retry. Read the header into a variable
    # instead of piping it: `grep -q` exits early, and the SIGPIPE that sends to
    # samtools would trip `pipefail` on a perfectly good file.
    local hdr
    hdr=$(samtools view -H "$1" 2>/dev/null) || return 1
    grep -q '^@PG.*ID:basecaller.*CL:' <<<"$hdr"
}

# Choose which checkpoint to continue from. Normally that is this rule's last
# attempt ($part) -- but dorado copies the resumed records through BEFORE it
# basecalls anything new, so an attempt killed inside that copy leaves a partial
# SHORTER than the file it resumed from. Taking the newest blindly would throw
# away the longer one, so take the bigger.
#
# Size stands in for record count: both files are the same basecall of the same
# reads, so they compress alike. It is not exact -- a resumed attempt's header
# carries an extra `--resume-from <path>` in its @PG CL, worth some tens of
# bytes -- so this only decides correctly when the two differ by more than a
# header. At the scale that matters (a partial flowcell against a nearly empty
# one) they differ by gigabytes, and when they don't, neither is worth keeping.
best=
for cand in "$part" "$prev"; do
    if resumable "$cand" &&
        { [ -z "$best" ] || [ "$(stat -c%s "$cand")" -gt "$(stat -c%s "$best")" ]; }; then
        best=$cand
    fi
done

resume=()
if [ -n "$best" ]; then
    if [ -r "$args" ] && [ "$(cat "$args")" = "$this_args" ]; then
        if [ "$best" = "$part" ]; then
            # "Do not reuse the filenames for --resume-from and the new output.
            #  If they are the same then the interrupted file will be deleted
            #  when Dorado is launched and the previous work will be lost."
            mv -f "$part" "$prev"
        fi
        resume=(--resume-from "$prev")
        echo "resume: continuing from $(basename "$prev") ($(stat -c%s "$prev") bytes)" >&2
    else
        echo "resume: discarding checkpoint -- basecall arguments changed since it was written" >&2
        best=
    fi
fi

if [ -z "$best" ]; then
    rm -f "$prev"
fi
rm -f "$part"

# Written before dorado starts, so it always describes the partial on disk.
printf '%s\n' "$this_args" >"$args"

dorado basecaller "${resume[@]}" "$@" >"$part"

# Only now is there something worth handing to Snakemake. Same directory, so the
# rename is atomic: the output never exists half-written.
mv -f "$part" "$final"
rm -f "$prev" "$args"
