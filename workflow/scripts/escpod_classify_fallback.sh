#!/usr/bin/env bash
#
# Run `escpod classify` so a sample too thin to detect the move-table frame is
# classified with a known frame instead of failing the whole run.
#
# WHY THIS EXISTS
#
# `escpod classify --orientation auto` decides the frame from the data, which
# needs >= 50 informative reads and a 95% consensus. Below that it is a hard
# ERROR, not a warning:
#
#   Error: orientation check underpowered: 19 informative reads (need 50);
#          pass --orientation to override for small batches
#
# A sample only gets that shallow when a run is split many ways -- at 224
# samples per flow cell a sparse condition lands ~20 anchored reads -- and since
# every sample owes a charging table, one of them fails `rule all` for the whole
# corpus. Retrying it with the frame the rest of the run measured turns that
# into a completed (if shallow) sample.
#
# WHY A SCRIPT AND NOT A SHELL BLOCK
#
# The same reason dorado_basecall_resume.sh is a script: the logic wants a
# conditional and a second invocation with one argument changed, and inlining a
# bash function in a rule's `shell:` makes snakefmt non-convergent -- it
# re-indents the function body on every run, so `snakefmt --check` can never
# pass. Keeping it here also makes it testable without a workflow.
#
# WHAT IT WILL NOT DO
#
# It retries ONLY on the underpowered-orientation error. Any other failure --
# OOM, a missing POD5, an unreadable model, a bad reference -- propagates
# unchanged and is never retried, so a real problem is never masked by a second
# attempt that hides it. The retry is announced in the log, so a run always
# records which samples had a frame supplied rather than detected.
#
# usage: escpod_classify_fallback.sh <fallback-frame|""> <log> <classify args...>

set -uo pipefail

if [ "$#" -lt 3 ]; then
    echo "usage: $(basename "$0") <fallback|''> <log> <escpod classify args...>" >&2
    exit 2
fi

fallback=$1
log=$2
shift 2

mkdir -p "$(dirname "$log")"

if escpod classify "$@" >"$log" 2>&1; then
    exit 0
fi

# Not the error this script is here for: hand back the original failure.
if [ -z "$fallback" ] || ! grep -q "orientation check underpowered" "$log"; then
    exit 1
fi

{
    echo
    echo "NOTE: orientation detection was underpowered for this sample, so the"
    echo "      frame was SUPPLIED rather than detected: --orientation $fallback"
    echo "      (charging.orientation_fallback). A sample this shallow is"
    echo "      complete, not necessarily meaningful -- exclude it downstream on"
    echo "      a depth guard rather than reading its charging fraction."
} >>"$log"

escpod classify "$@" --orientation "$fallback" >>"$log" 2>&1
