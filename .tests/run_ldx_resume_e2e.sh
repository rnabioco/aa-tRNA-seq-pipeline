#!/usr/bin/env bash
#
# End-to-end test of basecall resume (issue #149) against the committed LDX
# fixture. NEEDS A GPU, so it is not part of CI -- run it by hand, inside an
# allocation, before trusting a change to dorado_basecall_resume.sh:
#
#     pixi run test-ldx-resume
#
# WARNING: this DESTROYS and rebuilds .tests/outputs-ldx. It leaves that tree in
# a good, complete state, but anything you had there is gone.
#
# WHY THE TEST IS SHAPED LIKE THIS
#
# The fixture's 323 reads basecall in far less time than dorado takes to start,
# so no signal can be made to land mid-write. The two claims are therefore
# proved apart:
#
#   PHASE A   baseline run, snapshotted.
#   PHASE A2  CONTROL: re-run with no checkpoint anywhere. Whatever differs
#             between A and A2 is the pipeline's own nondeterminism floor and
#             cannot be laid at the resume path's door. This is measured, not
#             assumed: bwa picks arbitrarily among equal-scoring near-identical
#             tRNA isodecoders, which moves per-reference summary rows around
#             (a read lands on tRNA-Gln-CTG-1-1 in one run and tRNA-Gln-TTG-1-1
#             in the next, with an IDENTICAL charging score).
#   PHASE B   a genuine kill: SIGTERM to the whole process GROUP, which is what
#             Slurm does to a step that overruns. Proves Snakemake's failure
#             path takes {output} and leaves the checkpoint -- the blocker in
#             #149. The partial it leaves is empty, because the signal can only
#             land during dorado's startup on a fixture this small.
#   PHASE C   so the resume itself is exercised from a checkpoint holding real
#             dorado records: the Phase A uBAM truncated mid-record, byte for
#             byte what a kill leaves, carrying the wrapper's OWN .args from
#             Phase B so the argv match it performs is genuine.
#   PHASE D   resume, then compare everything. An artifact that moved in the
#             control is reported but cannot fail the run; every artifact that
#             held still must be identical.
#
# TWO TRAPS THIS SCRIPT EXISTS TO AVOID, both of which produced convincing false
# results while it was being written:
#
#   * SIGTERM to snakemake ALONE is not what Slurm does. Snakemake treats it as
#     a graceful stop and WAITS for the running job, so the basecall completes
#     and the checkpoint is cleaned up. Signal the process group.
#   * Deleting a rule's output is NOT enough to make it re-run. If the consumers
#     are up to date the DAG has nothing to do and snakemake exits 0 having done
#     nothing -- and every downstream comparison then passes against artifacts
#     nothing regenerated. Hence --forcerun, and the explicit "actually rebuilt"
#     assertions.

set -uo pipefail

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO_ROOT" || exit 1

CONFIG=config/config-ldx-test.yml
OUT=.tests/outputs-ldx
SNAP=${LDX_RESUME_SNAPDIR:-$(mktemp -d -t ldx-resume-XXXXXX)}
mkdir -p "$SNAP"

fail=0
step() { printf '\n=== %s ===\n' "$*"; }
check() {
    if [ "$2" = "$3" ]; then
        printf 'PASS  %-56s %s\n' "$1" "$2"
    else
        printf 'FAIL  %-56s got=%s want=%s\n' "$1" "$2" "$3"
        fail=1
    fi
}
note() { printf 'NOTE  %s\n' "$*"; }
die() { printf '\n%s\n=== RESULT: FAIL ===\n' "$*"; exit 1; }

# Read ids + sequences, sorted: BAM checksums differ on @PG lines alone.
ids_seqs() { samtools view "$1" | awk '{print $1"\t"$10}' | sort; }

snapshot() {
    local tag=$1 b s t
    ids_seqs "$UBAM" >"$SNAP/$tag.ubam.idseq"
    for b in "$OUT"/bam/final/*/*.bam; do
        s=$(basename "$b" .bam)
        samtools view "$b" | awk '{print $1}' | sort >"$SNAP/$tag.final.$s.ids"
        samtools view "$b" | grep -o 'cl:i:[0-9]*' | sort | uniq -c >"$SNAP/$tag.final.$s.cl"
    done
    for t in "$OUT"/summary/tables/*/*.tsv.gz; do
        zcat "$t" | sort >"$SNAP/$tag.table.$(basename "$t" .tsv.gz)"
    done
}
same() { diff -q "$1" "$2" >/dev/null 2>&1 && echo same || echo differ; }

step "environment"
echo "repo:      $REPO_ROOT ($(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout'))"
echo "snapshots: $SNAP"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -4

# The Snakefile's onstart is what puts the pinned tools on PATH inside a rule.
# This is only so the banner can report a version, and so the precondition check
# below can see the binary.
DORADO_VERSION=$(grep -Po '^dorado_version:\s*"?\K[^"\s]+' config/config-base.yml)
export PATH="$REPO_ROOT/resources/tools/dorado/$DORADO_VERSION/bin:$PATH"
dorado --version 2>&1 | head -2

# A worktree or fresh clone carries the COMMITTED half of resources/ only -- the
# charging and demux bundles are in git, the basecalling models are not
# (`pixi run setup` fetches those). Discovering that costs a GPU allocation.
step "preconditions"
BASE_MODEL=$(grep -Po '^base_calling_model:\s*"\K[^"]+' config/config-base.yml)
for p in "$BASE_MODEL" "resources/tools/dorado/$DORADO_VERSION/bin/dorado" \
    .tests/fixtures/ldx-demux/run/pod5 "$CONFIG"; do
    [ -e "$p" ] && printf 'ok    %s\n' "$p" || { printf 'MISSING %s\n' "$p"; fail=1; }
done
for c in snakemake samtools; do
    command -v "$c" >/dev/null && printf 'ok    %s\n' "$c" || { printf 'MISSING %s on PATH\n' "$c"; fail=1; }
done
[ "$fail" -ne 0 ] && die "preconditions unmet -- run under 'pixi run', and 'pixi run setup' for the models"

###########################################################################
step "PHASE A -- baseline"
###########################################################################
rm -rf "$OUT"
snakemake --configfile="$CONFIG" --cores 8 --rerun-incomplete >"$SNAP/phaseA.log" 2>&1
rc=$?
echo "snakemake rc=$rc"
[ "$rc" -ne 0 ] && { tail -30 "$SNAP/phaseA.log"; die "PHASE A FAILED -- nothing below is meaningful"; }

UBAM=$(find "$OUT/bam/rebasecall_run" -name '*.rbc.bam' | head -1)
[ -n "$UBAM" ] || die "no run-level uBAM produced"
RUN_DIR=$(dirname "$UBAM"); RUN_BASE=$(basename "$UBAM")
PART="$RUN_DIR/.$RUN_BASE.partial"
PREV="$RUN_DIR/.$RUN_BASE.resume"
ARGS="$RUN_DIR/.$RUN_BASE.args"
echo "run uBAM: $UBAM"
cp "$UBAM" "$SNAP/baseline.rbc.bam"
snapshot A
BASE_N=$(wc -l <"$SNAP/A.ubam.idseq")
echo "baseline uBAM records: $BASE_N"
check "phase A left no checkpoint" \
    "$([ -e "$PART" ] || [ -e "$PREV" ] || [ -e "$ARGS" ] && echo yes || echo no)" "no"

###########################################################################
step "PHASE A2 -- CONTROL: plain re-run, no checkpoint, to measure the floor"
###########################################################################
rm -f "$UBAM" "$PART" "$PREV" "$ARGS"
snakemake --configfile="$CONFIG" --cores 8 --rerun-incomplete \
    --forcerun rebasecall_ldx_run >"$SNAP/phaseA2.log" 2>&1
rc=$?
echo "snakemake rc=$rc"
check "control run completed" "$rc" "0"
check "control actually rebuilt the uBAM" "$([ -s "$UBAM" ] && echo yes || echo no)" "yes"
[ -s "$UBAM" ] || die "control did not run -- every comparison below would be vacuous"
snapshot A2

check "CONTROL: uBAM identical across a plain re-run" \
    "$(same "$SNAP/A.ubam.idseq" "$SNAP/A2.ubam.idseq")" "same"

UNSTABLE=""
for f in "$SNAP"/A.table.* "$SNAP"/A.final.*; do
    key=${f#"$SNAP"/A.}
    if [ "$(same "$f" "$SNAP/A2.$key")" = "differ" ]; then
        UNSTABLE="$UNSTABLE $key"
        note "unstable in control: $key"
    fi
done
[ -z "$UNSTABLE" ] && note "control is fully reproducible -- every artifact must match after resume"
is_unstable() { case " $UNSTABLE " in *" $1 "*) return 0;; *) return 1;; esac; }

###########################################################################
step "PHASE B -- genuine kill: SIGTERM to the process GROUP"
###########################################################################
# setsid puts the run in its own process group, so signalling the group can
# never touch anything else on the node -- which `pkill dorado` very much could,
# on a shared GPU node running someone's real basecalls.
rm -f "$UBAM" "$PART" "$PREV" "$ARGS"
setsid snakemake --configfile="$CONFIG" --cores 8 --rerun-incomplete \
    "$UBAM" >"$SNAP/phaseB.log" 2>&1 &
SNAKE=$!

# .partial is created by the wrapper's redirect at the moment dorado launches,
# just after .args is written. Waiting on it means both exist and we are inside
# the basecall -- deterministic, where a wall-clock guess is not.
for _ in $(seq 1 1200); do
    [ -e "$PART" ] && break
    kill -0 "$SNAKE" 2>/dev/null || break
    sleep 0.1
done
if [ -e "$PART" ]; then
    echo "phase B: SIGTERM to process group $SNAKE, dorado in flight"
    kill -TERM -"$SNAKE" 2>/dev/null
    sleep 2
    kill -KILL -"$SNAKE" 2>/dev/null
else
    echo "phase B: .partial never appeared -- the rule did not reach dorado"
fi
wait "$SNAKE" 2>/dev/null
echo "snakemake rc=$? (non-zero is the point)"

check "KILL: rule produced no {output}"     "$([ -e "$UBAM" ] && echo no || echo yes)" "yes"
check "KILL: .args SURVIVED the failure"    "$([ -e "$ARGS" ] && echo yes || echo no)" "yes"
check "KILL: .partial SURVIVED the failure" "$([ -e "$PART" ] && echo yes || echo no)" "yes"
[ -e "$PART" ] && echo "      .partial: $(stat -c%s "$PART") bytes, $(samtools view "$PART" 2>/dev/null | wc -l) complete records"

###########################################################################
step "PHASE C -- plant a checkpoint with real records, keep the real .args"
###########################################################################
[ -e "$ARGS" ] || die "FAIL  no .args survived phase B -- the resume path cannot be exercised"
SZ=$(stat -c%s "$SNAP/baseline.rbc.bam")
PLANTED=0
for frac in 90 75 60 50 40; do
    head -c $((SZ * frac / 100)) "$SNAP/baseline.rbc.bam" >"$SNAP/trunc.bam"
    n=$(samtools view "$SNAP/trunc.bam" 2>/dev/null | wc -l)
    echo "  truncate to ${frac}% -> $n complete records"
    if [ "$n" -ge 2 ] && [ "$n" -lt "$BASE_N" ]; then
        cp "$SNAP/trunc.bam" "$PART"; PLANTED=$n; break
    fi
done
check "planted a partial checkpoint with records" "$([ "$PLANTED" -gt 0 ] && echo yes || echo no)" "yes"
[ "$PLANTED" -gt 0 ] || die "could not build a truncation holding complete records"
echo "checkpoint holds $PLANTED of $BASE_N records"
echo "header readable: $(samtools view -H "$PART" >/dev/null 2>&1 && echo yes || echo no)"
echo "body truncated:  $(samtools view -c "$PART" >/dev/null 2>&1 && echo no || echo yes)"
rm -f "$UBAM"

###########################################################################
step "PHASE D -- resume, then compare against the control"
###########################################################################
snakemake --configfile="$CONFIG" --cores 8 --rerun-incomplete \
    --forcerun rebasecall_ldx_run >"$SNAP/phaseD.log" 2>&1
rc=$?
echo "snakemake rc=$rc"
[ "$rc" -ne 0 ] && tail -30 "$SNAP/phaseD.log"
check "phase D pipeline completed" "$rc" "0"
check "phase D actually rebuilt the uBAM" "$([ -s "$UBAM" ] && echo yes || echo no)" "yes"
[ -s "$UBAM" ] || die "the rule did not run -- the planted checkpoint was never read"

echo "--- what the wrapper and dorado said ---"
grep -h "resume:" "$SNAP/phaseD.log" "$OUT"/logs/rebasecall_ldx_run/* 2>/dev/null | head -3
CARRIED=$(grep -ho "> [0-9]* original read ids found in resume file" \
    "$SNAP/phaseD.log" "$OUT"/logs/rebasecall_ldx_run/* 2>/dev/null | grep -o '[0-9]*' | head -1)
echo "dorado carried through: ${CARRIED:-<not logged>} read ids"
check "dorado actually resumed"       "$([ "${CARRIED:-0}" -gt 0 ] 2>/dev/null && echo yes || echo no)" "yes"
check "carried count matches planted" "${CARRIED:-0}" "$PLANTED"

snapshot D
echo "--- uBAM equivalence ---"
check "resumed uBAM record count" "$(wc -l <"$SNAP/D.ubam.idseq")" "$BASE_N"
check "no duplicate read ids"     "$(cut -f1 "$SNAP/D.ubam.idseq" | uniq -d | wc -l)" "0"
check "uBAM ids+sequences identical to baseline" "$(same "$SNAP/A.ubam.idseq" "$SNAP/D.ubam.idseq")" "same"

# These records bypass the basecall pipeline entirely -- ResumeLoader pushes
# them straight into the writer -- so this is where a lossy copy-through shows.
echo "--- modbase tags and move tables survived the copy-through ---"
for tag in MM:Z: ML:B: mv:B:; do
    check "reads carrying ${tag%%:*} tags" \
        "$(samtools view "$UBAM" | grep -c "$tag")" \
        "$(samtools view "$SNAP/baseline.rbc.bam" | grep -c "$tag")"
done

echo "--- downstream: stable artifacts must match, unstable ones are reported ---"
for f in "$SNAP"/A.table.* "$SNAP"/A.final.*; do
    key=${f#"$SNAP"/A.}
    verdict=$(same "$f" "$SNAP/D.$key")
    if is_unstable "$key"; then
        note "$key: $verdict (moves in the control too -- not attributable to resume)"
    else
        check "$key" "$verdict" "same"
    fi
done

check "phase D left no checkpoint" \
    "$([ -e "$PART" ] || [ -e "$PREV" ] || [ -e "$ARGS" ] && echo yes || echo no)" "no"

printf '\n=== RESULT: %s ===\n' "$([ "$fail" -eq 0 ] && echo PASS || echo FAIL)"
printf 'snapshots and per-phase logs: %s\n' "$SNAP"
exit "$fail"
