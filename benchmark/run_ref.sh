#!/usr/bin/env bash
#
# Run the pipeline at a specific git ref into an isolated output directory and
# fingerprint the result. Used to produce the "baseline" (pre-migration) and
# "candidate" (post-migration) snapshots that compare.py then diffs.
#
# The ref is checked out into a throwaway git worktree, so your working tree is
# untouched and the OLD code runs the OLD pipeline. Fingerprinting always uses
# the CURRENT checkout's benchmark/ scripts (the old ref predates them).
#
# Usage:
#   benchmark/run_ref.sh --ref REF --config CONFIG --label NAME [options]
#
# Options:
#   --ref REF          git ref to run (branch, tag, or SHA)              [required]
#   --config FILE      snakemake --configfile                            [required]
#   --label NAME       snapshot name (dir under benchmark/snapshots/)    [required]
#   --profile DIR      snakemake --profile (e.g. cluster/slurm) for GPU
#   --cores N          snakemake --cores N (local run; default 4 if no --profile)
#   --setup            run `pixi run setup` in the worktree first (per-ref tools)
#   --targets "..."    extra snakemake targets/args (default: full `all`)
#   --keep-output      keep the run's output_directory (default: keep)
#   --keep-worktree    do not remove the git worktree afterward
#
# Example — generate both sides of a migration comparison:
#   benchmark/run_ref.sh --ref 403e755 --config config/config-test.yml \
#       --label old --profile cluster/slurm --setup
#   benchmark/run_ref.sh --ref migrate/dorado-2.0.1-escpod --config config/config-test.yml \
#       --label new --profile cluster/slurm --setup
#   pixi run python benchmark/compare.py \
#       benchmark/snapshots/old/fingerprint benchmark/snapshots/new/fingerprint \
#       --profile aggregate

set -euo pipefail

REF="" CONFIG="" LABEL="" PROFILE="" CORES="" DO_SETUP=0
TARGETS="" KEEP_WORKTREE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref) REF="$2"; shift 2 ;;
        --config) CONFIG="$2"; shift 2 ;;
        --label) LABEL="$2"; shift 2 ;;
        --profile) PROFILE="$2"; shift 2 ;;
        --cores) CORES="$2"; shift 2 ;;
        --setup) DO_SETUP=1; shift ;;
        --targets) TARGETS="$2"; shift 2 ;;
        --keep-output) shift ;;   # accepted for clarity; output is kept by default
        --keep-worktree) KEEP_WORKTREE=1; shift ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -n "$REF" && -n "$CONFIG" && -n "$LABEL" ]] || {
    echo "ERROR: --ref, --config and --label are required" >&2; exit 2; }

# Resolve paths relative to the CURRENT repo (where this script lives).
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAP_DIR="$REPO_ROOT/benchmark/snapshots/$LABEL"
WORKTREE="$REPO_ROOT/benchmark/.worktrees/$LABEL"
OUTPUT_DIR="$SNAP_DIR/output"
FINGERPRINT_DIR="$SNAP_DIR/fingerprint"

echo ">> ref=$REF  label=$LABEL"
echo ">> worktree=$WORKTREE"
echo ">> output=$OUTPUT_DIR"

mkdir -p "$REPO_ROOT/benchmark/.worktrees" "$SNAP_DIR"

# Fresh worktree at the requested ref (with submodules — leech lives in one).
if git -C "$REPO_ROOT" worktree list --porcelain | grep -q "worktree $WORKTREE"; then
    git -C "$REPO_ROOT" worktree remove --force "$WORKTREE"
fi
git -C "$REPO_ROOT" worktree add --force --detach "$WORKTREE" "$REF"
git -C "$WORKTREE" submodule update --init --recursive

cd "$WORKTREE"

if [[ "$DO_SETUP" == 1 ]]; then
    echo ">> pixi run setup (per-ref tools: dorado/escpod/leech)"
    pixi run setup
fi

# Assemble the snakemake invocation. output_directory is overridden so each ref
# writes to its own snapshot dir and the two runs never collide.
SMK=(snakemake --configfile "$CONFIG" --config "output_directory=$OUTPUT_DIR" --rerun-incomplete)
if [[ -n "$PROFILE" ]]; then
    SMK+=(--profile "$PROFILE")
else
    SMK+=(--cores "${CORES:-4}")
fi
[[ -n "$TARGETS" ]] && SMK+=($TARGETS)

echo ">> ${SMK[*]}"
pixi run "${SMK[@]}"

# Fingerprint with the CURRENT checkout's scripts (old ref lacks benchmark/).
echo ">> fingerprinting -> $FINGERPRINT_DIR"
pixi run python "$REPO_ROOT/benchmark/fingerprint.py" \
    "$OUTPUT_DIR" "$FINGERPRINT_DIR" --label "$LABEL" --git-ref "$REF"

if [[ "$KEEP_WORKTREE" == 0 ]]; then
    cd "$REPO_ROOT"
    git worktree remove --force "$WORKTREE" || true
fi

echo ">> done. fingerprint at: $FINGERPRINT_DIR"
