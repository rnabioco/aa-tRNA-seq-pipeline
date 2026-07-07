#!/usr/bin/env bash
#
# Build the "real" benchmark sample: a deterministic ~10k-read subset of a real
# E. coli aa-tRNA-seq run, frozen into a small pod5. Used by config-benchmark.yml
# to give the version-comparison rules statistically meaningful per-tRNA counts
# (the .tests/ data is only ~200 reads/sample).
#
# The pod5 output is NOT tracked in git (see .gitignore); it will live on S3.
# This script + the committed config make it reproducible: same source + N ->
# same reads (read-ids are sorted and the first N taken; UUIDs are random, so
# this is an unbiased subset with no RNG).
#
# Usage:
#   scripts/make_benchmark_sample.sh [SOURCE_POD5] [N_READS] [SAMPLE_DIR]
#
# Defaults target Batch1 Pool1 IVC (E. coli), 10k reads, benchmark/data/ecoli_ivc.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATA=/beevol/home/jhessel/devel/rnabioco/2026-aars-in-vitro/data
DEFAULT_SRC="$DATA/20260130_Batch1_IVC/Batch1_Pool1_IVC/20260130_1421_P2S-01617-B_PBG59663_3da4ea82/pod5/PBG59663_3da4ea82_629f05bb_5.pod5"

SRC="${1:-$DEFAULT_SRC}"
N="${2:-10000}"
SAMPLE_DIR="${3:-$REPO_ROOT/benchmark/data/ecoli_ivc}"

POD5_DIR="$SAMPLE_DIR/pod5"
IDS="$SAMPLE_DIR/read_ids.txt"
OUT="$POD5_DIR/ecoli_ivc.pod5"

[[ -f "$SRC" ]] || { echo "source pod5 not found: $SRC" >&2; exit 1; }
command -v escpod >/dev/null || { echo "escpod not on PATH (source scripts/setup-env.sh)" >&2; exit 1; }

mkdir -p "$POD5_DIR"

echo ">> source : $SRC"
echo ">> N reads: $N"
echo ">> sample : $SAMPLE_DIR"

# Deterministic read-id subset: list all ids, sort (C locale), take first N.
# Sort to a temp file rather than piping into `head` — `head` closing the pipe
# early would SIGPIPE `sort` and trip `set -o pipefail`.
echo ">> selecting $N read ids..."
allids="$SAMPLE_DIR/.allids.tmp"
pod5 view "$SRC" --ids 2>/dev/null | grep -v '^read_id' | LC_ALL=C sort >"$allids"
head -n "$N" "$allids" >"$IDS"
rm -f "$allids"
got=$(wc -l <"$IDS")
[[ "$got" -eq "$N" ]] || echo ">> WARNING: only $got ids available (< $N)"

echo ">> filtering to frozen pod5..."
escpod filter "$SRC" -i "$IDS" -o "$OUT" -f

# provenance
cat >"$SAMPLE_DIR/PROVENANCE.txt" <<EOF
source_pod5: $SRC
n_reads_requested: $N
n_reads_selected: $got
selection: read-ids sorted (C locale), first N
built_by: scripts/make_benchmark_sample.sh
EOF

echo ">> done: $OUT ($(du -h "$OUT" | cut -f1), $got reads)"
