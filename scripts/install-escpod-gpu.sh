#!/usr/bin/env bash
# Build escpod with the GPU features and install it as a distinct version.
#
# The published escpod release is built with the default `cli` feature set,
# which includes the CPU paths (cnn-detect, crf-decode via tract) but NOT
# `crf-gpu` / `cnn-gpu`. Those need onnxruntime's CUDA execution provider, which
# a portable static musl binary cannot assume, so the release binary has no GPU
# code at all and rejects `--gpu` outright. Running demux on a GPU therefore
# requires building from source.
#
# Installed under a `-gpu` version suffix rather than overwriting the release
# build, so both remain available and `escpod_version` in the config selects
# which one the pipeline uses (the Snakefile resolves
# resources/tools/escpod/<escpod_version>/bin).
#
# This is a large Rust build — run it inside a Slurm allocation, e.g.
#   srun -p rna -c 32 --mem=64G --time=1:00:00 pixi run install-escpod-gpu
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC="${REPO_ROOT}/resources/leech/escapepod-rs"

# The GPU features first exist in escapepod-rs 0.7.0; older checkouts build
# fine and then silently lack --gpu, which is worse than failing here.
ESCPOD_REF="${ESCPOD_REF:-v0.7.0}"

if [ ! -f "${SRC}/Cargo.toml" ]; then
    echo "escapepod-rs source not found at ${SRC}" >&2
    echo "Initialize it with:" >&2
    echo "  git -C ${REPO_ROOT} submodule update --init --recursive resources/leech" >&2
    exit 1
fi

if ! command -v cargo >/dev/null; then
    echo "Error: cargo not found. Install Rust to build escpod with GPU support." >&2
    exit 1
fi

current="$(git -C "${SRC}" describe --tags 2>/dev/null || echo unknown)"
if [ "${current}" != "${ESCPOD_REF}" ]; then
    echo "Checking out escapepod-rs ${ESCPOD_REF} (was ${current})..."
    git -C "${SRC}" fetch --tags --quiet origin
    git -C "${SRC}" checkout --quiet "${ESCPOD_REF}"
fi

if ! grep -qE '^crf-gpu =' "${SRC}/crates/escapepod-cli/Cargo.toml"; then
    echo "Error: ${ESCPOD_REF} has no crf-gpu feature — too old for GPU demux." >&2
    exit 1
fi

echo "Building escpod with crf-gpu,cnn-gpu (this takes ~10 min)..."
cargo build --release --manifest-path "${SRC}/Cargo.toml" \
    -p escapepod-cli --bin escpod --features crf-gpu,cnn-gpu

version="$("${SRC}/target/release/escpod" --version | awk '{print $2}')"
dest="${REPO_ROOT}/resources/tools/escpod/${version}-gpu/bin"
mkdir -p "${dest}"
install -m 0755 "${SRC}/target/release/escpod" "${dest}/escpod"

if ! "${dest}/escpod" demux --help 2>&1 | grep -q -- '--gpu'; then
    echo "Error: built binary has no --gpu flag; the features did not take." >&2
    exit 1
fi

echo "Installed ${dest}/escpod"
echo
echo "To use it, set in your run config:"
echo "  escpod_version: ${version}-gpu"
echo "  ldx: { gpu: true }"
echo "and run 'pixi run install-ort-gpu' for the CUDA onnxruntime it loads."
