#!/usr/bin/env bash
# Build and install the escpod CLI with the demux feature.
#
# The published escapepod-rs release binaries are built with default features
# only, so `escpod demux` is absent from them. Barcode demultiplexing therefore
# requires a source build, which needs a Rust toolchain (>=1.95).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ESCAPEPOD_VERSION="${ESCAPEPOD_VERSION:-$(awk '/^escapepod_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
PREFIX="${REPO_ROOT}/resources/tools/escapepod/${ESCAPEPOD_VERSION}"

# Cargo features. `cnn-detect` implies `demux` and additionally enables
# `--method cnn` adapter-boundary detection, which the shipped barcode GBM model
# was trained behind: the CNN and LLR boundaries agree only ~82% within +/-200
# samples, so LLR-only is a train/serve skew, not a free choice.
# `gpu` adds CUDA-accelerated DTW for demux classify (needs the CUDA driver and
# libnvrtc at runtime); `train` adds `escpod demux train-svm`.
ESCPOD_FEATURES="${ESCPOD_FEATURES:-cnn-detect}"

if ! command -v cargo >/dev/null 2>&1; then
    echo "Error: cargo not found. Install Rust >=1.95:" >&2
    echo "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh" >&2
    exit 1
fi

echo "Building escpod ${ESCAPEPOD_VERSION} (features: ${ESCPOD_FEATURES})..."
cargo install \
    --git https://github.com/rnabioco/escapepod-rs \
    --tag "${ESCAPEPOD_VERSION}" \
    --features "${ESCPOD_FEATURES// /,}" \
    --root "${PREFIX}" \
    --locked \
    escapepod-cli

# A default-features build still accepts `escpod demux` on the command line but
# errors at runtime from a stub, which would only surface mid-pipeline. Check now.
if "${PREFIX}/bin/escpod" demux --help >/dev/null 2>&1; then
    echo "escpod installed to ${PREFIX}/bin/escpod (demux available)"
else
    echo "Error: escpod built without a working demux subcommand." >&2
    echo "Check that 'demux' is in ESCPOD_FEATURES (got: ${ESCPOD_FEATURES})." >&2
    exit 1
fi
