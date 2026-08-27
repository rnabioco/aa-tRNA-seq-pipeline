#!/usr/bin/env bash
# Build escpod with the GPU features and install it as a distinct version.
#
# NOT NEEDED FOR A RELEASED VERSION since escapepod-rs 0.17.1, which publishes
# a GPU artifact (`escpod-v<ver>-x86_64-unknown-linux-gnu-gpu.tar.gz`) that
# `pixi run setup` downloads into the same `<version>-gpu` path this script
# writes. Use setup unless you need a ref that has no release — an unmerged
# branch, or a tag whose GPU artifact failed to build (as v0.17.0's did).
#
# The PORTABLE musl release is still built with the default `cli` feature set,
# which includes the CPU paths (cnn-detect, crf-decode via tract) but NOT the
# `gpu` feature. That needs onnxruntime's CUDA execution provider, which a
# static musl binary cannot dlopen — which is why the GPU artifact is a
# separate, dynamically linked, x86_64-Linux-only build rather than a flag on
# the usual one.
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
# escapepod-rs source. Cloned standalone under resources/, NOT taken from a
# submodule: it used to live nested inside resources/leech, which coupled the
# GPU build to an unrelated dependency and broke when that was dropped. Override
# with ESCPOD_SRC to build from a checkout you already have.
SRC="${ESCPOD_SRC:-${REPO_ROOT}/resources/escapepod-rs}"
ESCPOD_URL="${ESCPOD_URL:-https://github.com/rnabioco/escapepod-rs}"

# The GPU features first exist in escapepod-rs 0.7.0; older checkouts build
# fine and then silently lack any GPU path, which is worse than failing here.
#
# 0.8.0 is the floor for *full* GPU support: 0.7.0 runs the boundary CNN and the
# CRF encoder on the device but drops back to the CPU for the lattice decode,
# which then dominates. 0.8.0 adds the batched GPU lattice (`crf::lattice_gpu`),
# taking the CRF head end to end on the device (4.3x, escapepod-rs#186).
#
# The ref to build is DERIVED from `escpod_version` in config/config-base.yml,
# not hardcoded. The two have to match: the pipeline resolves the GPU binary as
# `<escpod_version>-gpu`, so a default that drifts from the config silently
# builds a version nothing will ever look for. Override ESCPOD_REF to build
# something else deliberately.
#
# Floors that still apply, in case someone points this backwards:
#   - 0.9.0: `demux --annotate`, which the LDX rule passes; absent before it.
#   - 0.8.1: the `boundary.margin` (#193) / `boundary.clamp_max_shift` (#194)
#     bundle keys the vendored model declares, plus the POD5 writer fix (#195).
#     An older escpod skips those keys silently and loses ~7% of a flowcell with
#     nothing to say so; one older than #195 writes POD5s dorado reads short.
ESCPOD_REF="${ESCPOD_REF:-v$(awk '/^escpod_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"

# Clone on demand rather than requiring a manual init step. Only the login node
# needs this; compute nodes build from the checkout it leaves behind.
if [ ! -f "${SRC}/Cargo.toml" ]; then
    if [ -n "${ESCPOD_SRC:-}" ]; then
        echo "Error: ESCPOD_SRC=${SRC} has no Cargo.toml." >&2
        exit 1
    fi
    echo "Cloning escapepod-rs into ${SRC}..."
    git clone --quiet "${ESCPOD_URL}" "${SRC}" || {
        echo "Error: could not clone ${ESCPOD_URL}." >&2
        echo "It is private; clone it yourself and point ESCPOD_SRC at it." >&2
        exit 1
    }
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

if ! grep -qE '^gpu =' "${SRC}/crates/escapepod-cli/Cargo.toml"; then
    echo "Error: ${ESCPOD_REF} has no gpu feature — too old for GPU demux." >&2
    exit 1
fi

# `gpu` is the meta-feature: demux + crf-gpu + cnn-gpu. The granular flags still
# exist, but naming them individually is how you end up with the encoder on the
# device and the detector left behind on the CPU.
echo "Building escpod with the gpu feature (this takes ~10 min)..."
cargo build --release --manifest-path "${SRC}/Cargo.toml" \
    -p escapepod-cli --bin escpod --features gpu

version="$("${SRC}/target/release/escpod" --version | awk '{print $2}')"
dest="${REPO_ROOT}/resources/tools/escpod/${version}-gpu/bin"
mkdir -p "${dest}"
install -m 0755 "${SRC}/target/release/escpod" "${dest}/escpod"

# Verified against the BINARY, not `--help`. Since 0.17.1 device placement is
# spelled `--device <auto|cpu|gpu>`, which clap renders identically whether or
# not the gpu features are compiled in — `escpod demux --help` is byte-for-byte
# the same for the musl release and the GPU build, so no amount of grepping the
# help text can tell them apart. The CUDA execution provider string is only
# linked in when the features actually took.
if ! grep -qa 'CUDAExecutionProvider' "${dest}/escpod"; then
    echo "Error: built binary has no CUDA execution provider;" >&2
    echo "       the gpu features did not take." >&2
    exit 1
fi

echo "Installed ${dest}/escpod"
echo
echo "To use it, set in your run config:"
echo "  escpod_version: ${version}-gpu"
echo "  ldx: { gpu: true }"
echo "and run 'pixi run install-ort-gpu' for the CUDA onnxruntime it loads."
