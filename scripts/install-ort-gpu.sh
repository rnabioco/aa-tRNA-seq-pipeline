#!/usr/bin/env bash
# Download the CUDA-enabled ONNX Runtime that `escpod demux --gpu` needs.
#
# Why this is not a conda/pixi dependency:
#
#   `ort 2.0.0-rc.13` (pinned by escapepod-rs) enables the `api-27` feature by
#   default, so ORT_API_VERSION is 27 and it refuses to dlopen any onnxruntime
#   older than 1.27.x — the failure is a hard `BadVersion` panic, not a
#   fallback. conda-forge's newest CUDA build is 1.26.0, so no conda package
#   satisfies it today. Hence a pinned tarball from the upstream release.
#
# The `gpu_cuda13` variant matches the CUDA 13 on this cluster's gpu nodes and
# bundles its own CUDA dependencies, so no system CUDA install is relied on.
# Switch to gpu_cuda12 if your nodes run CUDA 12.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ORT_VERSION="${ORT_VERSION:-1.27.1}"
ORT_VARIANT="${ORT_VARIANT:-gpu_cuda13}"
DEST="${REPO_ROOT}/resources/tools/onnxruntime"
NAME="onnxruntime-linux-x64-${ORT_VARIANT}-${ORT_VERSION}"
URL="https://github.com/microsoft/onnxruntime/releases/download/v${ORT_VERSION}/${NAME}.tgz"

if [ -f "${DEST}/${NAME}/lib/libonnxruntime.so.${ORT_VERSION}" ]; then
    echo "ONNX Runtime ${ORT_VERSION} (${ORT_VARIANT}) already installed at ${DEST}/${NAME}"
    exit 0
fi

mkdir -p "${DEST}"
tmp="$(mktemp -t ort.XXXXXX.tgz)"
trap 'rm -f "${tmp}"' EXIT

echo "Downloading ${NAME}..."
if ! curl -fL -o "${tmp}" "${URL}"; then
    echo "Error: failed to download ${URL}" >&2
    exit 1
fi

tar -xzf "${tmp}" -C "${DEST}"

lib="${DEST}/${NAME}/lib/libonnxruntime.so.${ORT_VERSION}"
if [ ! -f "${lib}" ]; then
    echo "Error: expected ${lib} after extraction" >&2
    exit 1
fi
if [ ! -f "${DEST}/${NAME}/lib/libonnxruntime_providers_cuda.so" ]; then
    echo "Error: CUDA execution provider missing — wrong variant downloaded?" >&2
    exit 1
fi

echo "Installed ${lib}"
echo
echo "The escapepod_demux rule sets ORT_DYLIB_PATH and LD_LIBRARY_PATH itself"
echo "when ldx.gpu is true, so nothing needs to be exported by hand."
