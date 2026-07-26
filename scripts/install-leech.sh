#!/usr/bin/env bash
# Install leech from the GitHub release wheels.
#
# leech lives in a private repo and is not published to PyPI, so wheels are
# pulled from the release with `gh`. leech-core is a separate Rust extension
# published as a per-interpreter wheel in the same release; it is optional
# (leech falls back to a numpy path) but enables `--backend rust`.
#
# Set LEECH_WHEEL_DIR to install from a local directory of wheels instead.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

LEECH_VERSION="${LEECH_VERSION:-$(awk '/^leech_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
LEECH_CORE_VERSION="${LEECH_CORE_VERSION:-$(awk '/^leech_core_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"

py_tag=$(python -c "import sys; print(f'cp{sys.version_info.major}{sys.version_info.minor}')")
arch=$(uname -m)

if [ -n "${LEECH_WHEEL_DIR:-}" ]; then
    dl_dir="${LEECH_WHEEL_DIR}"
    echo "Using local wheels from ${dl_dir}"
else
    if ! command -v gh >/dev/null 2>&1; then
        echo "Error: gh not found. leech wheels come from a private GitHub release." >&2
        echo "Install the GitHub CLI and run 'gh auth login', or set LEECH_WHEEL_DIR." >&2
        exit 1
    fi
    dl_dir="${REPO_ROOT}/resources/tools/leech/${LEECH_VERSION}"
    mkdir -p "${dl_dir}"
    echo "Downloading leech ${LEECH_VERSION} wheels (${py_tag}, ${arch})..."
    gh release download "${LEECH_VERSION}" \
        --repo rnabioco/leech \
        --dir "${dl_dir}" \
        --clobber \
        --pattern "leech-*-py3-none-any.whl" \
        --pattern "leech_core-${LEECH_CORE_VERSION}-${py_tag}-${py_tag}-*manylinux*_${arch}.whl"
fi

leech_whl=$(find "${dl_dir}" -name 'leech-*-py3-none-any.whl' | head -1)
core_whl=$(find "${dl_dir}" -name "leech_core-*-${py_tag}-*.whl" | head -1)

if [ -z "${leech_whl}" ]; then
    echo "Error: no leech wheel found in ${dl_dir}" >&2
    exit 1
fi

# --no-deps throughout: numpy/pysam/scikit-learn/rich-click come from conda and
# torch is installed against the CUDA index by scripts/setup-tools.sh. Only the
# leech deps conda does not provide are added explicitly.
if [ -n "${core_whl}" ]; then
    echo "Installing leech-core from ${core_whl##*/}..."
    uv pip install --no-deps "${core_whl}"
else
    echo "Warning: no leech_core wheel for ${py_tag}/${arch}; leech will use the" >&2
    echo "slower pure-python backend (--backend rust will be unavailable)." >&2
fi

echo "Installing leech from ${leech_whl##*/}..."
uv pip install --no-deps "${leech_whl}"
uv pip install polars

# Reconcile pyarrow: ensure pip hasn't overridden conda's version
CONDA_PYARROW=$(pixi list pyarrow 2>/dev/null | awk '/^pyarrow[[:space:]]/ {print $2}')
if [ -n "$CONDA_PYARROW" ]; then
    echo "Reconciling pyarrow to conda version ${CONDA_PYARROW}..."
    uv pip install --force-reinstall "pyarrow==${CONDA_PYARROW}"
fi

echo "Leech installed successfully"
