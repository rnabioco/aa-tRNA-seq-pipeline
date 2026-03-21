#!/usr/bin/env bash
# Install leech from submodule with release-optimized Rust build
set -euo pipefail

git submodule update --init --recursive resources/leech

echo "Installing leech-core (Rust, release build)..."
uv pip install resources/leech/rust

echo "Installing leech (Python, editable)..."
uv pip install --no-deps -e resources/leech

# Reconcile pyarrow: ensure pip hasn't overridden conda's version
CONDA_PYARROW=$(pixi list pyarrow 2>/dev/null | awk '/^pyarrow[[:space:]]/ {print $2}')
if [ -n "$CONDA_PYARROW" ]; then
    echo "Reconciling pyarrow to conda version ${CONDA_PYARROW}..."
    uv pip install --force-reinstall "pyarrow==${CONDA_PYARROW}"
fi

echo "Leech installed successfully"
