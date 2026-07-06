#!/usr/bin/env bash
# Install leech from submodule with release-optimized Rust build
set -euo pipefail

git submodule update --init --recursive resources/leech
# The superproject `--recursive` skips descending into leech when it is already
# checked out, leaving its nested escapepod-rs submodule (source of the
# escapepod-signal crate that leech-core builds against) uninitialized. Init it
# from within leech to be certain.
git -C resources/leech submodule update --init --recursive

echo "Installing leech-core (Rust, release build)..."
uv pip install resources/leech/rust

echo "Installing leech (Python, editable)..."
uv pip install --no-deps -e resources/leech

# leech imports `escapepod` (pyo3 pod5 reader) unconditionally at predict time,
# but declares it only as an optional extra, so --no-deps above does not pull it.
# Build and install it from the nested escapepod-rs submodule.
echo "Installing escapepod (Rust/pyo3 pod5 bindings for leech)..."
uv pip install resources/leech/escapepod-rs/crates/escapepod-python

# Reconcile pyarrow: ensure pip hasn't overridden conda's version
CONDA_PYARROW=$(pixi list pyarrow 2>/dev/null | awk '/^pyarrow[[:space:]]/ {print $2}')
if [ -n "$CONDA_PYARROW" ]; then
    echo "Reconciling pyarrow to conda version ${CONDA_PYARROW}..."
    uv pip install --force-reinstall "pyarrow==${CONDA_PYARROW}"
fi

echo "Leech installed successfully"
