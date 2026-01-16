#!/bin/bash
# Activation script for aa-tRNA-seq pipeline environment
# Sets up PATH and library paths - does NOT install packages
# For initial setup, run: pixi run setup

# Get the directory where this script is located, then go to repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ============================================================================
# Environment Variables
# ============================================================================
DORADO_VERSION="${DORADO_VERSION:-0.9.1}"
DORADO_DIR="${REPO_ROOT}/resources/tools/dorado/${DORADO_VERSION}"

# ============================================================================
# Library and PATH Setup
# ============================================================================
# Use conda/pixi libstdc++ instead of system version (fixes GLIBCXX version errors)
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"
export PATH="${DORADO_DIR}/bin:${PATH}"
