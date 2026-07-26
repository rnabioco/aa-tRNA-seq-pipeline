#!/bin/bash
# One-time setup script for aa-tRNA-seq pipeline
# Run this once before using the pipeline: pixi run setup
# Do NOT run this in parallel from multiple nodes

set -euo pipefail

# Get the directory where this script is located, then go to repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ============================================================================
# Configuration
# ============================================================================
DORADO_VERSION="${DORADO_VERSION:-$(awk '/^dorado_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
DORADO_MODEL="${DORADO_MODEL:-$(awk '/^dorado_model:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
ESCAPEPOD_VERSION="${ESCAPEPOD_VERSION:-$(awk '/^escapepod_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
LEECH_VERSION="${LEECH_VERSION:-$(awk '/^leech_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
LEECH_CORE_VERSION="${LEECH_CORE_VERSION:-$(awk '/^leech_core_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
CUDA_VERSION="${CUDA_VERSION:-cu124}"
DORADO_DIR="${REPO_ROOT}/resources/tools/dorado/${DORADO_VERSION}"
MODEL_DIR="${REPO_ROOT}/resources/models"
TOOLS_DIR="${REPO_ROOT}/resources/tools"

# ============================================================================
# Helper Functions
# ============================================================================
detect_platform() {
    local system machine
    system=$(uname -s | tr '[:upper:]' '[:lower:]')
    machine=$(uname -m | tr '[:upper:]' '[:lower:]')

    case "${machine}" in
        arm64|aarch64) arch="arm64" ;;
        x86_64|amd64|x64) arch="x64" ;;
        *) echo "Error: Unsupported architecture: ${machine}" >&2; return 1 ;;
    esac

    case "${system}" in
        linux) os_suffix="linux-${arch}"; file_ext="tar.gz" ;;
        darwin) os_suffix="osx-${arch}"; file_ext="zip" ;;
        *) echo "Error: Unsupported OS: ${system}" >&2; return 1 ;;
    esac

    echo "${os_suffix}|${file_ext}"
}

download_dorado() {
    local platform_info os_suffix file_ext url tmpfile
    platform_info=$(detect_platform) || return 1
    os_suffix="${platform_info%|*}"
    file_ext="${platform_info#*|}"
    url="https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VERSION}-${os_suffix}.${file_ext}"
    tmpfile="/tmp/dorado.${file_ext}"

    echo "Downloading dorado ${DORADO_VERSION} for ${os_suffix}..."
    mkdir -p "${DORADO_DIR}"

    if ! curl -L -o "${tmpfile}" "${url}"; then
        echo "Error: Failed to download dorado from ${url}" >&2
        return 1
    fi

    echo "Extracting dorado..."
    if [ "${file_ext}" = "tar.gz" ]; then
        tar -xzf "${tmpfile}" -C "${DORADO_DIR}" --strip-components=1
    elif [ "${file_ext}" = "zip" ]; then
        local tmpdir="${DORADO_DIR}_temp"
        unzip -o "${tmpfile}" -d "${tmpdir}"
        mv "${tmpdir}"/*/* "${DORADO_DIR}/"
        rm -rf "${tmpdir}"
    fi

    rm -f "${tmpfile}"
    chmod +x "${DORADO_DIR}/bin/dorado"
    echo "Dorado installed to ${DORADO_DIR}"
}

download_model() {
    local model_path="${MODEL_DIR}/${DORADO_MODEL}"
    echo "Downloading dorado model ${DORADO_MODEL}..."
    mkdir -p "${MODEL_DIR}"

    if ! "${DORADO_DIR}/bin/dorado" download --model "${DORADO_MODEL}" --models-directory "${MODEL_DIR}"; then
        echo "Error: Failed to download model ${DORADO_MODEL}" >&2
        return 1
    fi

    touch "${model_path}/.downloaded"
    echo "Model installed to ${model_path}"
}

# ============================================================================
# Dorado Setup
# ============================================================================
echo "=== Checking dorado ==="
if [ -x "${DORADO_DIR}/bin/dorado" ]; then
    echo "Dorado already installed at ${DORADO_DIR}"
else
    download_dorado
fi

echo "=== Checking dorado model ==="
if [ -d "${MODEL_DIR}/${DORADO_MODEL}" ]; then
    echo "Model already installed at ${MODEL_DIR}/${DORADO_MODEL}"
else
    download_model
fi

# ============================================================================
# Modification Models
# ============================================================================
download_mod_models() {
    local mod_bases="${MODIFIED_BASES:-m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG}"
    for mod in ${mod_bases}; do
        local mod_model="${DORADO_MODEL}_${mod}@v1"
        local mod_path="${MODEL_DIR}/${mod_model}"
        if [ -d "${mod_path}" ]; then
            echo "  Modification model ${mod_model} already exists"
            continue
        fi
        echo "  Downloading ${mod_model}..."
        if ! "${DORADO_DIR}/bin/dorado" download --model "${mod_model}" --models-directory "${MODEL_DIR}"; then
            echo "Error: Failed to download modification model ${mod_model}" >&2
            return 1
        fi
    done
}

echo "=== Checking modification models ==="
download_mod_models

# ============================================================================
# Remora Setup (via uv)
# ============================================================================
echo "=== Checking remora ==="
if python -c "import remora" 2>/dev/null; then
    echo "Remora already installed"
else
    echo "Installing PyTorch with CUDA support (${CUDA_VERSION})..."
    uv pip install torch --index-url "https://download.pytorch.org/whl/${CUDA_VERSION}"

    echo "Installing remora dependencies (excluding pyarrow/numpy to preserve conda versions)..."
    uv pip install plotnine statsmodels thop

    echo "Installing ont-remora from GitHub..."
    uv pip install --no-deps "git+https://github.com/nanoporetech/remora.git"

    echo "Remora installed successfully"
fi

# ============================================================================
# escapepod Setup — escpod CLI (Rust) + escapepod python bindings
#
# The pipeline uses escpod for all POD5 manipulation (merge, filter) and for
# barcode demultiplexing. The published release binaries are built with default
# features only, so `escpod demux` is absent from them — the CLI is therefore
# built from source. Requires a Rust toolchain (>=1.95).
# ============================================================================
ESCPOD_BIN="${TOOLS_DIR}/escapepod/${ESCAPEPOD_VERSION}/bin/escpod"

echo "=== Checking escpod CLI ==="
if [ -x "${ESCPOD_BIN}" ]; then
    echo "escpod already installed at ${ESCPOD_BIN}"
else
    bash "${SCRIPT_DIR}/install-escpod.sh"
fi

echo "=== Checking escapepod python bindings ==="
ESCAPEPOD_PY_VERSION="${ESCAPEPOD_VERSION#v}"
if python -c "import escapepod" 2>/dev/null; then
    echo "escapepod $(python -c 'import escapepod; print(escapepod.__version__)') already installed"
else
    echo "Installing escapepod==${ESCAPEPOD_PY_VERSION}..."
    uv pip install "escapepod==${ESCAPEPOD_PY_VERSION}"
    echo "escapepod installed successfully"
fi

# ============================================================================
# ONT pod5 (legacy) — required only by the optional remora signal-metrics QC
# path (workflow/scripts/extract_signal_metrics.py), which hands pod5 objects
# to remora's io API and so cannot use escapepod's reader. Skipped unless
# remora is present.
# ============================================================================
echo "=== Checking ONT pod5 (remora signal-metrics path) ==="
POD5_MIN_VERSION="0.3.36"
current_pod5=$(python -c "import pod5; print(pod5.__version__)" 2>/dev/null || echo "0.0.0")
if python -c "from packaging.version import Version; exit(0 if Version('${current_pod5}') >= Version('${POD5_MIN_VERSION}') else 1)" 2>/dev/null; then
    echo "pod5 ${current_pod5} already installed (>= ${POD5_MIN_VERSION})"
else
    echo "Installing pod5 >= ${POD5_MIN_VERSION}..."
    uv pip install --no-deps "pod5>=${POD5_MIN_VERSION}"
    # pod5 needs 'deprecated' but --no-deps skips it
    uv pip install deprecated
    echo "pod5 installed successfully"
fi

# ============================================================================
# escpod demux models
#
# Best-effort: these come from the private rnabioco/escapepod-models repo, so a
# fresh clone without access cannot fetch them. Only needed when demultiplexing
# with the default `escpod` backend, so a failure here is a warning.
# ============================================================================
if [ -f "${REPO_ROOT}/resources/models/demux/barcode_wdx4_rna004.gbm.json" ]; then
    echo "=== escpod demux models already installed ==="
elif ! bash "${SCRIPT_DIR}/install-demux-models.sh"; then
    echo "Warning: could not install escpod demux models." >&2
    echo "Only required for demultiplexing with warpdemux.backend=escpod." >&2
    echo "Re-run 'pixi run install-demux-models' once you have access to" >&2
    echo "rnabioco/escapepod-models, or set warpdemux.backend to 'warpdemux'." >&2
fi

# ============================================================================
# WarpDemuX Setup (via uv)
#
# Only needed for warpdemux.backend=warpdemux. The default escpod backend
# reimplements this in Rust and needs no python WarpDemuX install.
# ============================================================================
echo "=== Checking WarpDemuX ==="
if python -c "import warpdemux" 2>/dev/null; then
    echo "WarpDemuX already installed"
else
    echo "Cloning WarpDemuX..."
    [ -d resources/tools/WarpDemuX ] || git clone --recursive https://github.com/KleistLab/WarpDemuX.git resources/tools/WarpDemuX
    echo "Installing WarpDemuX..."
    uv pip install -e resources/tools/WarpDemuX
    echo "WarpDemuX installed successfully"
fi

# ============================================================================
# Leech Setup — from GitHub release wheels
#
# leech is not on PyPI (private repo), so the prebuilt wheels are pulled from
# the release. leech-core is a separate Rust extension published as a
# per-interpreter wheel in the same release; it is optional (leech falls back to
# a numpy path) but gives the accelerated `--backend rust` extraction.
# Requires `gh` auth with read access to rnabioco/leech.
# ============================================================================
install_leech() {
    local py_tag arch dl_dir leech_whl core_whl
    py_tag=$(python -c "import sys; print(f'cp{sys.version_info.major}{sys.version_info.minor}')")
    arch=$(uname -m)
    dl_dir="${TOOLS_DIR}/leech/${LEECH_VERSION}"

    if ! command -v gh >/dev/null 2>&1; then
        echo "Error: gh not found. leech lives in a private repo and its wheels" >&2
        echo "are fetched from the GitHub release. Install the GitHub CLI and run" >&2
        echo "'gh auth login', or set LEECH_WHEEL_DIR to a directory holding the" >&2
        echo "leech and leech_core wheels." >&2
        return 1
    fi

    mkdir -p "${dl_dir}"
    echo "Downloading leech ${LEECH_VERSION} wheels (${py_tag}, ${arch})..."
    gh release download "${LEECH_VERSION}" \
        --repo rnabioco/leech \
        --dir "${dl_dir}" \
        --clobber \
        --pattern "leech-*-py3-none-any.whl" \
        --pattern "leech_core-${LEECH_CORE_VERSION}-${py_tag}-${py_tag}-*manylinux*_${arch}.whl"

    leech_whl=$(find "${dl_dir}" -name 'leech-*-py3-none-any.whl' | head -1)
    core_whl=$(find "${dl_dir}" -name "leech_core-*-${py_tag}-*.whl" | head -1)

    if [ -z "${leech_whl}" ]; then
        echo "Error: no leech wheel found in ${dl_dir}" >&2
        return 1
    fi

    # --no-deps throughout: numpy/pysam/scikit-learn/rich-click come from conda
    # and torch was installed above against the CUDA index. Only leech deps that
    # conda does not provide are added explicitly.
    if [ -n "${core_whl}" ]; then
        echo "Installing leech-core from ${core_whl##*/}..."
        uv pip install --no-deps "${core_whl}"
    else
        echo "Warning: no leech_core wheel for ${py_tag}/${arch}; leech will use" >&2
        echo "the slower pure-python backend (--backend rust will be unavailable)." >&2
    fi

    echo "Installing leech from ${leech_whl##*/}..."
    uv pip install --no-deps "${leech_whl}"
    uv pip install polars
    echo "Leech installed successfully"
    echo "NOTE: pyarrow will be reconciled with conda in the next step"
}

echo "=== Checking leech ==="
if python -c "import leech" 2>/dev/null; then
    echo "Leech already installed"
elif [ -n "${LEECH_WHEEL_DIR:-}" ]; then
    echo "Installing leech from ${LEECH_WHEEL_DIR}..."
    uv pip install --no-deps "${LEECH_WHEEL_DIR}"/leech_core-*.whl 2>/dev/null || true
    uv pip install --no-deps "${LEECH_WHEEL_DIR}"/leech-*-py3-none-any.whl
    uv pip install polars
else
    install_leech
fi

# ============================================================================
# Reconcile pyarrow: ensure pip hasn't overridden conda's version
# ============================================================================
echo "=== Reconciling pyarrow with conda ==="
CONDA_PYARROW_VERSION=$(pixi list pyarrow 2>/dev/null \
    | awk '/^pyarrow[[:space:]]/ {print $2}')
if [ -n "${CONDA_PYARROW_VERSION}" ]; then
    PIP_PYARROW_VERSION=$(python -c "import pyarrow; print(pyarrow.__version__)" 2>/dev/null || echo "")
    if [ "${PIP_PYARROW_VERSION}" != "${CONDA_PYARROW_VERSION}" ]; then
        echo "pyarrow mismatch: pip has ${PIP_PYARROW_VERSION}, conda expects ${CONDA_PYARROW_VERSION}"
        echo "Force-reinstalling pyarrow==${CONDA_PYARROW_VERSION}..."
        uv pip install --force-reinstall "pyarrow==${CONDA_PYARROW_VERSION}"
    else
        echo "pyarrow ${PIP_PYARROW_VERSION} matches conda — OK"
    fi
fi

echo "=== Setup complete ==="
