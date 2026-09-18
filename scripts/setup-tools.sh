#!/bin/bash
# One-time setup script for aa-tRNA-seq pipeline
# Run this once before using the pipeline: pixi run setup
# Do NOT run this in parallel from multiple nodes

set -euo pipefail

# Get the directory where this script is located, then go to repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." >/dev/null && pwd)"

# ============================================================================
# Configuration
# ============================================================================
DORADO_VERSION="${DORADO_VERSION:-$(awk '/^dorado_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
DORADO_MODEL="${DORADO_MODEL:-$(awk '/^dorado_model:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
DORADO_DIR="${REPO_ROOT}/resources/tools/dorado/${DORADO_VERSION}"
MODEL_DIR="${REPO_ROOT}/resources/models"

ESCPOD_VERSION="${ESCPOD_VERSION:-$(awk '/^escpod_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
ESCPOD_DIR="${REPO_ROOT}/resources/tools/escpod/${ESCPOD_VERSION}"
# Pinned checksums for the release tarballs, from the release's SHA256SUMS.txt.
# Pinned rather than fetched alongside the tarball so that re-tagging the
# release upstream is caught here instead of being silently trusted.
ESCPOD_SHA256_x86_64_linux="0550461dd2e61476c80a39823e38fe5b70844a89f663bcc467aaecdce50c05db"
ESCPOD_SHA256_aarch64_linux="cc3211811addcbb3f5ebac5069f95a293aeecb18fc01b7280a8f30671096b0b1"
ESCPOD_SHA256_x86_64_darwin="d8a43b3bdaa813c4ac615e3aa17661a38acfa0756d0293472da4901b6983321f"
ESCPOD_SHA256_aarch64_darwin="7ee3583a2050337b8aafdbbf2590f11557c5a18cea75f9bf4fac4ba694054075"
# The GPU build, for `ldx.gpu: true`. It is a SEPARATE artifact and the only
# dynamically linked one (glibc >= 2.28), because the CUDA runtimes are
# dlopened and so cannot be static-musl. x86_64 Linux only — upstream publishes
# no other GPU target.
ESCPOD_SHA256_x86_64_linux_gpu="e54dbb8ee7b8112e027d2e236e50d67b5e133207257aceffd4b0763af6ca374c"
ESCPOD_GPU_TARGET="x86_64-unknown-linux-gnu-gpu"
ESCPOD_GPU_DIR="${REPO_ROOT}/resources/tools/escpod/${ESCPOD_VERSION}-gpu"

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

escpod_target() {
    # Rust target triple for the escpod release artifacts. The Linux builds are
    # static musl, so they run regardless of the host glibc.
    local system machine
    system=$(uname -s | tr '[:upper:]' '[:lower:]')
    machine=$(uname -m | tr '[:upper:]' '[:lower:]')

    case "${machine}" in
        arm64|aarch64) machine="aarch64" ;;
        x86_64|amd64|x64) machine="x86_64" ;;
        *) echo "Error: Unsupported architecture: ${machine}" >&2; return 1 ;;
    esac

    case "${system}" in
        linux) echo "${machine}-unknown-linux-musl" ;;
        darwin) echo "${machine}-apple-darwin" ;;
        *) echo "Error: Unsupported OS: ${system}" >&2; return 1 ;;
    esac
}

download_escpod() {
    local target url tmpfile expected actual
    target=$(escpod_target) || return 1

    case "${target}" in
        x86_64-unknown-linux-musl)  expected="${ESCPOD_SHA256_x86_64_linux}" ;;
        aarch64-unknown-linux-musl) expected="${ESCPOD_SHA256_aarch64_linux}" ;;
        x86_64-apple-darwin)        expected="${ESCPOD_SHA256_x86_64_darwin}" ;;
        aarch64-apple-darwin)       expected="${ESCPOD_SHA256_aarch64_darwin}" ;;
        *) echo "Error: no pinned checksum for target ${target}" >&2; return 1 ;;
    esac

    url="https://github.com/rnabioco/escapepod-rs/releases/download/v${ESCPOD_VERSION}/escpod-v${ESCPOD_VERSION}-${target}.tar.gz"
    tmpfile="$(mktemp -t escpod.XXXXXX.tar.gz)"

    echo "Downloading escpod ${ESCPOD_VERSION} for ${target}..."
    if ! curl -fL -o "${tmpfile}" "${url}"; then
        echo "Error: Failed to download escpod from ${url}" >&2
        rm -f "${tmpfile}"
        return 1
    fi

    actual=$(sha256sum "${tmpfile}" | awk '{print $1}')
    if [ "${actual}" != "${expected}" ]; then
        echo "Error: escpod checksum mismatch for ${target}" >&2
        echo "  expected ${expected}" >&2
        echo "  actual   ${actual}" >&2
        rm -f "${tmpfile}"
        return 1
    fi

    mkdir -p "${ESCPOD_DIR}/bin"
    # The tarball holds the bare `escpod` binary, so extract straight into bin/.
    tar -xzf "${tmpfile}" -C "${ESCPOD_DIR}/bin"
    rm -f "${tmpfile}"
    chmod +x "${ESCPOD_DIR}/bin/escpod"
    echo "escpod installed to ${ESCPOD_DIR}"
}

download_escpod_gpu() {
    # The `ldx.gpu: true` binary. Kept separate from download_escpod because it
    # is a different target triple, a different linkage, and optional: a run
    # with `ldx.gpu: false`, or on any host that is not x86_64 Linux, never
    # touches it. Skipping is therefore not an error here — get_escpod_bin in
    # demux.smk is what fails, and only if a run actually asks for the GPU.
    local url tmpfile actual
    if [ "$(uname -s)" != "Linux" ] || [ "$(uname -m)" != "x86_64" ]; then
        echo "Skipping escpod GPU build: published for x86_64 Linux only"
        return 0
    fi

    url="https://github.com/rnabioco/escapepod-rs/releases/download/v${ESCPOD_VERSION}/escpod-v${ESCPOD_VERSION}-${ESCPOD_GPU_TARGET}.tar.gz"
    tmpfile="$(mktemp -t escpod-gpu.XXXXXX.tar.gz)"

    echo "Downloading escpod ${ESCPOD_VERSION} (GPU) for ${ESCPOD_GPU_TARGET}..."
    if ! curl -fL -o "${tmpfile}" "${url}"; then
        echo "Error: Failed to download escpod GPU build from ${url}" >&2
        rm -f "${tmpfile}"
        return 1
    fi

    actual=$(sha256sum "${tmpfile}" | awk '{print $1}')
    if [ "${actual}" != "${ESCPOD_SHA256_x86_64_linux_gpu}" ]; then
        echo "Error: escpod GPU checksum mismatch" >&2
        echo "  expected ${ESCPOD_SHA256_x86_64_linux_gpu}" >&2
        echo "  actual   ${actual}" >&2
        rm -f "${tmpfile}"
        return 1
    fi

    mkdir -p "${ESCPOD_GPU_DIR}/bin"
    tar -xzf "${tmpfile}" -C "${ESCPOD_GPU_DIR}/bin"
    rm -f "${tmpfile}"
    chmod +x "${ESCPOD_GPU_DIR}/bin/escpod"
    echo "escpod (GPU) installed to ${ESCPOD_GPU_DIR}"
    echo "  ldx.gpu: true additionally needs a CUDA libonnxruntime"
    echo "  (pixi run install-ort-gpu) and cuDNN (pixi install -e gpu)."
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
# Pod5 Setup (via uv — bioconda version is outdated)
# ============================================================================
echo "=== Checking pod5 ==="
POD5_MIN_VERSION="0.3.36"
current_pod5=$(python -c "import pod5; print(pod5.__version__)" 2>/dev/null || echo "0.0.0")
if python -c "from packaging.version import Version; exit(0 if Version('${current_pod5}') >= Version('${POD5_MIN_VERSION}') else 1)" 2>/dev/null; then
    echo "Pod5 ${current_pod5} already installed (>= ${POD5_MIN_VERSION})"
else
    echo "Installing pod5 >= ${POD5_MIN_VERSION}..."
    uv pip install --no-deps "pod5>=${POD5_MIN_VERSION}"
    # pod5 needs 'deprecated' but --no-deps skips it
    uv pip install deprecated
    echo "Pod5 installed successfully"
fi

# ============================================================================
# escapepod (escpod) Setup
# ============================================================================
# Provides `escpod classify` (the tRNA charging classifier), `escpod
# merge`/`escpod filter` (POD5 handling), and `escpod demux` (the CTC-CRF
# barcode demultiplexer used for LDX/nbc barcodes). The models themselves are
# vendored in resources/models/charging/ and resources/models/demux/ rather
# than fetched — see the READMEs there.
echo "=== Checking escpod ==="
if [ -x "${ESCPOD_DIR}/bin/escpod" ]; then
    echo "escpod already installed at ${ESCPOD_DIR}"
else
    download_escpod
fi

# The GPU build is fetched alongside it, so that `ldx.gpu: true` (the default)
# works straight out of setup. Since escapepod-rs 0.17.1 this is a download
# rather than a source build of a private repo.
if [ -x "${ESCPOD_GPU_DIR}/bin/escpod" ]; then
    echo "escpod (GPU) already installed at ${ESCPOD_GPU_DIR}"
else
    download_escpod_gpu
fi

# ============================================================================
# WarpDemuX Setup (via uv)
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
