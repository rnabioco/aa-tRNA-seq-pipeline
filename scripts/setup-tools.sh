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
ESCPOD_VERSION="${ESCPOD_VERSION:-$(awk '/^escpod_version:/ {print $2}' "${REPO_ROOT}/config/config-base.yml")}"
CUDA_VERSION="${CUDA_VERSION:-cu124}"
DORADO_DIR="${REPO_ROOT}/resources/tools/dorado/${DORADO_VERSION}"
ESCPOD_DIR="${REPO_ROOT}/resources/tools/escpod/${ESCPOD_VERSION}"
MODEL_DIR="${REPO_ROOT}/resources/models"

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

escpod_target_triple() {
    # Map uname to the Rust target triple used in escapepod-rs release assets.
    local system machine
    system=$(uname -s | tr '[:upper:]' '[:lower:]')
    machine=$(uname -m | tr '[:upper:]' '[:lower:]')

    case "${system}" in
        linux)
            case "${machine}" in
                x86_64|amd64|x64) echo "x86_64-unknown-linux-musl" ;;
                arm64|aarch64)    echo "aarch64-unknown-linux-musl" ;;
                *) echo "Error: Unsupported architecture: ${machine}" >&2; return 1 ;;
            esac ;;
        darwin)
            case "${machine}" in
                x86_64|amd64|x64) echo "x86_64-apple-darwin" ;;
                arm64|aarch64)    echo "aarch64-apple-darwin" ;;
                *) echo "Error: Unsupported architecture: ${machine}" >&2; return 1 ;;
            esac ;;
        *) echo "Error: Unsupported OS: ${system}" >&2; return 1 ;;
    esac
}

download_escpod() {
    # escpod (rnabioco/escapepod-rs) — prebuilt CLI binary from GitHub releases.
    # The repo is private, so fetch release assets with authenticated `gh`
    # (falls back to curl with GITHUB_TOKEN, which also works once public).
    local triple tarball tmpdir tmpfile sumsfile
    triple=$(escpod_target_triple) || return 1
    tarball="escpod-v${ESCPOD_VERSION}-${triple}.tar.gz"
    tmpdir="$(mktemp -d)"
    tmpfile="${tmpdir}/${tarball}"
    sumsfile="${tmpdir}/SHA256SUMS.txt"

    echo "Downloading escpod ${ESCPOD_VERSION} for ${triple}..."
    mkdir -p "${ESCPOD_DIR}"

    if command -v gh >/dev/null 2>&1; then
        if ! gh release download "v${ESCPOD_VERSION}" \
            --repo rnabioco/escapepod-rs \
            --pattern "${tarball}" --pattern "SHA256SUMS.txt" \
            --dir "${tmpdir}" --clobber; then
            echo "Error: 'gh release download' failed for escpod v${ESCPOD_VERSION}." >&2
            echo "       Ensure 'gh auth login' has access to rnabioco/escapepod-rs." >&2
            rm -rf "${tmpdir}"
            return 1
        fi
    else
        local base="https://github.com/rnabioco/escapepod-rs/releases/download/v${ESCPOD_VERSION}"
        local auth=()
        [ -n "${GITHUB_TOKEN:-}" ] && auth=(-H "Authorization: token ${GITHUB_TOKEN}")
        if ! curl -fL "${auth[@]}" -o "${tmpfile}" "${base}/${tarball}"; then
            echo "Error: Failed to download escpod (repo is private)." >&2
            echo "       Install 'gh' and run 'gh auth login', or set GITHUB_TOKEN." >&2
            rm -rf "${tmpdir}"
            return 1
        fi
        curl -fsL "${auth[@]}" -o "${sumsfile}" "${base}/SHA256SUMS.txt" || true
    fi

    # Best-effort checksum verification against the release SHA256SUMS.txt
    if [ -f "${sumsfile}" ]; then
        local expected actual
        expected=$(awk -v f="${tarball}" '$2 == f || $2 == "*"f {print $1}' "${sumsfile}" | head -1)
        if [ -n "${expected}" ]; then
            actual=$(sha256sum "${tmpfile}" | awk '{print $1}')
            if [ "${expected}" != "${actual}" ]; then
                echo "Error: escpod checksum mismatch (expected ${expected}, got ${actual})" >&2
                rm -rf "${tmpdir}"
                return 1
            fi
            echo "escpod checksum verified"
        fi
    fi

    echo "Extracting escpod..."
    tar -xzf "${tmpfile}" -C "${tmpdir}"
    # The binary may sit at the archive root or inside a versioned subdir; locate it.
    local bin
    bin=$(find "${tmpdir}" -type f -name escpod -perm -u+x 2>/dev/null | head -1)
    [ -n "${bin}" ] || bin=$(find "${tmpdir}" -type f -name escpod 2>/dev/null | head -1)
    if [ -z "${bin}" ]; then
        echo "Error: escpod binary not found in ${tarball}" >&2
        rm -rf "${tmpdir}"
        return 1
    fi
    cp "${bin}" "${ESCPOD_DIR}/escpod"
    chmod +x "${ESCPOD_DIR}/escpod"
    rm -rf "${tmpdir}"
    echo "escpod installed to ${ESCPOD_DIR}/escpod"
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

echo "=== Checking escpod ==="
if [ -x "${ESCPOD_DIR}/escpod" ]; then
    echo "escpod already installed at ${ESCPOD_DIR}"
else
    download_escpod
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
# Leech Setup (via uv, from submodule)
# ============================================================================
echo "=== Checking leech ==="
if python -c "import leech" 2>/dev/null; then
    echo "Leech already installed"
else
    # Ensure the submodule is checked out. The directory exists as a mount
    # point even when uninitialized, so test for actual contents and init if
    # needed (requires access to the private rnabioco/leech repo).
    if [ ! -f "${REPO_ROOT}/resources/leech/rust/Cargo.toml" ]; then
        echo "Initializing leech submodule..."
        git -C "${REPO_ROOT}" submodule update --init --recursive resources/leech
    fi
    # The nested escapepod-rs submodule provides the escapepod-signal crate that
    # leech-core builds against. The superproject `--recursive` above does NOT
    # descend into leech when it is already checked out, so init it explicitly
    # from within leech (requires access to the private rnabioco/escapepod-rs).
    if [ -d "${REPO_ROOT}/resources/leech/.git" ] || [ -f "${REPO_ROOT}/resources/leech/.git" ]; then
        if [ ! -f "${REPO_ROOT}/resources/leech/escapepod-rs/crates/escapepod-signal/Cargo.toml" ]; then
            echo "Initializing nested leech submodules (escapepod-rs)..."
            git -C "${REPO_ROOT}/resources/leech" submodule update --init --recursive
        fi
    fi
    if [ -f "${REPO_ROOT}/resources/leech/rust/Cargo.toml" ]; then
        echo "Installing leech-core (Rust, release build)..."
        uv pip install "${REPO_ROOT}/resources/leech/rust"
        echo "Installing leech (Python, editable)..."
        uv pip install --no-deps -e "${REPO_ROOT}/resources/leech"
        # leech imports `escapepod` (pyo3 pod5 reader) unconditionally at predict
        # time but declares it only as an optional extra, so --no-deps above does
        # not pull it. Build it from the nested escapepod-rs submodule.
        echo "Installing escapepod (Rust/pyo3 pod5 bindings for leech)..."
        uv pip install "${REPO_ROOT}/resources/leech/escapepod-rs/crates/escapepod-python"
        echo "Leech installed successfully"
        echo "NOTE: pyarrow will be reconciled with conda in the next step"
    else
        echo "Leech submodule not found at resources/leech"
        echo "Run 'git submodule update --init --recursive resources/leech' to clone it"
        echo "(requires access to the private rnabioco/leech repository)"
    fi
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
