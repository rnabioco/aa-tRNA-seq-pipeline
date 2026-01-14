#!/bin/bash
# Activation script for aa-tRNA-seq pipeline environment
# Downloads dorado and installs remora on first activation

# Get the directory where this script is located, then go to repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ============================================================================
# Dorado Setup
# ============================================================================
DORADO_VERSION="${DORADO_VERSION:-0.9.1}"
DORADO_MODEL="${DORADO_MODEL:-rna004_130bps_sup@v5.1.0}"
DORADO_DIR="${REPO_ROOT}/resources/tools/dorado/${DORADO_VERSION}"
MODEL_DIR="${REPO_ROOT}/resources/models"

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

# Download dorado if not present
if [ ! -x "${DORADO_DIR}/bin/dorado" ]; then
    download_dorado || echo "Warning: Failed to download dorado"
fi

# Download model if not present
if [ -x "${DORADO_DIR}/bin/dorado" ] && [ ! -d "${MODEL_DIR}/${DORADO_MODEL}" ]; then
    download_model || echo "Warning: Failed to download dorado model"
fi

# ============================================================================
# Remora Setup (via uv)
# ============================================================================
# Check if remora is installed, install if not
if ! python -c "import remora" 2>/dev/null; then
    echo "Installing ont-remora from GitHub..."
    uv pip install git+https://github.com/nanoporetech/remora.git 2>/dev/null || \
        echo "Warning: Failed to install remora (may need to run manually)"
fi

# ============================================================================
# PATH Setup
# ============================================================================
export PATH="${DORADO_DIR}/bin:${PATH}"
