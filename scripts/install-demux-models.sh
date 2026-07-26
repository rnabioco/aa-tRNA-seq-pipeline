#!/usr/bin/env bash
# Install the escpod barcode-demultiplexing models into resources/models/demux.
#
# Why these models and not the WarpDemuX ones
# -------------------------------------------
# `escpod demux` cannot load a WarpDemuX kit directly. The tRNA kits
# (WDX4_tRNA_rna004_v1_0, WDX4b_tRNA_rna004_v1_0) are `warpdemux.models.fpt_boost.Fpt_Boost`
# objects wrapping a *CatBoost* classifier, and every converter shipped with
# escapepod-rs requires either a scikit-learn SVC (convert_warpdemux_model.py
# needs `_X`) or a scikit-learn HistGradientBoostingClassifier
# (export_gbm_model.py needs `_predictors`). Neither matches CatBoost, so
# conversion is not possible.
#
# The supported replacement is `barcode_wdx4_rna004`: a GBM *distilled* from the
# WDX4_tRNA_rna004_v1_0 teacher's high-confidence calls (conf >= 0.9), trained
# and served on the same Rust `--warpdemux-compat` fingerprint. It reports 0.971
# balanced recall against the teacher, and its label_mapper ({0:3, 1:4, 2:5, 3:7})
# is identical to the tRNA kit's, so barcode numbering is unchanged.
#
# Sources, in the order this script tries them
# --------------------------------------------
#   1. ESCAPEPOD_MODELS_DIR — a local rnabioco/escapepod-models checkout.
#   2. A GitHub release on rnabioco/escapepod-models (requires authenticated gh;
#      the repo is private). NOTE: as of escapepod-models @ main the barcode
#      models are NOT released — only adapter_rna004@v1.0.1 is. Publish with
#      `scripts/release_model.sh barcode_wdx4_rna004 1.0.0` from that checkout
#      to make this path work.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEST="${REPO_ROOT}/resources/models/demux"

# Model id@version pairs to install, and the artifact filename within each.
BARCODE_MODEL="${BARCODE_MODEL:-barcode_wdx4_rna004@v1.0.0}"
BARCODE_FILE="${BARCODE_MODEL%@*}.gbm.json"
# Boundary detector used by `--method cnn`. v1.0.0 is deprecated (static
# sequence axis; returns adapter_end=0 for every read) — v1.0.1 or newer only.
ADAPTER_MODEL="${ADAPTER_MODEL:-adapter_rna004@v1.0.1}"
ADAPTER_FILE="${ADAPTER_MODEL%@*}.onnx"

mkdir -p "${DEST}"

# sha256 values from escapepod-models/models/MANIFEST.json. Verified after copy
# so a truncated download or a stale vendored file fails here rather than
# producing silently wrong barcode calls.
declare -A EXPECTED_SHA256=(
    ["barcode_wdx4_rna004.gbm.json"]="807f5fe2"
    ["barcode_5class_rna004.gbm.json"]="d002b412"
    ["adapter_rna004.onnx"]="24232a82"
)

verify_sha256() {
    local path="$1" name expected actual
    name="$(basename "${path}")"
    expected="${EXPECTED_SHA256[${name}]:-}"
    [ -z "${expected}" ] && return 0
    actual=$(sha256sum "${path}" | cut -c1-8)
    if [ "${actual}" != "${expected}" ]; then
        echo "Error: ${name} sha256 starts ${actual}, expected ${expected}" >&2
        echo "Refusing to install a model that does not match MANIFEST.json." >&2
        return 1
    fi
    echo "  verified ${name} (sha256 ${actual}...)"
}

install_from_checkout() {
    local src_root="$1" id file src
    for spec in "${BARCODE_MODEL}:${BARCODE_FILE}" "${ADAPTER_MODEL}:${ADAPTER_FILE}"; do
        id="${spec%%:*}"
        file="${spec##*:}"
        src="${src_root}/models/${id}/${file}"
        if [ ! -f "${src}" ]; then
            echo "Error: ${src} not found in checkout" >&2
            return 1
        fi
        cp "${src}" "${DEST}/${file}"
        # provenance.json records the teacher, training data, and metrics; keep
        # it beside the model so a run's barcode calls stay attributable.
        [ -f "${src_root}/models/${id}/provenance.json" ] &&
            cp "${src_root}/models/${id}/provenance.json" "${DEST}/${file%.*}.provenance.json"
        verify_sha256 "${DEST}/${file}"
    done
}

install_from_release() {
    local id
    if ! command -v gh >/dev/null 2>&1; then
        echo "Error: gh not found and ESCAPEPOD_MODELS_DIR not set." >&2
        return 1
    fi
    for id in "${BARCODE_MODEL}" "${ADAPTER_MODEL}"; do
        echo "Downloading ${id} from rnabioco/escapepod-models..."
        if ! gh release download "${id}" \
            --repo rnabioco/escapepod-models \
            --dir "${DEST}" \
            --clobber \
            --pattern "*.zip"; then
            echo "Error: no release found for ${id}." >&2
            echo "The barcode GBM models are not published yet. Either set" >&2
            echo "ESCAPEPOD_MODELS_DIR to a local escapepod-models checkout, or" >&2
            echo "publish the model from that checkout with:" >&2
            echo "  scripts/release_model.sh ${id%@*} ${id#*@v}" >&2
            return 1
        fi
        (cd "${DEST}" && unzip -o "${id}.zip" && rm -f "${id}.zip")
    done
    verify_sha256 "${DEST}/${BARCODE_FILE}"
    verify_sha256 "${DEST}/${ADAPTER_FILE}"
}

echo "=== Installing escpod demux models into ${DEST} ==="

if [ -n "${ESCAPEPOD_MODELS_DIR:-}" ]; then
    echo "Using local escapepod-models checkout: ${ESCAPEPOD_MODELS_DIR}"
    install_from_checkout "${ESCAPEPOD_MODELS_DIR}"
elif [ -d "${REPO_ROOT}/../escapepod-models/models" ]; then
    echo "Using sibling escapepod-models checkout"
    install_from_checkout "$(cd "${REPO_ROOT}/../escapepod-models" && pwd)"
else
    install_from_release
fi

echo "=== Demux models installed ==="
echo "  barcode:  ${DEST}/${BARCODE_FILE}"
echo "  boundary: ${DEST}/${ADAPTER_FILE}"
echo
echo "NOTE: the shipped barcode GBM carries no per-class confidence thresholds,"
echo "so escpod assigns every read with a usable adapter boundary to some"
echo "barcode. WarpDemuX instead rejected low-confidence reads to 'unclassified'."
echo "Expect a smaller unclassified fraction. Filter on the confidence column of"
echo "the classifications table if you need WarpDemuX-like rejection."
