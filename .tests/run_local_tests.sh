#!/bin/bash
# Local test script for aa-tRNA-seq-pipeline
# This script runs the same checks as the CI pipeline locally

set -e

echo "========================================"
echo "Running local tests for aa-tRNA-seq-pipeline"
echo "========================================"

# Check if conda/mamba environment is activated
if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "Error: No conda environment is activated."
    echo "Please activate the aatrnaseqpipe environment first:"
    echo "  mamba activate aatrnaseqpipe"
    exit 1
fi

echo ""
echo "[1/5] Checking Snakemake installation..."
if ! command -v snakemake &> /dev/null; then
    echo "Error: Snakemake is not installed in the current environment"
    exit 1
fi
echo "✓ Snakemake version: $(snakemake --version)"

echo ""
echo "[2/5] Validating configuration files..."
python3 -c "import yaml; yaml.safe_load(open('config/config-base.yml'))" || exit 1
python3 -c "import yaml; yaml.safe_load(open('config/config-test.yml'))" || exit 1
python3 -c "import yaml; yaml.safe_load(open('config/config-preprint.yml'))" || exit 1
echo "✓ All config files are valid YAML"

echo ""
echo "[3/5] Checking sample files..."
if [ ! -f "config/samples-test.tsv" ] || [ ! -s "config/samples-test.tsv" ]; then
    echo "Error: samples-test.tsv is missing or empty"
    exit 1
fi
echo "✓ Sample files validated"

echo ""
echo "[4/5] Running Snakemake dry-run (syntax check)..."
snakemake -n --configfile=config/config-test.yml || exit 1
echo "✓ Snakemake dry-run successful"

echo ""
echo "[5/5] Checking required directories..."
for dir in workflow workflow/rules workflow/scripts config resources cluster; do
    if [ ! -d "$dir" ]; then
        echo "Error: Required directory $dir not found"
        exit 1
    fi
done
echo "✓ All required directories exist"

echo ""
echo "========================================"
echo "All local tests passed! ✓"
echo "========================================"
echo ""
echo "To run the full pipeline test:"
echo "  1. Download test data: bash .tests/dl_test_data.sh"
echo "  2. Setup tools: snakemake setup_dorado dorado_model setup_modkit"
echo "  3. Run test: snakemake --cores 2 --configfile=config/config-test.yml"
