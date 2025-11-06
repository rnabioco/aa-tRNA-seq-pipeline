# Testing and CI/CD Documentation

This document describes the testing infrastructure and continuous integration (CI/CD) setup for the aa-tRNA-seq-pipeline.

## Overview

The pipeline includes several automated checks to ensure code quality and functionality:

1. **CI Workflow** (`.github/workflows/ci.yml`) - Syntax validation and integration tests
2. **Lint Workflow** (`.github/workflows/lint.yml`) - Code quality and formatting checks
3. **Local Test Script** (`.tests/run_local_tests.sh`) - Run tests locally before pushing
4. **Pre-commit Hooks** (`.pre-commit-config.yaml`) - Automatic checks before git commits

## CI Workflow

The main CI workflow runs on every push and pull request. It includes three jobs:

### 1. Syntax Check
- Sets up Mambaforge and the conda environment
- Verifies Snakemake installation
- Runs `snakemake -n` (dry-run) to validate workflow syntax

### 2. Pipeline Integration Test
- Downloads test data
- Sets up dorado and modkit tools
- Downloads the dorado basecalling model
- Runs the `merge_pods` rule as a basic integration test
- Validates that output files are created

### 3. Configuration Validation
- Validates YAML syntax of all config files
- Checks that sample files exist and have content
- Verifies required directory structure

## Lint Workflow

The linting workflow checks code quality and formatting:

### 1. Snakemake Linting
- Uses `snakefmt` to check Snakemake file formatting

### 2. Python Linting
- **black**: Checks Python code formatting
- **flake8**: Checks for Python code quality issues
- Both checks are non-blocking (won't fail the build)

### 3. Markdown Linting
- Uses `markdownlint` to check Markdown file formatting
- Non-blocking

### 4. YAML Linting
- Uses `yamllint` to check YAML file formatting
- Non-blocking

## Running Tests Locally

### Quick Syntax Check

Run the local test script to verify basic functionality:

```bash
# Activate the conda environment first
mamba activate aatrnaseqpipe

# Run the test script
bash .tests/run_local_tests.sh
```

This script checks:
- Snakemake installation
- Config file validity
- Sample file existence
- Workflow syntax (dry-run)
- Required directory structure

### Full Integration Test

To run a complete test of the pipeline:

```bash
# 1. Download test data (first time only)
bash .tests/dl_test_data.sh

# 2. Set up tools (first time only)
snakemake setup_dorado dorado_model setup_modkit --configfile=config/config-test.yml

# 3. Run the pipeline with test data
snakemake --cores 2 --configfile=config/config-test.yml
```

### Running Specific Tests

```bash
# Syntax check only
snakemake -n --configfile=config/config-test.yml

# Run a specific rule
snakemake merge_pods --configfile=config/config-test.yml --cores 1

# Force rerun of a specific rule
snakemake <rule_name> --forcerun <rule_name> --configfile=config/config-test.yml
```

## Pre-commit Hooks

Pre-commit hooks automatically check your code before each commit.

### Installation

```bash
# Install pre-commit
pip install pre-commit

# Install the git hooks
pre-commit install
```

### Usage

Once installed, the hooks will run automatically on `git commit`. To run manually:

```bash
# Run on all files
pre-commit run --all-files

# Run on staged files only
pre-commit run
```

### What Gets Checked

- Trailing whitespace removal
- End-of-file fixer
- YAML syntax validation
- Large file detection (max 1MB)
- Merge conflict detection
- Line ending normalization
- Python code formatting (black)
- Python linting (flake8)
- Snakemake formatting (snakefmt)

## CI Status Badges

The README includes badges showing the status of CI workflows:

- [![CI](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/ci.yml) - Main CI tests
- [![Lint](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/lint.yml/badge.svg)](https://github.com/rnabioco/aa-tRNA-seq-pipeline/actions/workflows/lint.yml) - Code quality checks

## Troubleshooting

### CI Failures

If the CI workflow fails:

1. Check the GitHub Actions logs for detailed error messages
2. Run the same commands locally to reproduce the issue
3. Use the local test script to verify basic functionality
4. Ensure all config files are valid YAML
5. Verify that test data paths are correct in `config/samples-test.tsv`

### Common Issues

**Conda environment not activated:**
```bash
Error: No conda environment is activated.
Solution: mamba activate aatrnaseqpipe
```

**Missing test data:**
```bash
Error: Test data not found
Solution: bash .tests/dl_test_data.sh
```

**Snakemake syntax errors:**
```bash
Solution: Check the error message and fix the syntax in the indicated file
Run: snakemake -n --configfile=config/config-test.yml
```

**Pre-commit hook failures:**
```bash
Solution: Fix the issues reported by the hooks
Run: pre-commit run --all-files
Bypass (not recommended): git commit --no-verify
```

## GPU Tests

Note: The CI pipeline does not run GPU-intensive rules (`rebasecall`, `classify_charging`) as GitHub Actions runners don't have GPU access. These rules should be tested locally on a system with GPU access or on the HPC cluster.

To test GPU rules locally:

```bash
# Ensure CUDA is available
nvidia-smi

# Run specific GPU rule
snakemake rebasecall --configfile=config/config-test.yml --cores 1
snakemake classify_charging --configfile=config/config-test.yml --cores 1
```

## Contributing

When contributing to this repository:

1. Install pre-commit hooks: `pre-commit install`
2. Run local tests before pushing: `bash .tests/run_local_tests.sh`
3. Ensure all CI checks pass on your pull request
4. Fix any linting issues reported by the Lint workflow

## Future Enhancements

Potential improvements to the testing infrastructure:

- [ ] Add unit tests for Python scripts in `workflow/scripts/`
- [ ] Add integration tests for individual Snakemake rules
- [ ] Set up test data generation/validation
- [ ] Add code coverage reporting
- [ ] Add performance benchmarking
- [ ] Create Docker container for reproducible testing
- [ ] Add GPU-enabled CI runners for full pipeline testing
