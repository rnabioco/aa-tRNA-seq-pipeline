# Complete Test and Build Check Setup

This document describes the comprehensive test and build infrastructure for the aa-tRNA-seq-pipeline.

## ✅ What's Installed

### 1. GitHub Actions CI/CD

**Location**: `.github/workflows/`

#### CI Pipeline (`ci.yml`)
Runs on: push to main/master/develop/claude/**, pull requests, manual trigger

- **Syntax Check Job**
  - Sets up conda environment with Snakemake
  - Validates Snakemake syntax with dry-run

- **Pipeline Integration Test Job**
  - Downloads test data
  - Sets up dorado and modkit tools
  - Runs merge_pods rule (non-GPU test)
  - Validates output directory creation

- **Configuration Validation Job**
  - Validates YAML syntax in all config files
  - Checks samples.tsv files exist and have content
  - Verifies required directory structure

#### Lint Pipeline (`lint.yml`)
Runs on: push to main/master/develop/claude/**, pull requests, manual trigger

- **Snakemake Linting**: snakefmt format checking
- **Python Linting**: black and flake8 on workflow/scripts/
- **Markdown Linting**: markdownlint on all .md files
- **YAML Linting**: yamllint on config/ and .github/

### 2. Pre-commit Hooks

**Location**: `.pre-commit-config.yaml`

Automatically runs before each commit:
- Trailing whitespace removal
- End-of-file fixing
- YAML syntax validation
- Large file detection (>1MB warning)
- Merge conflict detection
- Line ending normalization
- Black formatting for Python scripts
- Flake8 linting for Python scripts
- Snakefmt formatting for Snakemake files

**Installation**:
```bash
pip install pre-commit
pre-commit install
```

### 3. Local Test Script

**Location**: `.tests/run_local_tests.sh`

Runs the same checks as CI locally:
1. Verifies Snakemake installation
2. Validates configuration files (YAML syntax)
3. Checks sample files exist
4. Runs Snakemake dry-run (syntax check)
5. Verifies required directory structure

**Usage**:
```bash
# Requires active conda environment
mamba activate aatrnaseqpipe
bash .tests/run_local_tests.sh
```

### 4. Claude Code Session Hooks

**Location**: `.claude/hooks/SessionStart`

Automatically runs when starting a Claude Code session:
- Verifies conda/mamba availability
- Creates aatrnaseqpipe environment if missing
- Activates the environment
- Validates Snakemake installation
- Checks config file syntax
- Installs pre-commit hooks if not present
- Runs quick Snakemake syntax check
- Displays helpful development commands

## 📋 Development Workflow

### Starting a New Session

When you start working with Claude Code, the SessionStart hook automatically:
1. Sets up your environment
2. Validates configurations
3. Shows available commands

### Before Making Changes

```bash
# Run full local test suite
bash .tests/run_local_tests.sh

# Or just syntax check
snakemake -n --configfile=config/config-test.yml
```

### During Development

```bash
# Test specific rule
snakemake <rule_name> -n --configfile=config/config-test.yml

# Run linting manually
pre-commit run --all-files

# Format Snakemake files
snakefmt workflow/

# Format Python scripts
black workflow/scripts/
```

### Before Committing

Pre-commit hooks run automatically. If they fail:
```bash
# Fix issues and re-stage
git add <fixed-files>

# Or skip hooks (not recommended)
git commit --no-verify
```

### After Pushing

GitHub Actions automatically runs:
- All syntax checks
- Integration tests
- All linters
- Configuration validation

Check status at: `https://github.com/<owner>/<repo>/actions`

## 🧪 Testing Levels

### Level 1: Quick Validation (< 30 seconds)
```bash
# Config and syntax only
python3 -c "import yaml; yaml.safe_load(open('config/config-test.yml'))"
snakemake -n --configfile=config/config-test.yml
```

### Level 2: Local Tests (~ 2 minutes)
```bash
# Full local test suite
bash .tests/run_local_tests.sh
```

### Level 3: Integration Test (~ 10-15 minutes)
```bash
# Download test data (first time only)
bash .tests/dl_test_data.sh

# Setup tools (first time only)
snakemake setup_dorado dorado_model setup_modkit --cores 1 --configfile=config/config-test.yml

# Run non-GPU pipeline rules
snakemake merge_pods --cores 2 --configfile=config/config-test.yml
```

### Level 4: Full Pipeline (requires GPU)
```bash
# Run complete pipeline with test data
snakemake --cores 12 --configfile=config/config-test.yml

# Or submit to LSF cluster
bsub < run-test.sh
```

## 🔍 Continuous Integration Details

### What Gets Tested on Every Push

1. **Snakemake syntax**: Dry-run validation
2. **Config files**: YAML syntax validation
3. **Sample files**: Existence and content checks
4. **Directory structure**: Required directories present
5. **Code formatting**: Black, flake8, snakefmt
6. **Documentation**: Markdown linting

### What Gets Tested on Pull Requests

All of the above, plus:
- Integration test with test data
- Tool setup (dorado, modkit)
- Pipeline execution (non-GPU rules)

### What's NOT Tested in CI

Due to GitHub Actions limitations:
- GPU-intensive rules (rebasecall, classify_charging)
- Full end-to-end pipeline
- LSF cluster execution

These should be tested locally or on your cluster before merging.

## 🛠️ Maintenance

### Updating Dependencies

```bash
# Update conda environment
mamba env update -n aatrnaseqpipe -f workflow/envs/aatrnaseqpipe-env.yml

# Update pre-commit hooks
pre-commit autoupdate
```

### Adding New Tests

1. **Local tests**: Edit `.tests/run_local_tests.sh`
2. **CI tests**: Edit `.github/workflows/ci.yml`
3. **Linting**: Edit `.pre-commit-config.yaml` and `.github/workflows/lint.yml`

### Troubleshooting

**Pre-commit hooks failing?**
```bash
# Run manually to see detailed errors
pre-commit run --all-files

# Update hooks to latest versions
pre-commit autoupdate
```

**CI failing but local tests pass?**
- Check GitHub Actions logs for specific errors
- Ensure all files are committed
- Verify config files are valid YAML

**SessionStart hook not running?**
- Check if Claude Code session hooks are enabled
- Run manually: `bash .claude/hooks/SessionStart`
- Ensure file is executable: `chmod +x .claude/hooks/SessionStart`

## 📚 Additional Resources

- **Project Overview**: `CLAUDE.md`
- **CI Workflows**: `.github/workflows/`
- **Pre-commit Config**: `.pre-commit-config.yaml`
- **Test Scripts**: `.tests/`
- **Snakemake Docs**: https://snakemake.readthedocs.io/
