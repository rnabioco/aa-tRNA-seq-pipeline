# Claude Code Configuration

This directory contains configuration and hooks for Claude Code development sessions.

## Session Start Hook

The `hooks/SessionStart` script runs automatically when you start a Claude Code session. It:

- Verifies conda/mamba is available
- Checks if the `aatrnaseqpipe` environment exists (creates it if not)
- Activates the environment
- Validates Snakemake installation
- Validates configuration files
- Sets up pre-commit hooks
- Runs a quick Snakemake syntax check

This ensures your development environment is properly configured before starting work.

## Manual Testing

You can run the session start hook manually:

```bash
bash .claude/hooks/SessionStart
```

## Development Workflow

1. **Before making changes**: Run `bash .tests/run_local_tests.sh` to verify everything works
2. **During development**: Use `snakemake -n --configfile=config/config-test.yml` for syntax checks
3. **Before committing**: Pre-commit hooks will automatically run linting and formatting
4. **After committing**: GitHub Actions will run full CI/CD checks

## Learn More

- Project documentation: `CLAUDE.md`
- Local test script: `.tests/run_local_tests.sh`
- CI/CD workflows: `.github/workflows/`
- Pre-commit config: `.pre-commit-config.yaml`
