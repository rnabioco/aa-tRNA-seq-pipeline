# Contributing

Guidelines for contributing to the aa-tRNA-seq pipeline.

## Getting Started

### Fork and Clone

1. Fork the repository on GitHub
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/aa-tRNA-seq-pipeline.git
   cd aa-tRNA-seq-pipeline
   ```

3. Add upstream remote:
   ```bash
   git remote add upstream https://github.com/rnabioco/aa-tRNA-seq-pipeline.git
   ```

### Set Up Development Environment

```bash
# Install dependencies
pixi install

# Download test data
pixi run dl-test-data

# Set up tools
pixi run setup-tools
```

### Create a Branch

```bash
git checkout -b feature/my-new-feature
```

## Development Workflow

### Make Changes

1. Edit code in `workflow/` directory
2. Add tests if applicable
3. Update documentation

### Format Code

Format Snakemake files:

```bash
pixi run fmt
```

### Run Tests

```bash
# Dry run
pixi run dry-run

# Full test (requires GPU)
pixi run test
```

### Check DAG

Verify workflow integrity:

```bash
pixi run dag
```

## Code Style

### Snakemake Rules

```python
rule example_rule:
    """
    Brief description of what the rule does.
    Include any important notes about behavior.
    """
    input:
        bam=rules.previous_rule.output.bam,
    output:
        tsv=os.path.join(outdir, "path", "{sample}", "{sample}.output.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "example_rule", "{sample}"),
    params:
        src=SCRIPT_DIR,
        param_name=config.get("param", "default"),
    threads: 4
    shell:
        """
        command --input {input.bam} --output {output.tsv}
        """
```

### Python Scripts

Follow PEP 8 style guidelines:

```python
#!/usr/bin/env python
"""
Script description.

Usage:
    python script.py --input file.bam --output file.tsv
"""
import argparse
import pysam
import pandas as pd


def process_bam(input_path, output_path):
    """
    Process BAM file and write results.

    Args:
        input_path: Path to input BAM file.
        output_path: Path to output TSV file.
    """
    # Implementation
    pass


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    process_bam(args.input, args.output)


if __name__ == "__main__":
    main()
```

### Configuration

Use descriptive keys with comments:

```yaml
# Description of this section
section:
    # What this parameter controls
    parameter: value
```

## Documentation

### Update Documentation

When adding features:

1. Update relevant docs in `docs/`
2. Add docstrings to rules
3. Update README if needed

### Docstring Format

```python
rule my_rule:
    """
    Short description of the rule.

    Longer description if needed, explaining:
    - Input requirements
    - Output format
    - Any special behavior
    """
```

## Testing

### Local Testing

Run the test pipeline:

```bash
pixi run test
```

### Specific Rule Testing

Test a specific rule:

```bash
pixi run snakemake <rule_name> --cores 4 --configfile=config/config-test.yml
```

### Dry Run

Verify workflow without execution:

```bash
pixi run snakemake -n --configfile=config/config-test.yml
```

## Submitting Changes

### Commit Messages

Use clear, descriptive commit messages:

```
feat: Add read length distribution rule

- Add get_read_lengths.py script
- Add read_length_distribution rule to aatrnaseq-qc.smk
- Update pipeline_outputs() to include new output
```

Prefixes:

- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation only
- `refactor:` - Code refactoring
- `test:` - Test additions/changes
- `chore:` - Maintenance tasks

### Push Changes

```bash
git push origin feature/my-new-feature
```

### Create Pull Request

1. Go to GitHub
2. Click "New Pull Request"
3. Select your branch
4. Fill in the template:

```markdown
## Description
Brief description of changes.

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Refactoring

## Testing
How was this tested?

## Checklist
- [ ] Code follows style guidelines
- [ ] Snakemake files formatted (`pixi run fmt`)
- [ ] Documentation updated
- [ ] Tests pass locally
```

## Review Process

### What Reviewers Check

1. Code quality and style
2. Test coverage
3. Documentation
4. Compatibility with existing pipeline

### Responding to Feedback

1. Address all comments
2. Push additional commits
3. Request re-review when ready

## Release Process

Releases are managed by maintainers:

1. Update version numbers
2. Update CHANGELOG
3. Create release tag
4. Build and test

## Getting Help

### Questions

- Open a GitHub Discussion
- Check existing issues

### Bug Reports

File an issue with:

1. Environment details
2. Steps to reproduce
3. Expected vs actual behavior
4. Error messages/logs

### Feature Requests

File an issue with:

1. Use case description
2. Proposed solution
3. Alternatives considered

## Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Help others learn

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
