# Adding Rules

Guide for adding new Snakemake rules to the pipeline.

## Overview

When adding new functionality:

1. Choose the appropriate rule file
2. Create the rule definition
3. Add any required scripts
4. Update output aggregation
5. Test thoroughly

## Choosing a Rule File

| Module | Use For |
|--------|---------|
| `aatrnaseq-process.smk` | Core data processing (POD5 → BAM) |
| `aatrnaseq-charging.smk` | Charging analysis rules |
| `aatrnaseq-qc.smk` | Quality control and statistics |
| `aatrnaseq-modifications.smk` | Modification calling and coverage |
| `warpdemux.smk` | Demultiplexing (conditional) |

## Rule Template

Basic rule structure:

```python
rule my_new_rule:
    """
    Brief description of what this rule does.
    """
    input:
        bam=rules.finalize_bam.output.bam,
        # Or use path patterns
        # bam=os.path.join(outdir, "bam", "final", "{sample}.bam"),
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.my_output.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "my_new_rule", "{sample}"),
    params:
        # Pass config values
        some_param=config.get("my_param", "default"),
        # Pass script directory for Python scripts
        src=SCRIPT_DIR,
    threads: 4  # Optional: number of threads
    shell:
        """
        python {params.src}/my_script.py \
            --input {input.bam} \
            --output {output.tsv} \
            --param {params.some_param}
        """
```

## Input Patterns

### Reference Previous Rules

```python
input:
    bam=rules.finalize_bam.output.bam,
    bai=rules.finalize_bam.output.bai,
```

### Use Path Patterns

```python
input:
    bam=os.path.join(outdir, "bam", "final", "{sample}.bam"),
```

### Use Input Functions

```python
def get_my_input(wildcards):
    return some_logic(wildcards.sample)

rule my_rule:
    input:
        get_my_input,
```

## Output Patterns

### Standard Output

```python
output:
    tsv=os.path.join(
        outdir, "summary", "tables", "{sample}", "{sample}.output.tsv.gz"
    ),
```

### Protected Output

Prevent accidental deletion:

```python
output:
    protected(os.path.join(outdir, "bam", "final", "{sample}.bam")),
```

### Temporary Output

Deleted after downstream rules complete:

```python
output:
    temp(os.path.join(outdir, "temp", "{sample}.tmp")),
```

### Directory Output

For rules that produce directories:

```python
output:
    directory(os.path.join(outdir, "my_dir", "{sample}")),
```

## Adding Python Scripts

### Create the Script

Place in `workflow/scripts/`:

```python
#!/usr/bin/env python
"""
my_script.py - Brief description.
"""
import argparse
import pysam
import pandas as pd

def main(input_bam, output_tsv, param):
    """Main processing function."""
    # Your logic here
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--param", default="default")
    args = parser.parse_args()

    main(args.input, args.output, args.param)
```

### Call from Rule

```python
params:
    src=SCRIPT_DIR,
shell:
    """
    python {params.src}/my_script.py \
        --input {input.bam} \
        --output {output.tsv}
    """
```

## Using `run:` Directive

For simple Python logic inline:

```python
rule process_inline:
    input:
        tsv="input.tsv",
    output:
        tsv="output.tsv",
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t")
        df["new_column"] = df["old_column"] * 2
        df.to_csv(output.tsv, sep="\t", index=False)
```

## Updating Output Aggregation

If your rule produces a final output, add to `pipeline_outputs()` in `common.smk`:

```python
def pipeline_outputs():
    outputs = []
    for sample in samples:
        outputs.extend([
            # Existing outputs...
            f"{outdir}/summary/tables/{sample}/{sample}.charging.cpm.tsv.gz",

            # Add your new output
            f"{outdir}/summary/tables/{sample}/{sample}.my_output.tsv.gz",
        ])
    return outputs
```

## Adding Configuration Options

### Add to config-base.yml

```yaml
# My new feature
my_feature:
    enabled: true
    threshold: 0.5
```

### Access in Rules

```python
params:
    threshold=config.get("my_feature", {}).get("threshold", 0.5),
```

## GPU Rules

For rules requiring GPU:

```python
rule gpu_rule:
    input:
        pod5="pod5/{sample}.pod5",
    output:
        bam="output/{sample}.bam",
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            export CUDA_VISIBLE_DEVICES
        fi

        my_gpu_tool --input {input.pod5} --output {output.bam}
        """
```

Add cluster resources:

```yaml
# In cluster/lsf/config.yaml
set-resources:
  - gpu_rule:lsf_queue="gpu"
  - gpu_rule:lsf_extra="-gpu num=1:j_exclusive=yes"
  - gpu_rule:ngpu=1
  - gpu_rule:mem_mb=24
```

## Conditional Rules

### Based on Config

```python
if config.get("my_feature", {}).get("enabled", False):
    rule optional_rule:
        ...
```

### Based on Input Existence

```python
def get_conditional_input(wildcards):
    if some_condition():
        return "optional_input.txt"
    return []

rule conditional:
    input:
        required="required.txt",
        optional=get_conditional_input,
```

## Testing New Rules

### Dry Run

```bash
pixi run snakemake -n --configfile=config/config-test.yml my_new_rule
```

### Run Specific Rule

```bash
pixi run snakemake --cores 4 --configfile=config/config-test.yml \
    results/summary/tables/sample1/sample1.my_output.tsv.gz
```

### Force Rerun

```bash
pixi run snakemake --forcerun my_new_rule --configfile=config/config-test.yml
```

### Print Commands

```bash
pixi run snakemake -n -p --configfile=config/config-test.yml my_new_rule
```

## Best Practices

### Rule Design

1. **Single responsibility**: Each rule does one thing
2. **Clear inputs/outputs**: Use explicit paths
3. **Proper logging**: Always include `log:` directive
4. **Document**: Add docstring describing purpose

### Output Paths

1. Use `os.path.join()` for cross-platform compatibility
2. Include sample wildcard in paths
3. Follow existing directory structure

### Error Handling

1. Check input validity in scripts
2. Provide informative error messages
3. Clean up on failure (use temp files)

### Performance

1. Set appropriate `threads:` for parallelizable tasks
2. Use `temp()` for intermediate files
3. Consider memory requirements for cluster config

## Example: Complete New Rule

```python
# In aatrnaseq-qc.smk

rule read_length_distribution:
    """
    Calculate read length distribution from final BAM.
    """
    input:
        bam=rules.finalize_bam.output.bam,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.read_lengths.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "read_length_distribution", "{sample}"),
    params:
        src=SCRIPT_DIR,
    shell:
        """
        python {params.src}/get_read_lengths.py \
            {input.bam} \
            {output.tsv} \
            2>&1 | tee {log}
        """
```

Then add to `pipeline_outputs()`:

```python
f"{outdir}/summary/tables/{sample}/{sample}.read_lengths.tsv.gz",
```

## Next Steps

- [Architecture](architecture.md) - Understand the codebase
- [Contributing](contributing.md) - Submit your changes
