# SLURM Setup

Configure the pipeline for SLURM cluster execution.

## Overview

The pipeline includes a pre-configured SLURM profile at `cluster/slurm/config.yaml` that handles:

- Job submission via the native Snakemake SLURM executor
- Memory and partition allocation
- GPU job routing
- Runtime limits

## Quick Start

```bash
# Run with SLURM profile
pixi run snakemake --profile cluster/slurm --configfile=config/config.yml

# Or use the test shortcut
pixi run test-slurm
```

## Profile Configuration

### Location

```
cluster/slurm/config.yaml
```

### Full Configuration

```yaml
executor: slurm
jobs: 300
latency-wait: 15

slurm-logdir: "logs/slurm"
slurm-keep-successful-logs: false
slurm-delete-logfiles-older-than: 10

default-resources:
  slurm_partition: "cpu"
  slurm_account: "aatrnaseq"
  runtime: 120
  mem_mb: 8000
  cpus_per_task: 1

set-resources:
  rebasecall:
    slurm_partition: "gpu"
    gres: "gpu:1"
    runtime: 480
    mem_mb: 24000
    cpus_per_task: 4

  classify_charging:
    slurm_partition: "gpu"
    gres: "gpu:1"
    runtime: 240
    mem_mb: 24000
    cpus_per_task: 4

  bwa_align:
    runtime: 240
    mem_mb: 24000
    cpus_per_task: 8

  remora_signal_stats:
    runtime: 180
    mem_mb: 24000
    cpus_per_task: 4

  modkit_extract_calls:
    runtime: 360
    mem_mb: 96000
    cpus_per_task: 4

  warpdemux:
    runtime: 360
    mem_mb: 32000
    cpus_per_task: 4

  parse_warpdemux:
    runtime: 60
    mem_mb: 8000
    cpus_per_task: 1

  merge_pods_for_demux:
    runtime: 120
    mem_mb: 16000
    cpus_per_task: 1

rerun-incomplete: true
keep-going: true
printshellcmds: true
show-failed-logs: true
```

## Configuration Options

### Global Settings

| Option | Value | Description |
|--------|-------|-------------|
| `executor` | `slurm` | Use native SLURM executor |
| `jobs` | `300` | Maximum concurrent jobs |
| `latency-wait` | `15` | Seconds to wait for file sync |

### Default Resources

Applied to all rules unless overridden:

| Resource | Value | Description |
|----------|-------|-------------|
| `slurm_partition` | `cpu` | Default partition (customize for your cluster) |
| `slurm_account` | `aatrnaseq` | Account for job submission |
| `runtime` | `120` | Default runtime in minutes |
| `mem_mb` | `8000` | Memory in MB (8GB) |
| `cpus_per_task` | `1` | CPUs per task |

## Per-Rule Resources

### GPU Rules

These rules are automatically submitted to the GPU partition:

| Rule | Partition | Memory | Runtime | GPU |
|------|-----------|--------|---------|-----|
| `rebasecall` | gpu | 24 GB | 8 hours | 1 |
| `classify_charging` | gpu | 24 GB | 4 hours | 1 |

### Memory-Intensive Rules

| Rule | Memory | Runtime |
|------|--------|---------|
| `modkit_extract_calls` | 96 GB | 6 hours |
| `warpdemux` | 32 GB | 6 hours |
| `remora_signal_stats` | 24 GB | 3 hours |
| `bwa_align` | 24 GB | 4 hours |
| `merge_pods_for_demux` | 16 GB | 2 hours |

## Customization

### Change Default Partition

Edit `slurm_partition` in default-resources:

```yaml
default-resources:
  slurm_partition: "your_partition"
```

### Change Account

For job accounting:

```yaml
default-resources:
  slurm_account: "your_account"
```

### Adjust Max Jobs

Reduce if you're hitting queue limits:

```yaml
jobs: 100
```

### Increase Memory for a Rule

Add or modify in `set-resources`:

```yaml
set-resources:
  your_rule:
    mem_mb: 64000
```

### Change GPU Partition Name

If your GPU partition has a different name:

```yaml
set-resources:
  rebasecall:
    slurm_partition: "your_gpu_partition"
  classify_charging:
    slurm_partition: "your_gpu_partition"
```

### Request Specific GPU Type

For clusters with multiple GPU types:

```yaml
set-resources:
  rebasecall:
    gres: "gpu:v100:1"
```

## Monitoring Jobs

### View Your Jobs

```bash
squeue -u $USER
```

### View Job Details

```bash
scontrol show job <job_id>
```

### Cancel a Job

```bash
scancel <job_id>
```

### Cancel All Your Jobs

```bash
scancel -u $USER
```

### View Partition Status

```bash
sinfo
```

### View Job History

```bash
sacct -j <job_id>
```

## Submit Scripts

For long-running pipelines, use a submit script:

=== "run-slurm.sh"

    ```bash
    #!/bin/bash
    #SBATCH --job-name=aa-tRNA-seq
    #SBATCH --output=logs/pipeline.%j.out
    #SBATCH --error=logs/pipeline.%j.err
    #SBATCH --partition=cpu
    #SBATCH --mem=4G
    #SBATCH --cpus-per-task=1
    #SBATCH --time=24:00:00

    mkdir -p logs

    pixi run snakemake --profile cluster/slurm \
        --configfile=config/config.yml
    ```

Submit:

```bash
sbatch run-slurm.sh
```

## Troubleshooting

### Jobs Pending Too Long

Check partition limits:

```bash
sinfo -p cpu
```

Reduce concurrent jobs:

```yaml
jobs: 50
```

### Memory Errors

Increase memory for the failing rule:

```yaml
set-resources:
  failing_rule:
    mem_mb: 128000
```

### GPU Jobs Not Starting

Check GPU partition availability:

```bash
sinfo -p gpu
```

Verify GPU resource specification for your cluster:

```yaml
set-resources:
  rebasecall:
    gres: "gpu:1"
```

### File Sync Errors

Increase latency wait:

```yaml
latency-wait: 60
```

### Runtime Exceeded

Increase runtime for the failing rule (in minutes):

```yaml
set-resources:
  slow_rule:
    runtime: 720  # 12 hours
```

## Comparison with LSF

| Feature | LSF | SLURM |
|---------|-----|-------|
| Executor | `lsf` | `slurm` |
| Queue/Partition | `lsf_queue` | `slurm_partition` |
| Account | `lsf_project` | `slurm_account` |
| Memory units | GB (e.g., `24`) | MB (e.g., `24000`) |
| GPU request | `lsf_extra="-gpu num=1"` | `gres: "gpu:1"` |
| Runtime | Not specified | `runtime` (minutes) |

## Next Steps

- [LSF Setup](lsf-setup.md) - For LSF clusters
- [GPU Configuration](gpu-configuration.md) - GPU-specific settings
