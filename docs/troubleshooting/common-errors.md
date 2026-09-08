# Common Errors

Solutions for frequently encountered errors in the aa-tRNA-seq pipeline.

## POD5 File Errors

### No POD5 Files Found

**Error:**
```
Error: No input files found for sample: sample1
```

**Causes:**

1. Incorrect path in samples file
2. POD5 files not in expected subdirectories
3. Wrong file extension

**Solutions:**

1. Verify the path exists:
   ```bash
   ls -la /path/from/samples/file/
   ```

2. Check for POD5 files in expected locations:
   ```bash
   ls /path/to/run/pod5_pass/
   ls /path/to/run/pod5_fail/
   ls /path/to/run/pod5/
   ```

3. Ensure files have `.pod5` extension (not `.fast5`)

### POD5 Merge Fails

**Error:**
```
pod5 merge: error
```

**Solutions:**

1. Check disk space:
   ```bash
   df -h .
   ```

2. Verify POD5 files are valid:
   ```bash
   pod5 inspect summary input.pod5
   ```

3. Check for corrupted files and exclude them

---

## GPU Errors

### CUDA Out of Memory

**Error:**
```
RuntimeError: CUDA out of memory
```

**Solutions:**

1. Ensure exclusive GPU access in cluster profile:
   ```yaml
   set-resources:
     - rebasecall:lsf_extra="-gpu num=1:j_exclusive=yes"
   ```

2. Reduce concurrent GPU jobs:
   ```yaml
   resources:
     - ngpu=4
   ```

3. Check for other GPU processes:
   ```bash
   nvidia-smi
   ```

### No CUDA GPUs Available

**Error:**
```
No CUDA GPUs are available
```

**Solutions:**

1. Verify CUDA installation:
   ```bash
   nvidia-smi
   ```

2. Check CUDA_VISIBLE_DEVICES:
   ```bash
   echo $CUDA_VISIBLE_DEVICES
   ```

3. Verify job is running on GPU node (for cluster execution)

### GPU Driver Mismatch

**Error:**
```
CUDA driver version is insufficient for CUDA runtime version
```

**Solutions:**

1. Check driver version:
   ```bash
   nvidia-smi | head -3
   ```

2. Update GPU drivers (contact system admin)

3. Use a compatible Dorado version

---

## Memory Errors

### Out of Memory

**Error:**
```
MemoryError
# or
Killed (signal 9)
```

**Solutions:**

1. Increase memory for the rule in cluster profile:
   ```yaml
   set-resources:
     - failing_rule:mem_mb=64
   ```

2. For local execution, close other applications

3. Check system memory:
   ```bash
   free -h
   ```

### Modkit Memory Issues

**Error:**
```
modkit extract calls: memory allocation failed
```

**Solution:**

The `modkit_extract_calls` rule requires significant memory (default 96 GB):

```yaml
set-resources:
  - modkit_extract_calls:mem_mb=128
```

---

## Alignment Errors

### BWA Index Missing

**Error:**
```
[bwa_idx_load_from_disk] fail to locate the index
```

**Solution:**

The index should be built automatically. If it fails, build manually:

```bash
bwa index resources/ref/sacCer3-mature-tRNAs-dual-adapt-v2.fa
```

### No Reads Aligned

**Error:**
```
Warning: 0 reads aligned
```

**Causes:**

1. Wrong reference sequence
2. Incompatible read format
3. Data quality issues

**Solutions:**

1. Verify reference matches your samples
2. Check FASTQ quality:
   ```bash
   zcat results/fq/sample.fq.gz | head -20
   ```

---

## Charging Classifier Errors

### Model Bundle Not Found

**Error:**
```
No such file or directory ... charging_feature_nn_sup6_rna004@v0.1.0
```

**Solution:**

`charging.model` is a **directory**, not a file. Verify it in config:

```yaml
charging:
  model: "resources/models/charging/charging_feature_nn_sup6_rna004@v0.1.0"
```

The bundle is vendored in the repository, so it should already be present:

```bash
pixi run verify-charging-model
```

### `missing field \`gbm\``

**Error:**
```
Error: missing field `gbm`
```

**Solution:**

The `escpod` on your PATH predates the per-base-feature bundle format. The
runtime and the model are pinned together — `escpod_version` must be >= 0.19.0 (the floor the charging bundle's `basecaller` block enforces).

```bash
pixi run setup
which escpod   # should be under resources/tools/escpod/<version>/bin
```

### Kmer Table Checksum Mismatch

**Error:**
```
kmer table sha256 does not match the bundle
```

**Solution:**

The bundle's `9mer_levels_v1.txt` is a symlink into `resources/kmers/`. If that
file was replaced, the residual feature is no longer the one the model was
trained against. Restore it, or give the bundle its own copy:

```bash
pixi run verify-charging-model
```

### Most Reads Get No `cl` Tag

This is usually not an error. The model **abstains** on reads whose common arm
did not align, and those reads carry no `cl` tag rather than a default class.
Check the rate and the reason:

```bash
zcat results/summary/tables/sample1/sample1.charging_calls.tsv.gz \
  | awk -F'\t' 'NR>1 {print ($5=="" ? "called" : $5)}' | sort | uniq -c
```

A high `no_aligned_arm` rate is a real signal, not a bug — but it also biases
the charging fraction low, so report it alongside. If instead nearly every read
is missing, check that the reference records carry `...CCA` followed by the
common arm `GGCTTCTTCTTGCTCTT`, and that `charging.min_mapq` is 0.

---

## Snakemake Errors

### Locked Directory

**Error:**
```
Directory cannot be locked
```

**Solution:**

```bash
pixi run snakemake --unlock --configfile=config/config.yml
```

### Missing Input Files

**Error:**
```
MissingInputException
```

**Solutions:**

1. Check if prerequisite rules completed
2. Verify file paths in config
3. Run a dry-run to check DAG:
   ```bash
   pixi run snakemake -n --configfile=config/config.yml
   ```

### Rule Failed

**Error:**
```
Error in rule <rule_name>
```

**Solution:**

1. Check the rule log:
   ```bash
   cat results/logs/<rule_name>/<sample>
   ```

2. Re-run with verbose output:
   ```bash
   pixi run snakemake -p --configfile=config/config.yml
   ```

---

## Configuration Errors

### Sample File Parse Error

**Error:**
```
Error parsing samples file
```

**Solutions:**

1. For TSV: Ensure tab-separated (not spaces)
2. For YAML: Check indentation
3. Validate with:
   ```bash
   # TSV
   cat -A config/samples.tsv  # Shows ^I for tabs

   # YAML
   python -c "import yaml; yaml.safe_load(open('config/samples.yml'))"
   ```

### Config Key Missing

**Error:**
```
KeyError: '<key>'
```

**Solution:**

Ensure your config inherits from base:

```yaml
# Your config should be used with --configfile
# Base config is loaded automatically by Snakefile
```

---

## WarpDemuX Errors

### WarpDemuX Not Found

**Error:**
```
warpdemux: command not found
```

**Solution:**

WarpDemuX is opt-in since v0.3.0, so `pixi run setup` does not install it.
Install it into its own environment:

```bash
pixi install -e warpdemux
pixi run -e warpdemux install-warpdemux
```

### Invalid Barcode Kit

**Error:**
```
Invalid barcode kit name
```

**Solution:**

Use a valid kit name:

- `WDX4_tRNA_rna004_v1_0`
- `WDX4b_tRNA_rna004_v1_0`

### No Reads for Barcode

**Error:**
Sample has 0 reads after demultiplexing.

**Solutions:**

1. Verify barcode assignment in YAML file
2. Check demux summary:
   ```bash
   zcat results/demux/read_ids/run_id/demux_summary.tsv.gz
   ```
3. Ensure barcode kit matches library prep

---

## File System Errors

### Disk Full

**Error:**
```
No space left on device
```

**Solutions:**

1. Check disk usage:
   ```bash
   df -h .
   ```

2. Clean intermediate files:
   ```bash
   rm -rf results/bam/rebasecall results/bam/aln results/fq
   ```

3. Use a different output directory

### Permission Denied

**Error:**
```
PermissionError: [Errno 13] Permission denied
```

**Solutions:**

1. Check file permissions:
   ```bash
   ls -la <file>
   ```

2. Verify write access to output directory
