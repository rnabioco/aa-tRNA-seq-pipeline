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
escpod merge: error
```

**Solutions:**

1. Check disk space:
   ```bash
   df -h .
   ```

2. Verify POD5 files are valid:
   ```bash
   escpod inspect summary input.pod5
   ```

3. Check for corrupted files and exclude them

!!! note "Interrupted merges are safe"
    `escpod` stages output to a temp file and renames it into place, so an interrupted
    or killed job never leaves a partially written POD5 at the destination. A truncated
    output file is not an expected failure mode; look at disk space and the inputs
    instead.

### Empty Read-ID File

**Error:**
```
No read IDs found
```

**Cause:**

`escpod filter` (used by `split_pod5` and `filter_pod5_by_edx`) treats an empty read-ID
list as a hard error. This happens when a barcode or an EDX 3' adapter matched zero
reads. The ONT `pod5 filter --missing-ok` this replaced would instead have written an
empty POD5, so this is a behavior change: the job now fails rather than silently
producing an empty file that downstream rules process.

**Solutions:**

1. Check the barcode distribution for the run:
   ```bash
   zcat results/demux/read_ids/{run_id}/demux_summary.tsv.gz
   ```

2. For EDX samples, check the adapter detection table:
   ```bash
   zcat results/demux/edx/{sample}/{sample}.edx_adapters.tsv.gz | cut -f2 | sort | uniq -c
   ```

3. Fix the barcode / `edx` assignment in the samples file, or drop the sample

!!! note "Missing read IDs are only a warning"
    `escpod filter` has no `--missing-ok` flag because it does not need one — read IDs
    present in the list but absent from the input POD5 produce a warning, not an error.
    Only a completely empty list fails.

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

## Remora Errors

### Model Not Found

**Error:**
```
FileNotFoundError: remora model not found
```

**Solution:**

Verify the model path in config:

```yaml
remora_cca_classifier: "resources/models/cca_classifier.pt"
```

Ensure the file exists:

```bash
ls -la resources/models/cca_classifier.pt
```

### Kmer Table Error

**Error:**
```
Error loading kmer table
```

**Solution:**

Verify kmer table path:

```yaml
remora_kmer_table: "resources/kmers/9mer_levels_v1.txt"
```

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

## Leech Errors

### Leech Not Found

**Error:**
```
leech: command not found
```

**Solution:**

leech is installed from the release wheels of the private `rnabioco/leech` repo, which
needs an authenticated `gh`:

```bash
gh auth login
pixi run install-leech
```

Or install from a local directory of wheels:

```bash
LEECH_WHEEL_DIR=/path/to/wheels pixi run install-leech
```

leech is only needed for `classifier: leech` and the amino-acid classification rules —
the default Remora path does not require it. Note that there is no longer a
`resources/leech` git submodule; do not run `git submodule update`.

### Leech Requires Python 3.12

**Error:**

Wheel install fails, or `import leech` raises a syntax/ABI error on an older
interpreter.

**Solution:**

leech requires Python >= 3.12. Check the environment's interpreter:

```bash
pixi run python --version
```

### `--backend rust` Unavailable / Slow Extraction

**Warning during setup:**
```
no leech_core wheel for cpXY/<arch>; leech will use the slower pure-python backend
```

`leech-core` is a separate optional Rust extension wheel published per interpreter and
architecture. Without it, leech still works but falls back to a slower pure-python
extraction backend and `--backend rust` is unavailable. Check that
`leech_core_version` in `config/config-base.yml` has a wheel matching your Python
version and architecture in the release.

---

## Demultiplexing Errors

### `escpod demux` Unavailable

**Error:**
```
error: unrecognized subcommand 'demux'
```

or, at pipeline start:
```
'escpod demux' is unavailable — this escpod was built without the demux feature
```

**Cause:**

The installed `escpod` is a default-features build. All published escapepod-rs release
binaries are default-features only and do **not** include a working `escpod demux`.

**Solution:**

Rebuild from source (requires a Rust toolchain >= 1.95):

```bash
pixi run install-escpod
```

Or use the python backend instead:

```yaml
warpdemux:
    backend: "warpdemux"
```

### Barcode Model Not Found

**Error:**
```
warpdemux.barcode_model is not set. Run 'bash scripts/install-demux-models.sh' and check config.
```

or a missing `resources/models/demux/barcode_wdx4_rna004.gbm.json`.

**Solution:**

```bash
pixi run install-demux-models
```

The barcode GBM models are not yet published as releases on `rnabioco/escapepod-models`
(only `adapter_rna004@v1.0.1` is), so a local checkout is currently required:

```bash
ESCAPEPOD_MODELS_DIR=/path/to/escapepod-models pixi run install-demux-models
```

With `warpdemux.method: cnn` (the default, and the setting the barcode model was trained
behind), `warpdemux.adapter_model` must also point at an existing ONNX file.

### WarpDemuX Not Found

**Error:**
```
warpdemux: command not found
```

**Solution:**

Only the `warpdemux` backend needs the python package. Install it:

```bash
pixi run install-warpdemux   # or: pixi run setup
```

Or use the default `escpod` backend.

### Invalid Barcode Kit

**Error:**
```
Invalid barcode kit name
```

**Solution:**

`barcode_kit` is only read by the `warpdemux` backend. Use a valid kit name:

- `WDX4_tRNA_rna004_v1_0`
- `WDX4b_tRNA_rna004_v1_0`

On the `escpod` backend the barcode set comes from `warpdemux.barcode_model` instead; the
shipped `barcode_wdx4_rna004` GBM covers barcodes 03, 04, 05 and 07. If a sample's
barcode is not one of the model's classes, `collect_escpod_pod5` fails with a pointer to
`warpdemux.barcode_model`.

### No Reads for Barcode

**Error:**
Sample has 0 reads after demultiplexing (on the `warpdemux` backend this surfaces as
`escpod filter`'s "No read IDs found" — see [Empty Read-ID File](#empty-read-id-file)).

**Solutions:**

1. Verify barcode assignment in YAML file
2. Check demux summary:
   ```bash
   zcat results/demux/read_ids/run_id/demux_summary.tsv.gz
   ```
3. Ensure barcode kit / barcode model matches library prep

### Unclassified Fraction Much Smaller Than Expected

Expected on the `escpod` backend: the shipped barcode GBM carries no per-class confidence
thresholds, so every read with a usable adapter boundary is assigned a barcode, whereas
WarpDemuX rejected low-confidence reads to `unclassified`. Filter on the `confidence`
column of `demux/escpod_output/{run_id}/classifications.csv` to restore WarpDemuX-like
rejection. See [Demultiplexing](../workflow/demultiplexing.md#choosing-a-backend).

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
