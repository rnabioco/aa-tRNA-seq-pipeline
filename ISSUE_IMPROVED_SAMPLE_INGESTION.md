# Issue: Improve sample ingestion to support dorado demux and WarpDemuX workflows

## Problem Statement

The current sample ingestion approach is limited to raw pod5 files in a rigid directory structure. This doesn't accommodate modern demultiplexing workflows like dorado demux and WarpDemuX (especially WarpDemuX-tRNA), which are increasingly important for multiplexed tRNA sequencing experiments.

## Current Limitations

The existing implementation (`workflow/rules/common.smk:10-87`) has several constraints:

1. **Rigid TSV format**: Only accepts 2 columns (sample_id, data_path)
2. **Raw pod5 only**: Expects specific directory structure (pod5_pass/pod5_fail/pod5 subdirectories)
3. **No demux support**: Cannot handle pre-demultiplexed data from dorado or WarpDemuX
4. **Always rebasecalls**: Cannot skip rebasecalling even when data is already basecalled
5. **No barcode awareness**: Cannot selectively process specific barcodes from multiplexed runs
6. **No metadata capture**: Cannot track experimental conditions or sample relationships

## Proposed Solution

### 1. Flexible CSV Sample Sheet Format

Replace the 2-column TSV with a comprehensive CSV format:

```csv
sample_id,input_type,input_path,barcode,flow_cell_id,experiment_id,metadata
sample1,raw_pod5,/path/to/run1,,,231212_exp1,
sample2,raw_pod5,/path/to/run2,,,231212_exp2,condition=treated
sample3,dorado_demux,/path/to/demux_output,SQK-NBD114-96_barcode01,PAO12345,231212_exp3,
sample4,dorado_demux,/path/to/demux_output,SQK-NBD114-96_barcode02,PAO12345,231212_exp3,
sample5,warpdemux,/path/to/warpdemux_output,1,,,
sample6,merged_pod5,/path/to/merged.pod5,,,231212_exp4,
sample7,basecalled_bam,/path/to/basecalled.bam,,,231212_exp5,
```

**Column Definitions:**

- **sample_id** (required): Unique sample identifier
- **input_type** (required): One of:
  - `raw_pod5`: Raw sequencing run directory (current behavior)
  - `dorado_demux`: Pre-demultiplexed BAM files from dorado demux
  - `warpdemux`: WarpDemuX prediction CSV output
  - `merged_pod5`: Pre-merged single pod5 file
  - `basecalled_bam`: Pre-basecalled BAM file with move tables
- **input_path** (required): Path to input data
- **barcode** (optional): Barcode identifier for demuxed samples
- **flow_cell_id** (optional): Flow cell ID for tracking
- **experiment_id** (optional): Experiment identifier
- **metadata** (optional): Key=value pairs for additional metadata

### 2. Implementation Architecture

#### A. Enhanced Sample Parsing

```python
def parse_samples_v2(sample_sheet):
    """
    Parse flexible CSV sample sheet supporting multiple input types.

    Returns:
        Dict structure:
        {
            'sample1': {
                'input_type': 'raw_pod5',
                'input_path': '/path/to/run1',
                'barcode': None,
                'flow_cell_id': 'PAO12345',
                'experiment_id': '231212_exp1',
                'metadata': {},
                'raw_files': [...],
                'start_rule': 'merge_pods'  # Entry point in workflow
            }
        }
    """
```

#### B. Input Type Handlers

Each input type needs a specific handler:

```python
def find_raw_inputs_v2(sample_dict):
    """Route to appropriate handler based on input_type"""
    handlers = {
        'raw_pod5': handle_raw_pod5,
        'dorado_demux': handle_dorado_demux,
        'warpdemux': handle_warpdemux,
        'merged_pod5': handle_merged_pod5,
        'basecalled_bam': handle_basecalled_bam
    }

    for sample, info in sample_dict.items():
        handler = handlers[info['input_type']]
        handler(sample, info)

    return sample_dict
```

**Handler Functions:**

- `handle_raw_pod5()`: Current behavior (find pod5 files in subdirectories)
- `handle_dorado_demux()`: Find demuxed BAM file matching barcode, skip rebasecalling
- `handle_warpdemux()`: Parse WarpDemuX predictions CSV and filter pod5 by barcode
- `handle_merged_pod5()`: Single pre-merged pod5 file, start at rebasecall
- `handle_basecalled_bam()`: Pre-basecalled BAM, start at ubam_to_fastq

#### C. Workflow Entry Point Logic

Modify rules to support conditional entry points:

```python
def get_fastq_input(wildcards):
    """
    Determine input for ubam_to_fastq based on input_type.
    """
    input_type = samples[wildcards.sample]['input_type']

    if input_type in ['raw_pod5', 'merged_pod5', 'warpdemux']:
        return rules.rebasecall.output
    elif input_type in ['dorado_demux', 'basecalled_bam']:
        return samples[wildcards.sample]['input_path']
```

### 3. New Rules for Demux Workflows

#### A. WarpDemuX Integration

```python
rule filter_warpdemux_pods:
    """
    Filter pod5 files based on WarpDemuX barcode predictions.
    """
    input:
        predictions=lambda wc: samples[wc.sample]['predictions_file'],
        pod5_dir=lambda wc: samples[wc.sample]['pod5_dir']
    output:
        os.path.join(outdir, "pod5", "{sample}", "{sample}.filtered.pod5")
    params:
        barcode=lambda wc: samples[wc.sample]['barcode']
    shell:
        """
        python {SCRIPT_DIR}/filter_warpdemux_reads.py \
            --predictions {input.predictions} \
            --pod5-dir {input.pod5_dir} \
            --barcode {params.barcode} \
            --output {output}
        """
```

#### B. Dorado Demux Input

```python
rule link_dorado_demux:
    """
    Create symbolic link to dorado demuxed BAM.
    """
    input:
        lambda wc: samples[wc.sample]['demux_bam']
    output:
        os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam")
    shell:
        "ln -s {input} {output}"
```

### 4. Backward Compatibility

Auto-detect sample sheet format:

```python
def parse_samples(sample_file):
    """
    Auto-detect sample sheet format and parse accordingly.
    - .csv extension: use parse_samples_v2()
    - .tsv extension: use legacy format
    """
    ext = os.path.splitext(sample_file)[1]

    if ext == '.csv':
        return parse_samples_v2(sample_file)
    else:
        # Legacy TSV: convert to v2 structure with defaults
        samples = parse_samples_legacy(sample_file)
        for sample, info in samples.items():
            info.update({
                'input_type': 'raw_pod5',
                'barcode': None,
                'start_rule': 'merge_pods'
            })
        return samples
```

### 5. Helper Scripts

**New scripts to implement:**

1. **`workflow/scripts/filter_warpdemux_reads.py`**
   - Parse WarpDemuX predictions CSV
   - Filter read IDs by barcode confidence threshold
   - Extract matching reads from pod5 files using pod5 Python API

2. **`workflow/scripts/validate_sample_sheet.py`**
   - Validate sample sheet format
   - Check required columns
   - Verify file/directory paths exist
   - Check for duplicate sample_ids
   - Validate barcode formats

3. **`workflow/scripts/convert_dorado_samplesheet.py`**
   - Convert dorado demux sample sheet to pipeline format
   - Map aliases to sample_ids
   - Auto-detect demuxed BAM file paths

### 6. Configuration Updates

Add to `config/config-base.yml`:

```yaml
# Sample sheet format version
sample_sheet_version: 2

# Demux-related options
demux:
  # WarpDemuX confidence threshold (0-1)
  warpdemux_confidence: 0.99

  # Skip rebasecalling for pre-basecalled inputs
  skip_rebasecall_when_possible: true

  # Dorado demux file naming pattern
  dorado_demux_pattern: "{barcode}.bam"
```

## Implementation Plan

### Phase 1: Core Infrastructure
- [ ] Implement `parse_samples_v2()` with CSV parsing
- [ ] Add input type handler framework
- [ ] Implement backward compatibility layer
- [ ] Add unit tests for sample parsing

### Phase 2: Input Type Handlers
- [ ] Implement `handle_raw_pod5()` (refactor existing code)
- [ ] Implement `handle_merged_pod5()`
- [ ] Implement `handle_basecalled_bam()`
- [ ] Add conditional input functions to rules

### Phase 3: Demux Support
- [ ] Implement `handle_dorado_demux()`
- [ ] Add `link_dorado_demux` rule
- [ ] Create `convert_dorado_samplesheet.py` script
- [ ] Test with dorado demux data

### Phase 4: WarpDemuX Support
- [ ] Implement `handle_warpdemux()`
- [ ] Add `filter_warpdemux_pods` rule
- [ ] Create `filter_warpdemux_reads.py` script
- [ ] Test with WarpDemuX-tRNA data

### Phase 5: Validation & Documentation
- [ ] Create `validate_sample_sheet.py` script
- [ ] Add validation to workflow startup
- [ ] Update CLAUDE.md with new sample sheet format
- [ ] Create example sample sheets for each input type
- [ ] Add migration guide from TSV to CSV format

## Benefits

1. **Flexibility**: Supports multiple workflow entry points
2. **Efficiency**: Skips unnecessary steps (e.g., rebasecalling when already done)
3. **Dorado Integration**: Native support for dorado demux sample sheets
4. **WarpDemuX Support**: Handles tRNA-optimized demultiplexing
5. **Backward Compatible**: Existing TSV files continue to work
6. **Extensible**: Easy to add new input_type handlers
7. **Metadata**: Captures experimental metadata for downstream analysis
8. **Validation**: Sample sheets validated before workflow execution
9. **Resource Efficiency**: Can process only specific barcodes from large multiplexed runs

## Testing Strategy

Create test configs for each input type:
- `config/samples-test-raw.csv` (current behavior)
- `config/samples-test-dorado-demux.csv`
- `config/samples-test-warpdemux.csv`
- `config/samples-test-mixed.csv` (multiple input types in one run)
- `config/samples-test-legacy.tsv` (backward compatibility)

## References

- [Dorado Sample Sheets Documentation](https://software-docs.nanoporetech.com/dorado/latest/barcoding/sample_sheet/)
- [WarpDemuX GitHub](https://github.com/KleistLab/WarpDemuX)
- [WarpDemuX-tRNA Preprint](https://www.biorxiv.org/content/10.1101/2025.03.21.644602v1)
