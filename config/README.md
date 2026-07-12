# Configuring the pipeline with config.yml

Edit config.yml to specify the following parameters.

## Sample File Configuration

The pipeline supports two sample file formats depending on whether you need barcode demultiplexing.

### Standard TSV Format (No Demultiplexing)

For non-multiplexed sequencing runs, use a two-column TSV file:

```
sample1    /path/to/run1
sample2    /path/to/run2
sample2    /path/to/run2_replicate
```

- Column 1: Unique sample ID
- Column 2: Path to sequencing run folder containing `pod5_pass`, `pod5`, `pod5_fail`, `fast5_pass`, or `fast5_fail` subdirectories

If multiple rows share the same sample ID, the reads will be merged before processing. See `samples-test.tsv` for an example.

### YAML Format (With WarpDemuX Demultiplexing)

For multiplexed/pooled sequencing runs using WarpDemuX barcodes, use a YAML file:

```yaml
runs:
  - path: /path/to/pooled/sequencing/run
    barcode_kit: "WDX4_tRNA_rna004_v1_0"  # optional, uses config default if omitted
    samples:
      charged_sample: "barcode03"
      uncharged_sample: "barcode04"
      control_sample: "barcode05"

  - path: /path/to/another/pooled/run
    samples:
      experimental_bc03: "barcode03"
      experimental_bc04: "barcode04"

  # Non-multiplexed run within a demux config (skip demultiplexing)
  - path: /path/to/non-pooled/run
    samples:
      direct_sample: ~  # null barcode skips demultiplexing
```

See `samples-demux-example.yml` for a complete example with comments.

**Available barcode kits for Nano-tRNAseq:**
- `WDX4_tRNA_rna004_v1_0` (recommended) with barcodes: `barcode03`, `barcode04`, `barcode05`, `barcode07`

**Note:** WarpDemuX-tRNA models do NOT work with Thomas splint adapter data.

### Enabling WarpDemuX Demultiplexing

To use demultiplexing, add the following to your config file:

```yaml
samples: config/samples-demux.yml  # YAML format sample file

warpdemux:
    enabled: true
    barcode_kit: "WDX4_tRNA_rna004_v1_0"  # default kit if not specified per-run
    save_boundaries: true  # optional, saves adapter boundary information
    threads: 8
```

Run with the demux environment:

```bash
pixi run snakemake --configfile=config/config-demux.yml --cores 8
```

See `config-demux-test.yml` for a complete example.

## Other Configuration Parameters

- `base_calling_model`: Path to the dorado basecalling model to use for rebasecalling. We use `rna004_130bps_sup@v5.0.0` for now, will evaluate newer model soon.

- `input_format`: A string, either "FAST5" or "POD5". If FAST5, files will be converted to POD5 before rebasecalling.

- `output_directory`: Path where pipeline outputs will be written.

- `cleanup_intermediates`: Controls automatic `temp()` deletion of large,
  regenerable intermediates *during* a run. Accepts a boolean or a list of tier
  names (default: off / opt-in):
  - `false` (or omitted): nothing is auto-deleted (all intermediates retained).
  - `true`: all tiers enabled.
  - a list: only the named tiers are deleted. Tiers:
    - `cascade` — `bam/aln`, `bam/tagged`, `bam/charging`, `bam/classified`,
      `bam/adapter_tagged` (redundant near-copies; `bam/final` hardlinks the last one)
    - `basecall` — `bam/rebasecall` (GPU-hours to regenerate)
    - `fastq` — `fq/`, `demux/edx/fq`
    - `merged_pod5` — `pod5/` (pre-demux merged)
    - `demux_scratch` — `demux/warpdemux_output`, `demux/read_ids`, EDX read-id lists
    - `split_pod5` — `demux/pod5` (split, pre-EDX-filter)

  Always kept regardless of tiers: `bam/final`, `demux/edx/pod5` (the per-sample
  EDX-filtered POD5 used as the classification input — keeping it lets
  `classify_charging` / `classify_aa_identity` be re-run without redoing rebasecall
  or demux), plus `summary/`, `bam/aa_classified/`, `reference/`, and `logs/`.

  **Constraint:** only enable `split_pod5` for **all-EDX** runs. In non-EDX or
  mixed runs, `demux/pod5` is the classification input for non-EDX samples and must
  be kept. The on-demand `clean` rule remains the catch-all for reclaiming space on
  runs that completed with intermediates retained.

- `fasta`: Path to the reference FASTA file for BWA alignment. A BWA index will be built automatically if it doesn't exist.

- `remora_kmer_table`: Path to a table of expected normalized signal intensities for each kmer, provided by ONT at [nanoporetech/kmer_models](https://github.com/nanoporetech/kmer_models).

- `trna_table`: Path to a table with tRNA isodecoder + sequencing adapter annotation from the FASTA reference file.

  The format is four whitespace-delimited columns (no header):
  1. **uncharged tRNA name**: Name of the uncharged tRNA sequence, must match FASTA entry (e.g., `tRNA-Ala-AGC-1-1-uncharged`)
  2. **charged tRNA name**: Name of the charged tRNA sequence, must match FASTA entry (e.g., `tRNA-Ala-AGC-1-1-charged`)
  3. **isodecoder**: The isodecoder family (e.g., `Ala-AGC`)
  4. **tRNA gene name**: Representative name for the tRNA (can be any string)

  This table is currently optional since charging classification uses Remora signal analysis rather than adapter sequences.

- `opts`: Customized command-line options for pipeline tools. The `bam_filter` option controls full-length read filtering parameters. 
