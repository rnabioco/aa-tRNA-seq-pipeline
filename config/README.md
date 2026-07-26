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

### YAML Format (With Barcode Demultiplexing)

For multiplexed/pooled sequencing runs using WDX barcodes, use a YAML file:

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

### Enabling Demultiplexing

To use demultiplexing, add the following to your config file:

```yaml
samples: config/samples-demux.yml  # YAML format sample file

warpdemux:
    enabled: true
    backend: "escpod"                     # "escpod" (default) or "warpdemux"
    barcode_kit: "WDX4_tRNA_rna004_v1_0"  # warpdemux backend: default kit if not set per-run
    save_boundaries: true  # warpdemux backend: saves adapter boundary information
    threads: 8

    # escpod backend only
    method: "cnn"  # "cnn" (default, matches how the barcode model was trained) or "llr"
    barcode_model: "resources/models/demux/barcode_wdx4_rna004.gbm.json"
    adapter_model: "resources/models/demux/adapter_rna004.onnx"  # needed by method: cnn
```

The config section is still called `warpdemux` for backward compatibility, but
`backend` now selects between two implementations that write the same downstream
files (`demux/read_ids/{run_id}/barcode_mapping.tsv.gz`, `demux/pod5/{sample}/{sample}.pod5`):

- `escpod` (default) — one fused `escpod demux` pass that classifies barcodes *and*
  writes one POD5 per barcode, so no separate read-ID/filter pass is needed.
- `warpdemux` — the original python implementation, followed by an `escpod filter` pass.

**The two do not produce identical barcode calls.** `escpod` cannot load a WarpDemuX
kit, so it uses `barcode_wdx4_rna004`, a GBM distilled from the
`WDX4_tRNA_rna004_v1_0` teacher. Barcode numbering is unchanged, but the shipped GBM
has no per-class confidence thresholds, so every read with a usable adapter boundary
gets a barcode and the `unclassified` fraction drops sharply. Install the escpod
models with `pixi run install-demux-models`. See
[docs/workflow/demultiplexing.md](../docs/workflow/demultiplexing.md) for the full
comparison, caveats, and how to validate a backend switch.

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
    - `demux_scratch` — `demux/warpdemux_output`, the per-run `barcode_mapping.tsv.gz`
      and per-sample read-id lists under `demux/read_ids`, the escpod
      `classifications.csv`, EDX read-id lists
    - `split_pod5` — `demux/pod5` and `demux/escpod_output` (split, pre-EDX-filter)

  Always kept regardless of tiers: `bam/final`, `demux/edx/pod5` (the per-sample
  EDX-filtered POD5 used as the classification input — keeping it lets
  `classify_charging` / `classify_aa_identity` be re-run without redoing rebasecall
  or demux), `demux/read_ids/{run_id}/demux_summary.tsv.gz` (the durable record of
  what the demultiplexer called), plus `summary/`, `bam/aa_classified/`,
  `reference/`, and `logs/`.

  **Constraint:** only enable `split_pod5` for **all-EDX** runs. In non-EDX or
  mixed runs, `demux/pod5` is the classification input for non-EDX samples and must
  be kept. The on-demand `clean` rule remains the catch-all for reclaiming space on
  runs that completed with intermediates retained.

- `fasta`: Path to the reference FASTA file for BWA alignment. A BWA index will be built automatically if it doesn't exist.

- `remora_kmer_table`: Path to a table of expected normalized signal intensities for each kmer, provided by ONT at [nanoporetech/kmer_models](https://github.com/nanoporetech/kmer_models). Optional (default `null`); when set, the `remora_signal_stats` QC path runs, which is the only place the pipeline still uses the ONT `pod5` python package.

- `classifier`: Charging classifier, `remora` (default) or `leech`. `leech` is a
  GPU-accelerated alternative installed by `pixi run setup` / `pixi run install-leech`.
  Note leech writes the `ML` and `MP` tags where remora writes `ML` and `MM`.

- Pinned tool versions, read by the setup scripts as well as the workflow:
  - `dorado_version` / `dorado_model`
  - `escapepod_version`: the `escpod` CLI tag built from source *and* the `escapepod`
    PyPI package version. escpod handles all POD5 merge/filter in the pipeline and
    must be a source build — the published escapepod-rs release binaries are
    default-features only and do not contain a working `escpod demux`. Building
    requires a Rust toolchain >= 1.95.
  - `leech_version` / `leech_core_version`: the leech release wheel and the optional
    `leech-core` Rust extension wheel. Fetched from the private `rnabioco/leech`
    release with an authenticated `gh`, or from `LEECH_WHEEL_DIR`. Without
    `leech-core`, leech falls back to a slower pure-python extraction backend.

- `trna_table`: Path to a table with tRNA isodecoder + sequencing adapter annotation from the FASTA reference file.

  The format is four whitespace-delimited columns (no header):
  1. **uncharged tRNA name**: Name of the uncharged tRNA sequence, must match FASTA entry (e.g., `tRNA-Ala-AGC-1-1-uncharged`)
  2. **charged tRNA name**: Name of the charged tRNA sequence, must match FASTA entry (e.g., `tRNA-Ala-AGC-1-1-charged`)
  3. **isodecoder**: The isodecoder family (e.g., `Ala-AGC`)
  4. **tRNA gene name**: Representative name for the tRNA (can be any string)

  This table is currently optional since charging classification uses Remora signal analysis rather than adapter sequences.

- `opts`: Customized command-line options for pipeline tools. The `bam_filter` option controls full-length read filtering parameters. 
