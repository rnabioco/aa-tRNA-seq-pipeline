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

### YAML Format (With LDX Demultiplexing)

LDX is the successor barcode set, and `escpod demux` is the successor demux
backend. Upstream (escapepod-models) names these barcodes `nbc01`..`nbc16`;
**LDX is what we call them.** Rather than classifying boundary-gated
fingerprints, escapepod basecalls the barcode out of the raw adapter signal
with a CTC-CRF model and matches the decode to references by edit distance.

Assign barcodes with the `ldx:` key instead of `wdx:`:

```yaml
runs:
  - path: /path/to/pooled/sequencing/run
    samples:
      sample_a: { ldx: "nbc01" }
      sample_b: { ldx: "nbc02" }
      # `edx:` may still be combined with `ldx:` to filter a library down to a
      # single 3' adapter, exactly as with `wdx:`.
      sample_c: { ldx: "nbc03", edx: "edx01" }
```

and enable the backend:

```yaml
ldx:
    enabled: true
    model: "resources/models/demux/barcode_crf_nbc16_rna004@v0.2.0"
    min_margin: 0   # unclassify calls whose edit-distance margin is below this
    threads: 32
```

`warpdemux.enabled` and `ldx.enabled` are mutually exclusive — they populate
the same per-sample barcode field and their rules write the same outputs, so
turning on both is rejected at parse time rather than producing an ambiguous
DAG.

**No barcode kit is configured.** The model is a self-describing bundle
*directory* that carries its own barcode references and pins the boundary
detector it was calibrated against, so neither `--barcodes` nor `--method` is
passed. Inspect one with:

```bash
escpod demux --model resources/models/demux/barcode_crf_nbc16_rna004@v0.2.0 --info
```

Do not override the boundary detector. LLR boundaries cost 17.2 points of
balanced recall against the same classifier and the failure is silent — it runs
and produces plausible output.

**This backend is fused.** A single pass detects, basecalls, matches and routes
each read straight into its barcode's POD5, so unlike the WarpDemuX path there
is no separate read-ID extraction or `pod5 filter` split; the per-sample POD5
already exists when the command returns. The per-read classifications CSV
(`demux/read_ids/<run>/classifications.csv`) is kept because it is the only
record of each call's confidence margin.

**Performance.** The released `escpod` binary has no CUDA execution provider, so
the CRF encoder runs on CPU. Measured on 20k RNA004 reads: **59 ms of CPU per
read** for detect + encode + decode. A single 561k-read POD5 is therefore ~9
CPU-hours — about 20 minutes at 32 cores, and most of a day at 1. Size
`cpus_per_task` for `escapepod_demux` accordingly, and keep `ldx.threads` in
step with it, since that is the value `escpod` is actually launched with.

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

- `trna_table`: Path to a table with tRNA isodecoder + sequencing adapter annotation from the FASTA reference file.

  The format is four whitespace-delimited columns (no header):
  1. **uncharged tRNA name**: Name of the uncharged tRNA sequence, must match FASTA entry (e.g., `tRNA-Ala-AGC-1-1-uncharged`)
  2. **charged tRNA name**: Name of the charged tRNA sequence, must match FASTA entry (e.g., `tRNA-Ala-AGC-1-1-charged`)
  3. **isodecoder**: The isodecoder family (e.g., `Ala-AGC`)
  4. **tRNA gene name**: Representative name for the tRNA (can be any string)

  This table is currently optional since charging classification uses signal analysis rather than adapter sequences.

- `charging`: Charged vs uncharged classification, run by `escpod signal
  classify`. See `resources/models/charging/README.md` for the model bundle
  itself.

  - `model`: The model bundle **directory** (not a file). It is
    self-describing — it carries the anchor definition, the feature recipe, the
    k-mer table the features are defined against (pinned by sha256), the
    abstain rule and the recommended operating point — so no motif, offsets or
    threshold are passed as flags. Computing the features differently gives a
    wrong answer rather than an error, which is why they are not configurable.

    The bundle is vendored in this repository rather than fetched: upstream
    (`rnabioco/escapepod-models`) is private, and compute nodes have no route
    to GitHub. It also pins the escpod version — `escpod_version` must be
    >= 0.10.0 or the binary refuses the bundle outright.

  - `min_mapq`: Minimum MAPQ for a read to be classified. **0**, deliberately,
    rather than escpod's own default of 1: tRNA references are highly
    redundant, so a read mapping equally well to two isodecoders gets MAPQ 0
    from bwa and is still a perfectly good read. On the test data `--min-mapq
    1` drops 118 of 209 records (56%).

  - `ml_threshold`: `cl` at or above this is called charged. Must match the
    bundle's declared `operating_point.cl` (200, i.e. P(charged) >= 0.7824)
    unless you intend to move it. It is a recommendation measured on held-out
    data, not a property of the model — and precision depends on the *sample's*
    charged fraction, so 95% precision needs `cl >= 205` at f=0.25 but
    `cl >= 254` at f=0.05.

  **Reads the model abstains on get no `cl` tag**, rather than a default class.
  Abstention is charging-correlated (the aminoacyl adduct is what stops the
  aligner reaching the common arm), so a charging fraction over called reads
  alone is an **underestimate**. Report the no-call rate beside it:
  `summary/tables/{sample}/{sample}.charging_calls.tsv.gz` has a per-read
  `reason`, and `summary/read_attrition.tsv.gz` has the run-level breakdown.

- `escpod_version`: Version of the `escpod` binary that `pixi run setup`
  downloads. It is on the critical path of every run — POD5 merge/filter,
  charging classification, and LDX demux — and is pinned alongside the charging
  model bundle, not independently of it.

- `opts`: Customized command-line options for pipeline tools. The `bam_filter` option controls full-length read filtering parameters. 
