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
backend. The current bundle names these barcodes `ldx01`..`ldx16`. Older
`barcode_crf_nbc16_rna004` bundles emit `nbc01`..`nbc16` for the same physical
barcodes; the pipeline canonicalises those to `ldx` and records the upstream
name in an `@CO` line, so either bundle can be configured. Rather than classifying boundary-gated
fingerprints, escapepod basecalls the barcode out of the raw adapter signal
with a CTC-CRF model and matches the decode to references by edit distance.

Assign barcodes with the `ldx:` key instead of `wdx:`:

```yaml
runs:
  - path: /path/to/pooled/sequencing/run
    samples:
      sample_a: { ldx: "ldx01" }
      sample_b: { ldx: "ldx02" }
      # `edx:` may still be combined with `ldx:` to filter a library down to a
      # single 3' adapter, exactly as with `wdx:`.
      sample_c: { ldx: "ldx03", edx: "edx01" }
```

and enable the backend:

```yaml
ldx:
    enabled: true
    model: "resources/models/demux/barcode_crf_ldx16_rna004@v0.1.0"
    gpu: true            # needs the published GPU build; see Performance below
    min_margin: 0        # retired; the lattice gate supersedes it
    ref_scores: true     # record the lattice's own log P(barcode | signal)
    min_crf_margin: 1.0  # the false-positive control; swept, see below
    boundary_margin: 0   # NOT declared by any bundle; unset means escpod's 200
    clamp_max_shift: 300 # likewise. See resources/models/demux/README.md
    threads: 32
```

`warpdemux.enabled` and `ldx.enabled` are mutually exclusive — they populate
the same per-sample barcode field and their rules write the same outputs, so
turning on both is rejected at parse time rather than producing an ambiguous
DAG.

**No barcode kit is configured.** The model is a self-describing bundle
*directory* that carries its own barcode references and pins the boundary
detector it was calibrated against, so neither `--barcodes` nor `--method` is
passed. It does **not** declare `boundary.margin` or `boundary.clamp_max_shift`
— no upstream bundle does — so those two are set in config and passed as flags;
leaving them unset silently costs reads. Inspect one with:

```bash
escpod demux --model resources/models/demux/barcode_crf_ldx16_rna004@v0.1.0 --info
```

Do not override the boundary detector. LLR boundaries cost 17.2 points of
balanced recall against the same classifier and the failure is silent — it runs
and produces plausible output.

**This backend writes no POD5.** A single pass detects, basecalls and matches
each read, and `--annotate` records the assignment in a `.p5s` sidecar written
next to the POD5 it describes — in the raw data directory, not under
`output_directory`. Nothing is copied: where the WarpDemuX path leaves a second
full set of the run's reads on disk, an LDX run leaves a few MB.

The split happens later and one level up, on the basecalled reads. The run is
basecalled whole (`rebasecall_ldx_run`), restricted with dorado's `-l` to the
reads demux actually assigned to a sample, and then cut into per-sample uBAMs
(`split_ldx_ubam`) at `bam/rebasecall/<sample>/`. That is the same path the
WarpDemuX path's `rebasecall` produces, so everything downstream is identical.

The per-read classifications CSV (`demux/read_ids/<run>/classifications.csv`) is
kept for both of its jobs: it is the read→barcode source those two rules read,
and it is the only record of how each call went — `confidence` always, and,
under `ref_scores` (on by default), `crf_logp`, `crf_margin`, `crf_best` and
`mean_logpost`. The sidecar carries those same scores from escpod 0.12.0, which
is when `.p5s` columns became numeric.

Consequences worth knowing:

- **No per-sample POD5 exists.** The signal consumers — `classify_charging` and
  the signal-metrics QC — are pointed at the raw run instead; they walk the
  sample's BAM and look each read's signal up by id, so reads outside the sample
  are never touched. To get one anyway — for squiggy, or to inspect a barcode by
  hand — cut it on demand:

  ```bash
  escpod filter <run>/pod5 --annotation barcode=ldx05 -o ldx05.pod5
  ```

- **Re-demuxing is cheap and safe.** Running demux again over the same POD5
  replaces the sidecar's `barcode` column in place, so changing the model or a
  gate costs one pass and no cleanup.

- **A sidecar is bound to its POD5.** escpod checks the POD5's footer UUID and
  size before reading one, so a POD5 replaced under the same name makes the
  sidecar fail loudly rather than silently describe reads that are no longer
  there. Recover by deleting the `.p5s` and re-running.

- **One sidecar per POD5, one `barcode` column in it.** Two configs demuxing the
  same run with different models would overwrite each other's assignments. If
  you need to compare models on one run, copy the POD5 or compare the
  classifications CSVs instead.

**Performance.** The released `escpod` binary has no CUDA execution provider, so
the CRF encoder runs on CPU. Measured on 20k RNA004 reads: **59 ms of CPU per
read** for detect + encode + decode. A single 561k-read POD5 is therefore ~9
CPU-hours — about 20 minutes at 32 cores, and most of a day at 1. Size
`cpus_per_task` for `escapepod_demux` accordingly, and keep `ldx.threads` in
step with it, since that is the value `escpod` is actually launched with.

`ref_scores: true` adds **+3.6%** to that, which is why it is on by default. It
costs more under `ldx.gpu: true`: the constrained scan reads the raw scores, so
the decode comes back to the host while the encoder stays on the device.

**GPU is worth it, and by more than upstream's headline figure.** Measured here
on the 2026-08-06 flowcell, `--ref-scores` on, one GPU + 8 cores against the
released CPU binary on 32 cores:

```
                          wall        reads/s   yield
CPU (released, 32 cores)  ~2.6 h          108   92.22%
GPU (0.12.0-gpu, 1 GPU)   458.9 s       2,182   92.22%
```

**20.3x**, at identical yield and using a quarter of the cores. Upstream quotes
2.4x because that compares the GPU encoder against the CPU encoder *within the
same GPU-capable binary*; this compares what you would actually switch between.
The stage trace shows the device saturated (workers busy 94% of wall, producer
blocked on send 347.7 s), so a second GPU would push it further.

**`ldx.gpu` defaults to `true`**, which needs the GPU build of escpod. Since
escapepod-rs 0.17.1 that is a published artifact rather than a source build, so
`pixi run setup` installs it along with everything else:

```bash
pixi run setup             # musl default + <version>-gpu, both checksum-pinned
pixi run install-ort-gpu   # CUDA onnxruntime (not on conda-forge; see pixi.toml)
pixi install -e gpu        # cuDNN
```

The GPU artifact is x86_64 Linux only, and is the one **dynamically linked**
build (glibc >= 2.28) because the CUDA runtimes are dlopened — so `ldx.gpu:
true` is simply unavailable on macOS and aarch64. `pixi run setup` skips it
there rather than failing.

`escapepod_demux` resolves that `<escpod_version>-gpu` binary itself and fails
loudly if it is missing. **`escpod_version` is shared with the musl default** —
demux is the only rule with a GPU path, so pointing the global pin at `-gpu`
would hand every other rule a dynamically linked, single-platform binary just
to run `escpod merge` and `escpod signal classify`, and would break `pixi run
setup`, which derives its download URL from that string.

`scripts/install-escpod-gpu.sh` still exists for building an unreleased ref out
of the private repo, but is no longer on anyone's normal path.

Without cuDNN the CUDA provider fails to register and onnxruntime falls back to
CPU with **only a warning**, so confirm the log says
`CRF encoder: N worker(s) on GPU [0]`.

Set `ldx.gpu: false` on a CPU-only host — the test configs do exactly that,
since `ldx.gpu: true` is resolved while the DAG is built and CI has no GPU.

### Testing the LDX path

`config-ldx-test.yml` runs the whole demux path end to end against a committed
fixture — 415 reads of a real pooled barcoded run, in
`.tests/fixtures/ldx-demux`:

```bash
pixi run dry-run-ldx    # DAG only; no GPU, no download
pixi run test-ldx       # full run (dorado needs a GPU)
```

The fixture holds three claimed barcodes plus two populations that must *not*
become samples — a barcode no sample claims, and reads the CRF could not call —
and each filtered barcode deliberately carries wrong-adapter reads so
`filter_{fastq,pod5}_by_edx` has something real to remove. See the README beside
it for the full composition and how to rebuild it.

Two things about this config are easy to get wrong when copying it:

- **The `adapters.three_prime` override is load-bearing.** `config-base.yml`
  also defines names `edx01`/`edx02`, but with the sacCer3 dual-adapter
  *sequences* — same names, different molecules. Inherit those and every read
  detects as `none`, which is exactly the failure mode of issue #120.
- **All seven adapters are listed**, though only two are filtered on. Detection
  can only assign a read to an adapter in the list, so omitting the rest would
  collapse the off-target reads to `none` and make the concordance table
  meaningless.

### `min_margin` is retired to 0, and this is why it was 12

Measured on a 1,001,307-read run against a second, independently trained bundle
(two models disagreeing is a floor on error):

```
threshold  drops      of calls   of which disagreements
     1     1,295       0.14%          96.2%
    12     8,817       0.95%          85.2%   <- the default
    13   431,747      46.76%           8.5%   <- cliff, do not
```

The residual error rate barely moves, which makes 12 look pointless — but what it
*discards* is 85% wrong calls, removing 10.5% of all disagreements for 0.95% of
reads. That is the trade to make when a misassignment means cross-sample
contamination and an unclassified read only costs yield.

13 is a trap: 99.0% of reads sit on a margin plateau of 12/13/14 (the references
are >=12 apart by design, so a wrong decode scores like a right one), so 13 eats
half the run at 8.5% precision.

**Set 0 when the measurement IS the crosstalk** — an adapter-ligation QC run
exists to count misassignment, and gating removes exactly the reads it is
counting.

### `crf_margin` is the false-positive control, and 1.0 is the operating point

`min_margin` can only ever reach ~1% of reads, because it measures the wrong
thing. `confidence` is the edit-distance margin to the runner-up, and on a
*designed* panel that measures how far apart the references are, not how sure
the model is — 90% of the reads two independently trained bundles disagree
about are **exact** matches to a reference, so no edit-distance threshold can
see them.

Since escpod 0.12.0, `ref_scores: true` (the default) asks the lattice instead.
Restricting the CRF forward recursion to the paths that emit a given reference
and normalising by the full partition function is a real probability:

```text
crf_logp = logZ_target(reference) - logZ_full = log P(reference | signal)
```

`classifications.csv` gains four columns, appended after `confidence` (so
anything parsing it positionally, `summarize_demux.awk` included, is unaffected):

| column | meaning |
|---|---|
| `crf_logp` | log-probability of the **called** barcode |
| `crf_best` | the reference the lattice itself prefers — need not be the one edit distance called |
| `crf_margin` | log-odds in nats against the runner-up; **negative** when `crf_best` disagrees with the call |
| `mean_logpost` | the decoded path's mean per-timestep log-posterior |

The resolution difference is the whole point: over 20k RNA004 reads
`confidence` takes 15 distinct values with 98.7% of reads in three of them,
while `crf_margin` takes 14,818 over 16,747 reads. Even within the 98.4% of
reads that match a reference exactly, `P(barcode | signal)` still spans from
below 0.1 up to 0.9–0.99 — reads a clean decode cannot tell apart and the
lattice can.

**Swept on the 2026-08-06 donor flowcell** (1,001,307 reads, 923,431 calls at
margin 0), scoring precision against a second, independently trained bundle —
the same floor-on-error methodology as the `min_margin` table above:

```
gate                 recall%   of all errors removed   discard precision%
min_margin 12         91.34            9.52                  83.78
min_crf_margin 0.5    89.53           29.87                  86.12
min_crf_margin 0.7    88.74           38.59                  85.84
min_crf_margin 1.0    87.85           47.69                  84.56   <- set
min_crf_margin 2.3    85.59           64.48                  75.38
min_crf_margin 4.6    80.18           77.71                  50.02
```

**1.0 is where the dial stops being free.** Its discards are 84.56% real errors —
the same quality as `min_margin 12`'s 83.78% — while removing 47.69% of all
error against that gate's 9.52%. Five times the error caught at the same discard
precision, for 3.5 points of recall. Past ~1.3 discard precision decays; by 4.6
half of what is dropped is correct.

**`min_margin` goes to 0 because it is not complementary.** On top of the crf
gate it removes only 728 further reads, **90.5% of which are correct**, moving
total error removed by +0.09 points. Since the crf gate's own discards are ~85%
real errors, keeping 12 lowers the average quality of what the pipeline throws
away. Raise it again only for a run with no lattice gate.

`min_crf_prob` is left unset — a different cut (absolute confidence rather than
separability), and not the one that was swept.

These margins were measured on the **GPU encoder**, unlike the CPU-measured
`min_margin` table. Calls agree with the CPU encoder on 99.76% of fixture reads
and `crf_margin` differs at a median of 0.0016 nats, so the threshold is not
sensitive to that — but the two tables are not strictly one measurement.

`config-demux-test.yml` (WarpDemuX) is a **dry-run target only** and cannot
complete a run; its header explains why, and there is no committed WDX fixture.

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
    - `basecall` — `bam/rebasecall`, `bam/rebasecall_run` (GPU-hours to regenerate)
    - `fastq` — `fq/`, `demux/edx/fq`
    - `merged_pod5` — `pod5/` (pre-demux merged per-sample, or per-run on an LDX
      run whose reads span several POD5 directories)
    - `demux_scratch` — `demux/warpdemux_output`, `demux/read_ids`, EDX read-id lists
    - `split_pod5` — `demux/pod5` (WarpDemuX split, pre-EDX-filter)

  Always kept regardless of tiers: `bam/final`, `demux/edx/pod5` (the per-sample
  EDX-filtered POD5 used as the classification input — keeping it lets
  `classify_charging` be re-run without redoing rebasecall or demux), plus
  `summary/`, `reference/`, and `logs/`.

  **Constraint (WarpDemuX only):** only enable `split_pod5` for **all-EDX** runs.
  In non-EDX or mixed runs, `demux/pod5` is the classification input for non-EDX
  samples and must be kept. LDX runs produce no `demux/pod5` at all — their
  classification input is the raw POD5 plus the sidecar — so the tier is inert
  there. The on-demand `clean` rule remains the catch-all for reclaiming space on
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
