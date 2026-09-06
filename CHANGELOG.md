# Changelog

All notable changes to the aa-tRNA-seq pipeline are documented in this file.

## [Unreleased]

### Changed

- **`verify-demux-model.py --record` merges into the manifest instead of
  replacing it.** It rewrote the file from the bundles named on the command
  line, so recording one bundle by path dropped the digests of every other
  bundle in that directory -- and since a missing digest fails closed, the next
  `verify-demux-model` refused all of them. That is the failure the
  tracked/local split was built to prevent, arriving through a different door:
  a deployment carrying several hand-vendored bundles would lose the lot the
  moment it recorded a new one by path. Entries whose bundle is no longer on
  disk are still dropped, so the manifest keeps describing what is there.

- **`escpod` bumped 0.20.0 -> 0.21.0.** No demux call moves: no signal
  orientation change, no `.p5s` format change (still `p5s_version` 3), no GPU
  behaviour change, and the CPU boundary CNN emits identical sorted
  classifications.

  What it unblocks is `fdx.fused: true`. escapepod-rs#323 scopes
  `--boundary-margin` / `--clamp-max-shift` to the heads that actually have a
  boundary detector: read-end heads ignore them and are named as doing so, and
  only a run of read-end models ALONE is refused. The blocker this repo
  documented under `fdx:` -- "Flip to `true` once the pinned escpod accepts the
  flags in a fused run" -- is gone. The default stays `false` until the fused
  path is measured here; the win when it is taken is collapsing the FDX second
  pass into one POD5 sweep, hours of IO-bound wall on a large flow cell.

  Two diagnostics also arrive. The classifications CSV carries **`adapter_end`**
  as a trailing column wherever a boundary detector ran, so calls can be split
  by `adapter_end` band against an independent label per run -- present on the
  LDX pass, absent on the FDX pass, so the two axes' CSVs now differ in column
  count. Both consumers read by name (`summarize_demux.awk` resolves by header
  with a `c = 2` positional fallback a trailing column cannot move;
  `select_demux_reads.py` uses `csv.DictReader`), so nothing breaks. And
  **`Demux summary:` breaks `unclassified` down by reason** per axis -- boundary
  margin, window clamp, the `adapter_end = 0` sentinel, short read, encoder
  error, no-match decode, and the gates, only the last of which is a classifier
  decision. Nothing parses that output; it reaches logs only.

  Also corrects the documented escpod floor in `docs/`, which still said
  `>= 0.10.0` where the real floor is 0.19.0, and the stale `"0.12.0"` fallback
  in `workflow/Snakefile` (inert, since config-base.yml always sets the key).

- **`verify-demux-model.py` no longer refuses a bundle that declares
  `boundary.margin` / `boundary.clamp_max_shift`.** It refused any bundle that
  did, on the reasoning that no upstream CRF bundle ever had -- true until
  escapepod-models#127/#128 (2026-09-06), which reissued
  `barcode_crf_wdx4_rna004@v0.2.1` and `barcode_crf_ldx32_rna004@v0.2.2`
  declaring both keys deliberately, stating escpod's own fallbacks (200 / 0) so
  a run cannot silently differ from the values the bundle was characterised at.
  The check had become unsatisfiable: `verify-demux-model` failed on the current
  `latest` of two families while `check-models` reported those same versions as
  BEHIND.

  The refusal was compensating for a real gap, so it is replaced rather than
  dropped. A bundle's internal hashes cover the ONNX files only -- `metadata.json`
  is *where those hashes live*, so nothing upstream can attest it without being
  self-referential, and an edit to it (adding window rules, as #107/#108 did to
  nbc16) passes every other check. Both sidecars are now checked against digests
  recorded at vendoring time in `resources/models/demux/sidecars.sha256`, and a
  bundle with no recorded digest **fails**, rather than skipping. That is
  tamper-evidence rather than authenticity -- it proves a bundle has not changed
  since it was vendored, not that it is upstream's -- but it is what the old
  check was approximating, it is version-agnostic, and it stays offline, which
  matters because compute nodes have no route to GitHub.

  A deployment that vendors bundles this repo does not ship records them in
  `sidecars.local.sha256` (gitignored) via `--record --local`, so pulling an
  updated tracked manifest cannot drop their digests and fail them closed.

- **`barcode_crf_wdx4_rna004` vendored at v0.2.1**, the reissue described above.
  Weights are byte-identical to v0.2.0 (`a33459d1...` CRF, `b59f8667...`
  boundary) and every published metric is unchanged, so no demux call moves.

- **`get_bcerror_freqs.py` reads each reference sequence once** instead of
  calling `faidx.fetch()` for every aligned base. The per-base fetch was inside
  the CIGAR loop, so it ran once per base of every read in the sample. On a real
  flow-cell sample (46,577 reads, the 2026-09-05 GlnRS pilot) this takes the
  rule from 12.4 s to 3.7 s, a 3.3x speedup, with byte-identical output. The
  committed fixtures are too small to show it (1.3-1.6x at 200-400 reads);
  per-reference setup dominates there.

- **`get_mismatch_calls.py` reads each reference sequence once**, the same
  per-base `faidx.fetch()` in the same shape, in the script that is the
  per-read counterpart of `get_bcerror_freqs.py`. On the same flow-cell sample:
  11.5 s to 6.5 s, 1.8x, output identical. These two were the only per-base
  reference fetches left in `workflow/scripts/`.

- **`get_align_stats.py` no longer holds every query name in a set.** Reads are
  deduplicated by skipping secondary and supplementary records, which is what
  "count each read once" means, rather than by remembering names.

  This is a correctness fix as much as a memory one. The set counted whichever
  record for a query name came FIRST, so on a BAM carrying supplementary
  alignments the read was measured by whichever alignment the file happened to
  list first -- its length, base quality and MAPQ could all come from a partial
  supplementary record. It now always measures the primary.

  Read counts do not move for BAMs this pipeline writes today: `bwa_align` keeps
  primary forward alignments only (`-F 2324`), and a v0.7.2 flow-cell sample
  measures 46,577 records against 46,577 distinct query names, with no secondary
  and no supplementary records. **Re-running `align_stats` over a tree written
  before that filter existed will move its `aligned` row** -- same read count,
  but length, quality and MAPQ now taken from the primary alignment.

  The set cost a measured 122 bytes per read: 5.7 MB on a 46.6k-read sample, but
  ~1.2 GB on a 9.6M-read whole-run BAM against the rule's 8 GB default.


## [v0.7.3] - 2026-09-06

### Added

- **`charging.orientation`**, forcing the move-table frame for samples too thin
  for `escpod classify` to detect it. `auto` (the default) needs >= 50
  informative reads and a 95% consensus, and below that floor it is a hard
  ERROR, not a warning -- so on a heavily multiplexed run a sparse condition
  (~20 anchored reads at 224 samples per flow cell) fails, and since every
  sample owes a charging table, one such sample fails `rule all` for the whole
  corpus. `auto` emits no flag at all, so an existing run's command line is
  byte-identical to before; `time` / `reversed` force the frame, and anything
  else stops the DAG rather than 448 individual jobs.

  The frame is a property of the RUN -- chemistry plus basecaller -- not of a
  sample, so the value must be re-derived per run rather than carried between
  projects. Every sample deep enough for `auto` logs what it found, so a run
  can be asked directly:

  ```
  grep -h "move-table frame" <outdir>/logs/classify_charging/* | sort | uniq -c
  ```

  In practice it is unanimous: on the 2026-09-05 GlnRS flow cells, 446 samples
  and 8,057,646 votes returned `reversed` with not one dissenter. Forcing it for
  the few samples below the floor inherits that measurement rather than guessing.
  A forced frame that is WRONG does not error -- it mis-anchors every feature and
  the calls silently become noise -- so if a run's deep samples disagree with
  each other, do not set it.

  Note what it buys: completeness, not signal. A sample that fails `auto` has
  tens of reads and its charging fraction carries a Wilson interval too wide to
  use. It makes the run finish; keep excluding those samples on a depth guard.

- **`charging.orientation_fallback`** (default `reversed`), which supplies that
  frame automatically for a sample `auto` could not decide, instead of failing
  it. This is the part that matters by default: shipping `auto` alone meant
  every high-multiplex run met the failure, diagnosed it and opted out, which
  is a footgun rather than a safeguard.

  It is a fallback rather than a forced default because detection is worth
  keeping wherever it works. The frame's two inputs are already pinned and
  enforced -- `base_calling_model` is `rna004_sup@v6.0.0` and
  `charging.basecaller_check` is a hard error on a mismatch -- so the frame is
  effectively constant, and forcing it everywhere would stop escpod computing
  the per-sample consensus that would notice if that ever stopped being true.
  Deep samples detect; only the thin ones inherit.

  `escpod classify` now runs through `workflow/scripts/escpod_classify_fallback.sh`,
  which retries ONLY on `orientation check underpowered`. Any other failure --
  OOM, a missing POD5, an unreadable model -- propagates unchanged and is never
  retried, so a real problem is never masked. Every supplied frame is announced
  in the sample's log, so a run records which samples had one rather than
  detecting their own. `orientation_fallback: none` restores the v0.7.2
  failure.

  The measurement behind `reversed` is specific to `rna004_sup@v6.0.0`; the
  other vendored bundle has not been measured, and the config says so.

## [v0.7.2] - 2026-09-05

### Fixed

- **One LDX code fanned across several EDX adapters no longer dies mid-run.**
  v0.7.1 fixed the parse-time uniqueness check (#161) but left two more copies
  of the same EDX-blind reasoning, and the run got further before failing --
  past `escapepod_demux`, so the GPU pass was paid for first (#163).

  `select_demux_reads.py` refused the design outright. Its `--sample` specs
  carry only the SIGNAL axes (`ldx`, `fdx`), because `edx` is read off the uBAM
  by `detect_edx_adapters` and has no classifications CSV at that point, so
  samples differing only by adapter arrived as identical tuples. `edx` is now
  passed as a `--downstream-axis`: it counts towards the uniqueness checks and
  is explicitly skipped when selecting reads. That second half is load-bearing
  -- a downstream axis has no calls, so matching on it would compare every read
  against `None` and hand back empty samples rather than an error.

  The all-or-none rule that exists for `fdx` was also missing for `edx`, in both
  `_check_barcode_tuples` and `run_config.py`. A sample naming the LDX code
  alone beside one naming `ldx01` + `edx01` is never EDX-filtered, so it
  swallows the other's reads. Both now refuse it, and `is_edx_enabled()` joins
  its `is_ldx_enabled` / `is_fdx_enabled` siblings.

  Verified against the run that hit this: 224 samples and a real 387 MB
  `classifications.csv`, EXIT=0 and 4,200,833 read ids, per-sample counts
  non-zero.

## [v0.7.1] - 2026-09-05

### Fixed

- **An LDX code repeated across EDX adapters no longer dies at parse time.**
  `_check_barcode_tuples` (new in v0.7.0, #156) built its uniqueness tuple from
  `(ldx, fdx)` and left `edx` out, so any samples file fanning one LDX code
  across several 3' adapters was refused with "have identical barcode
  assignments" -- a design the pipeline supports and its own validator accepts,
  since the sibling check in `run_config.py` keys on `(fdx, edx)`. The two
  checks disagreed, so `pixi run check-run-config` passed configs the pipeline
  then rejected. `edx` is now in the tuple. `tests/unit/test_barcode_tuples.py`
  covers both axes and asserts the two checks return the same verdict; no
  shipped test config repeated an LDX code, which is why CI missed it (#161).

## [v0.7.0] - 2026-09-05

### Added

- **Pins are checked against upstream on a schedule.** `scripts/check-currency.py`
  (was `check-model-currency.py`; `pixi run check-currency`, with
  `check-models` and `check-escpod` for either half) now also reports
  `escpod_version` and `dorado_version` against the escapepod-rs and dorado
  releases, and compares the escpod checksums in `scripts/setup-tools.sh` with
  the pinned release's `SHA256SUMS.txt` -- a pin and its checksums have to move
  together, and that is what a hand bump forgets. `.github/workflows/currency.yml`
  runs it every Monday (plus `verify-demux-model` / `verify-charging-model`);
  the private escapepod-models half needs the repository secret
  `ESCAPEPOD_MODELS_TOKEN`, without which that step fails visibly rather than
  passing silently. Holds for tools live in `resources/models/pins.yml` under
  `tools:`; dorado 2.1.1 is recorded there, held to the basecaller the charging
  bundle declares while upstream is at 2.1.2.
- **The FDX 5' index, as a second demux axis.** `fdx.enabled` (with
  `ldx.enabled`) calls the 5' barcode of a dual-index library with the vendored
  `barcode_crf_fdx4_rna004@v0.2.0` bundle — four codes `fdx01`..`fdx04`, anchored
  on the read end, Run2 of the pilot held out of its training — and a sample is
  now `{ldx: ldx01, fdx: fdx01}`, owning the reads on which both calls agree.
  New `select_demux_reads.py` does that join for the run's basecall list
  (`ldx_run_read_ids`) and per sample, replacing the awk in `rebasecall_ldx_run`;
  it reads escpod's two-pass CSVs and its fused per-axis CSV alike, applies each
  axis's lattice gate itself, and refuses sample tuples that would overlap.
  Final BAMs of dual-index samples carry `BC:Z:ldx01-fdx01`.

  escpod can call both models in one sweep and `fdx.fused: true` asks for it,
  but escpod 0.19.0 refuses `--boundary-margin`/`--clamp-max-shift` whenever a
  read-end model is in the run and the LDX axis cannot give those up, so the
  default is a second pass (`escapepod_demux_fdx`) that annotates the same
  sidecar (`fdx` column beside `barcode`). The bundle's declared 3.5-nat gate
  ships as `fdx.min_crf_margin`; measured on the LDX fixture, whose reads carry
  no 5' index, it cuts confident mis-assignments from 17% to 2.2%.

  A second committed fixture, `.tests/fixtures/fdx-demux` (365 reads of Run2,
  E. coli, plain 3' adapter), with `config-fdx-test.yml`, `pixi run dry-run-fdx`
  (in CI) / `test-fdx`, and `tests/integration/test_fdx_demux.py`.
  `summarize_demux.awk` takes `-v col=` to tally an axis of a fused CSV.
- **`pixi run new-run-config` / `check-run-config`** (`workflow/scripts/
  run_config.py`): write a run's `config-<name>.yml` + samples file from
  `--run`/`--sample NAME=SPEC` (`ldx01`, `ldx01+fdx01`, `ldx01/edx01`,
  `wdx:barcode03`, `-`)/`--reference-raw`/`--three-prime`, or validate an
  existing pair, refusing what the pipeline would refuse: codes the vendored
  bundles do not emit, fdx without ldx, overlapping tuples, undeclared edx
  names, 3' adapters of unequal length or not starting with GGC, run dirs
  without POD5, missing references. `--json` for callers that are programs;
  `--dry-run` builds the DAG. `.claude/skills/new-run-config/SKILL.md` is the
  agent-facing walkthrough.

### Changed

- **`ldx.clamp_max_shift` is off and `ldx.boundary_margin` is back at 200.**
  The `0` / `300` set in v0.2.0 (#123) rested on "median edit distance 0, 98% within 2
  edits" for the recovered reads — a metric a window holding no barcode passes,
  because such a window still decodes to a *perfect* reference. Scored instead
  against the independent 5′ index on the same molecule (ldx32 v0.2.1, FDX
  Run1, rnabioco/escapepod-rs#323), the band the clamp reaches is 28%
  design-consistent with 45% of its calls on codes absent from the pool, and
  the margin band is 55% (79% under this pipeline's lattice gate of 1.0, at 63%
  of the band) against 85% for full-window reads. The full table and the
  method are in `resources/models/demux/README.md`; a run config that gates at
  `min_crf_margin: 2.0` may set `boundary_margin: 0` for ~1% more reads at 84%.
  Expect demux yield to drop by roughly 3.8 points on an RNA004 tRNA run, and
  the reads that leave to have been mostly wrong.
- **`escpod` bumped 0.19.0 -> 0.20.0.** Nothing the pipeline calls changed:
  every flag the rules pass exists with the same meaning, the sidecar and CSV
  shapes are 0.19.0's, and the vendored bundles load as before. What it brings
  is performance and hardening -- the CPU boundary CNN is ~1.6x faster (input
  shape pinned, padding hoisted; per-read classifications identical), per-read
  work on the LLR path is bounded so one stalled-pore read cannot hold a worker
  for seconds, a failed GPU `--ref-scores` upload is reported once rather than
  retried every batch, and a reads-table row that fails to decode is warned
  about rather than dropped. `pixi run setup` installs it; the checksums in
  `scripts/setup-tools.sh` are the release's. Measured, not assumed: the fused
  `--model ldx=… --model fdx=…` pass still refuses `--boundary-margin` with a
  read-end model in the run, so `fdx.fused` stays off by default.
- **`read_attrition` names the join's loss.** On the escpod path the
  `barcode assigned -> basecalled` row is now two: `barcode assigned -> sample
  assigned` (a code no sample claims, or a dual-index pair the axes did not
  agree on, read from the new per-run `assigned_summary.tsv`) and `sample
  assigned -> basecalled` (what the basecaller could not read). The two were
  one number before, labelled as basecaller loss -- which is how the LDX
  fixture's 40 decoy reads had always appeared.

- **Alignment carries dorado's tags through bwa; `inject_ubam_tags` and
  `ubam_to_fastq` are gone.** `bwa_align` now streams the uBAM through
  `samtools fastq -T '*'` into `bwa mem -C`, so the move table (`mv`/`ns`/`ts`),
  the MM/ML modbase calls and `RG` land on the aligned records directly. That
  removes one full BAM copy per sample (`bam/tagged`), the on-disk FASTQ, and a
  Python job budgeted at 48 GB and 6 h that re-read the uBAM to put the same
  tags back. The dangling-@RG repair (#121) moved with it: a small
  `stamp_read_groups.py` emits the uBAM's `@RG` lines with SM/LB/BC stamped for
  `bwa mem -H`, kept beside the BAM as `{sample}.rg.sam`; the per-read `BC:Z:`
  on demux runs is written into the FASTQ comment. Output is primary forward
  alignments (`-F 2324`), which is what the tagged BAM used to contain.

  This is the `-C` that v0.1.1 (#86) removed over an OOM at 48 GB. The tag
  payload is real — the FASTQ comment is ~13x the read, the move table dominates
  — but bwa only holds one input batch of it, and the rule now pins that batch
  with `-K 100000000` (bwa's default is 10 Mbases *per thread*), which bounds
  the comment text at a few GB. The 160 GB this rule is budgeted for today is
  set by `-k 6` on a dense reference and dwarfs it. Re-measure MaxRSS on the
  first full flowcell through this path; the profile comment says so too.
- **No POD5 is copied any more.** `merge_pods` — a full second copy of every
  run per sample, kept for the life of the analysis — is replaced by
  `stage_pod5`, which lays a directory of symlinks over the raw files
  (`pod5/{sample}/<run>/<pod5_pass|pod5_fail|pod5>/`). dorado basecalls it
  with `--recursive` and `escpod classify` walks it, exactly as the LDX path
  has always handed both tools the raw run. `filter_pod5_by_edx` (an
  EDX-filtered POD5 per sample, kept forever) and `filter_fastq_by_edx` are
  gone for the same reason: the classifier only touches reads the BAM names,
  and `bwa_align` bounds the alignment with `samtools view -N` on the uBAM.
  The `fastq` and `merged_pod5` cleanup tiers are retired (accepted and
  ignored; an unknown tier name is now an error), and `split_pod5` now
  describes the WarpDemuX split POD5 as that path's classification input.
- **The Slurm profile no longer clamps every rule to one CPU.** The plugin uses
  `cpus_per_task` when the resource is set and falls back to `threads:` only
  when it is not, so the `cpus_per_task: 1` in `default-resources` overrode
  every rule without its own `set-resources` entry — `bamCoverage -p 4` and
  `compute_seq_similarity --threads 4` ran on a single core. Removed; the
  per-rule entries are unchanged.
- The Squiggy session lists each sample's real POD5 files (new
  `generate_squiggy_session.py --pod5 SAMPLE=PATH`). Before, it wrote
  `pod5/{sample}/{sample}.pod5` for every sample, which on an LDX run never
  existed.

### Removed

- `transfer_tags.py` and its tests, superseded by `stamp_read_groups.py`.

## [v0.6.0] - 2026-09-04

### Added

- **An end-to-end test for basecall resume** — `pixi run test-ldx-resume`
  (`.tests/run_ldx_resume_e2e.sh`). Needs a GPU, so it is not in CI; it is the
  only thing that exercises `dorado_basecall_resume.sh` against real dorado,
  and it rebuilds `.tests/outputs-ldx`.

  It kills a basecall the way Slurm does — SIGTERM to the whole process group —
  and asserts Snakemake's failure path takes `{output}` while the checkpoint and
  its `.args` survive, which is the blocker in #149. Because the fixture's 323
  reads basecall in far less time than dorado takes to start, no signal can land
  mid-write, so the resume itself is driven from a checkpoint holding real dorado
  records: the baseline uBAM truncated mid-record, carrying the wrapper's own
  `.args`. It then checks the resumed run against a clean one — read ids and
  sequences, no duplicates, MM/ML/`mv` tag counts (those records bypass the
  basecall pipeline, so a lossy copy-through shows there), and every downstream
  BAM, `cl` call and summary table.

  It carries a **control** run whose only job is to measure what the pipeline
  does *not* reproduce on its own, because a first cut of this test blamed the
  resume path for six differing summary tables that turned out to be bwa picking
  arbitrarily among equal-scoring near-identical tRNA isodecoders — the same read
  landing on `tRNA-Gln-CTG-1-1` one run and `tRNA-Gln-TTG-1-1` the next, with an
  identical charging score. Artifacts that move in the control are reported;
  artifacts that hold still must be identical after a resume. The control also
  confirms the uBAM itself is bit-reproducible, so that wobble is downstream of
  dorado.

  Two traps are recorded in the script header because each produced a convincing
  false pass while it was written: SIGTERM to snakemake *alone* is not what Slurm
  does — snakemake treats it as a graceful stop and waits for the running job, so
  the basecall completes — and deleting a rule's output does not make it re-run
  when its consumers are up to date, so the DAG exits 0 having done nothing and
  every downstream comparison passes against artifacts nothing regenerated.

- **`rebasecall` and `rebasecall_ldx_run` resume** (#149). Both now run through
  `workflow/scripts/dorado_basecall_resume.sh`, which basecalls into a sibling
  dot-file and renames it onto `{output}` only on success, then hands that
  partial back to dorado as `--resume-from` on the next attempt. A job killed on
  wall clock costs the tail of a basecall and a requeue instead of the whole
  flowcell.

  Two things had to be worked around. The partial cannot be the rule's output:
  Snakemake deletes a failed job's outputs ("since they might be corrupted") and
  the basecall is additionally `maybe_temp(tier="basecall")`, so the failure path
  destroyed the exact artifact `--resume-from` exists to consume. And dorado will
  not read and write the same file, so the checkpoint is rotated to a second name
  before it is passed back. `clean` still reclaims both — it removes the whole
  `bam/rebasecall{,_run}` directory.

  What `--resume-from` does was read off dorado's sources rather than assumed,
  because guessing it wrong loses reads silently: `ResumeLoader::copy_completed_reads()`
  pushes every record of the resume file into the writer, so the **resumed BAM is
  self-contained** and the two files are never concatenated; the read ids it saw
  are then excluded from the new basecall, via the `pi` parent tag for split
  reads; and its record loop treats a read failure as end-of-records, so a BAM
  cut off mid-write is the designed input, not a corruption. Confirmed against
  dorado 2.1.1 on the LDX fixture.

  The wrapper refuses a checkpoint it should not continue. dorado re-parses the
  checkpoint's `@PG CL` and errors if the models differ — correct, but leaning on
  it would wedge the rule on every retry after a model change, so the argv is
  compared first and a superseded checkpoint is dropped. A checkpoint whose
  header never landed is dropped the same way, and of two survivors the larger
  wins, because an attempt killed inside dorado's copy-through leaves a partial
  shorter than the file it resumed from.

### Changed

- **`rebasecall_ldx_run` walltime 285m -> 600m**, and the sizing note inverted
  with it. 285m was derived on `rna004_130bps_sup@v5.3.0` at `gpu:4` (1,632
  reads/s wall); the rule now runs v0.5.0's `rna004_sup@v6.0.0` at `gpu:2`, where
  the measured rate is 931 reads/s (job 256074, 8,646,049 reads in 02:34:44),
  putting a 30.7M read flowcell at ~9.2 h. The old advice — size to finish,
  because an overrun destroys everything — is what forced walls to be padded
  rather than measured; with resume the number can be the measured need, and an
  overrun costs a requeue. `escapepod_demux` still has no resume (escpod has no
  `--resume-from`), and its note now says so on its own rather than by reference
  to this rule.

- **Right-sized the two GPU rules against cgroup measurements, and taught the
  Slurm profile to poll with `squeue`.** No behaviour changes, only what the
  rules ask the scheduler for.

  `escapepod_demux` drops from 32 to **16 cpus** and `rebasecall_ldx_run` from
  32 to **12**. Both are measured from the jobs' own `cpu.stat usage_usec`
  rather than from MaxRSS: demux draws a mean of 1.86 cores (peak 4.24) and is
  bound by filesystem read latency, while dorado draws a flat 4.67 (peak 4.84)
  and is GPU-bound. Neither was using what it held.

  `rebasecall_ldx_run` drops from `gpu:4` to **`gpu:2`**, which is a scheduling
  decision and not a throughput one. Asked of the scheduler with `--test-only`,
  every other field held constant, `gpu:1/2/3` all started immediately and
  `gpu:4` started **12.5 h later** — a whole-node request waits for a node with
  nothing else on it. On an idle partition 4 is strictly better; re-derive per
  site.

  `slurm-status-command: squeue` is new. Without it the executor's status check
  needs slurmdbd, and when that is unreachable the controller does not error —
  it silently stops making progress after its first submission batch, which is
  what a controller driven from inside an allocation hits.

- **Documented why `escapepod_demux` memory must NOT be cut**, after a cut to
  24000 cgroup-OOM-killed two flowcells. escpod's arena grows to fill the
  cgroup, so its 99.7%-of-request MaxRSS states no requirement — but every
  surviving job also sat at 24.44 GB, because with an arena the limit *becomes*
  the working set: lowering it moves the cliff rather than reclaiming slack. The
  check that separates a safe cut from a fatal one is `memory.stat` `anon` vs
  `file`, now recorded at both keys — escpod is ~24 GB of anon (a cliff), dorado
  is 5.3 GiB anon against 111.1 GiB of reclaimable page cache (not one), which
  is why its 128 G looks equally oversized and is left alone.


## [v0.5.0] - 2026-09-03

### Changed

- **BREAKING: the pipeline basecalls with `rna004_sup@v6.0.0` and scores charging
  with `charging_feature_nn_sup6_rna004@v0.1.0`.** Previously
  `rna004_130bps_sup@v5.3.0` + `charging_feature_nn_rna004@v0.1.1`. Both halves
  move together because they are one pairing, not two settings — see below.
  **Charging calls, modification calls and base-calling error rates all change**;
  results are not comparable across this switch without re-running.

  The v6 charging bundle was built with dorado `2.1.1+d66c17c` on
  `rna004_sup@v6.0.0`, which is exactly the binary `pixi run setup` installs and
  exactly the model `base_calling_model` names — verified byte-for-byte (190
  files, `9cab42f3…bbadb1`). So the pairing is now exact where before the
  vendored bundle declared dorado 1.4.0 against a pinned 2.1.1, and the
  major-version warning added in the previous release is silent.

  It is a **trade**, not a free upgrade. Against the v5.3.0 pairing the v6 one
  loses 0.0034 test AUROC (0.9906 -> 0.9872) and ~1 pp of balanced accuracy
  (0.9625 -> 0.9530), and cuts the fraction of CHARGED reads excluded as
  unreadable from 4.03% to **1.52%**. That last row is why it is the default:
  abstention is charging-correlated and both bundles' own `coverage_note` says a
  charging fraction over called reads alone is an UNDERESTIMATE, so a 2.6x cut
  in the charged-class exclusion rate reduces a systematic bias in the headline
  number — which matters more here than 0.0034 of AUROC. Report the no-call rate
  beside the fraction either way.

  The v6 bundle's operating point is **measured on the v6 corpus with its own
  checkpoint**, not carried over: `operating_point.cl` is 200, the same value
  `charging.ml_threshold` already used, but arrived at independently.

  The four modification models follow automatically —
  `download_mod_models` builds `"{dorado_model}_{mod}@v1"`, and
  `rna004_sup@v6.0.0_{m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG}@v1` all
  exist. Demux is untouched: `escpod demux` runs its own CTC-CRF over raw signal
  and never sees a dorado basecall, so no demux bundle declares a basecaller and
  none needed to move.

  `charging_feature_nn_rna004@v0.1.1` **stays vendored** for data already
  basecalled with v5.3.0 — re-basecalling an existing run to reach v6 is a GPU
  cost, not a correctness fix. Using it means setting `base_calling_model`,
  `dorado_model` and `charging.model` together; the basecaller check added last
  release is what stops you setting only some of them.

### Fixed

- The dorado major-version warning claimed "the basecalling model matches" even
  when it was firing alongside a model-identity error that said the opposite.

### Added

- **The basecaller a charging bundle was trained on is now checked, not just
  printed.** Since escapepod-models#106 a charging bundle declares
  `{model, model_sha256, dorado_version}`, and escpod states it at load without
  enforcing it. Two rules, deliberately of different strength:

  - **Model identity is an error.** The charging feature set is
    `mean + z-scored k-mer residual` with the expected level predicted from the
    read's own basecall, so the model substantially detects *how the basecaller
    fails* at the aminoacyl adduct — a different basecalling model is a domain
    shift, not a detail. Upstream measured ~0.0097 AUROC, 3.0–3.2 pp of TPR and
    **3.9% of per-read calls flipped**, while the aggregate charged fraction
    moved 0.04 pp: the one number anyone would check reads "no change".
    `charging.basecaller_check` (`error` by default, `warn`, `off`) governs it;
    `warn` is the setting for a deliberate cross-basecaller run, where
    arm-to-arm contrasts survive the shift and absolute charged fractions do not.
  - **Dorado version only ever warns**, at a major difference. The same weights
    run by a later dorado are still the same weights, so it cannot justify
    blocking a run — but a major bump is a different implementation and can move
    basecalls with the weights unchanged. **This fires today**: the vendored
    bundle was built with dorado 1.4.0 and `dorado_version` pins 2.1.1.

  Both run while the DAG is built, so a mismatch costs a dry-run rather than a
  basecall plus a classification pass. A bundle that declares nothing produces
  no findings — "cannot tell" is not "invalid".

- **`pixi run verify-basecaller`** proves byte identity against the declared
  `model_sha256`, reproducing upstream's hashing scheme exactly (sha256 over
  every file, name-sorted, each contributing its relative path then its bytes).
  That catches what a name comparison cannot: a model directory named right that
  is not the same bytes — a partial download, a re-fetch of a retagged upstream
  model, an edited `config.toml`. It reads ~300 MB, so it is a task rather than
  part of DAG construction. The vendored `rna004_130bps_sup@v5.3.0` verifies:
  159 files, `2f3e0926…3672f9`, byte-identical to what the bundle names.

## [v0.4.0] - 2026-09-03

### Changed

- **Charging model repinned `charging_feature_nn_rna004@v0.1.0` -> `@v0.1.1`,
  which needs escpod >= 0.19.0.** A sidecar-only reissue: the ONNX graph and the
  k-mer table are byte-identical (`cmp` clean; the graph's sha256 is unchanged),
  `metadata.json` differs only by the version string and a new top-level
  `basecaller` block, and the fixture scores **identically, read for read**
  under either bundle. No call moves.

  What it buys is that the block is readable. The charging feature set is
  `mean + z-scored k-mer residual` with the expected level predicted from the
  read's own basecall, so the model substantially detects *how the basecaller
  fails* at the aminoacyl adduct — swapping basecaller costs ~0.0097 AUROC and
  flips **3.9% of per-read calls** while the aggregate charged fraction moves
  0.04 pp, so the one number anyone would check reads "no change" while one read
  in 26 answers differently. The bundle now names what it was trained on
  (`rna004_130bps_sup@v5.3.0`, dorado 1.4.0) and escpod states it at load.
  `dorado_model` matches; escpod does **not** enforce it.

  This is a hard escpod floor, not a preference: a charging bundle's schema is
  `deny_unknown_fields`, so every escpod before 0.19.0 refuses the file outright
  with ``unknown field `basecaller` ``. Verified both directions — 0.19.0 scores
  the fixture, 0.18.1 will not open it. `pixi run check-models` now reports every
  vendored bundle as current or deliberately pinned.

### Fixed

- **LDX demux was broken for every sample on the shipped default bundle.** A
  sample configured `ldx01` was rewritten to `nbc01` before being matched
  against `classifications.csv` — but the default bundle since v0.3.0 is
  `barcode_crf_ldx16_rna004@v0.1.0`, whose references *are* named
  `ldx01`..`ldx16`, so the rewritten name appeared nowhere in the file. Every
  barcoded sample on the run then failed, and it failed **late**: the demux pass
  and the whole-run dorado basecall had already been paid for by the time
  anything looked at a barcode name. The rename was correct for the nbc16 bundle
  it was written against and became wrong when the default moved; a dry-run
  could not catch it, because both call sites merely produced a mapping that
  matched nothing rather than raising.

### Changed

- **`escpod` bumped 0.18.1 -> 0.19.0, and the LDX path now writes ONE `.p5s`
  per POD5 directory instead of one per file.** `escpod demux --annotate`
  pointed at a *directory* writes a single **collection** sidecar beside it
  (`<run>/pod5/` gets `<run>/pod5.p5s`); pointed at files it still writes one
  each. So the release alone does not get you the new shape — the argument
  does, and `escapepod_demux` now passes the POD5 directories rather than a
  `*.pod5` glob. A run that produced fifty POD5s produced one set of barcode
  calls, not fifty, and this is the layout that says so. Declared as
  `escapepod:p5s_version` 3; an older escpod refuses a collection *by name*
  rather than misreading it, so a mixed-version tree fails loudly. Consumers are
  unchanged — a POD5 is matched against its own `.p5s` first and then the
  collection beside its directory, and only when footer UUID and byte size match
  its member entry.

  One migration note: a per-file `*.pod5.p5s` left by a pre-0.19.0 demux
  **shadows** the collection, because the file's own sidecar wins per column.
  Nothing here reads a sidecar (the split is driven by `classifications.csv`),
  so pipeline output is unaffected — but delete the old per-file sidecars before
  pointing `escpod demux split --sidecar` or the Python `Reader` at a run
  demuxed under both versions.

- **`merge_run_pods` is gone: an LDX run no longer merges a second full copy of
  its signal.** The rule existed because `classify` takes one path and a run
  split over `pod5_pass`/`pod5_fail` is two — but the directory argument is
  walked **recursively**, so naming the run directory covers both. Verified
  byte-identical: the same 89 reads, same calls, whether `classify` is handed
  the flat POD5 directory or a run root containing it. Naming the run is a
  superset of its POD5 directories and that is harmless here, because `classify`
  is driven by the BAM and looks each aligned read up by id, so signal it is
  never asked for is never touched. This removes a full duplicate of the run
  from `pod5/runs/`, which is the exact cost the LDX path exists to avoid.
  (Recursion is not new in 0.19.0; the pipeline simply never used it.)

- **`escpod signal classify` is spelled `escpod classify` again**, upstream's
  current name. Not a floor: 0.18.1 already carried `escpod classify` as an
  alias in the other direction, so both spellings work on both releases. The
  only user-visible difference is the `@PG` line on the output BAM.

  0.19.0 also accepts a top-level `basecaller` block in a charging bundle's
  metadata, which every escpod before it refused outright. That matters for the
  *next* charging bundle rather than the vendored one — `charging_feature_nn_rna004@v0.1.1`
  is a sidecar-only reissue with byte-identical weights that declares it — and
  is why taking that bundle needs this bump first.

- **The `nbc` barcode vocabulary is retired, and nbc bundles are refused.**
  `barcode_crf_nbc16_rna004@v0.2.0` is no longer vendored, `ldx.model` rejects
  any `barcode_crf_nbc*` bundle by name, and `nbc` is gone from the prefix table
  in `workflow/scripts/barcode_names.py` — the LDX panel is `ldx01`..`ldx16` end
  to end, in the samples file, in the bundle and in the BAM, with no translation
  step to get wrong. This is not a naming preference: nbc16 and ldx16 are
  separate retrains whose calls differ on ~8% of reads, so keeping the retired
  bundle selectable meant keeping alive the crosswalk that caused the failure
  above, in exchange for an option nobody should take.

- **Barcode names are resolved against the bundle, not guessed from a prefix.**
  `resolve_to_bundle()` asks the configured bundle's `metadata.json` what it
  calls its references and raises `BarcodeVocabularyError` — naming the bundle's
  own vocabulary — when a configured name matches none of them. Both sites that
  compare against `classifications.csv` now use it, so a samples file and a model
  that disagree fail during DAG construction, at `pixi run dry-run-ldx`, instead
  of after the GPU hours. `bundle_barcode_names()` already existed for this and
  was called only from tests.

### Added

- **`pixi run check-models`** reports vendored model bundles that upstream has
  since superseded. `verify-demux-model` answers "is this copy what upstream
  released?"; this answers "is that still the release we want?", which a bundle
  can fail while verifying perfectly. It is a report and never an upgrade —
  newest is not best here, and `resources/models/pins.yml` records deliberate
  holds with their reasoning, so an undecided drift is distinguishable from a
  decided one. Needs network and a GitHub login, so it is a task rather than a
  workflow step: upstream is private and compute nodes have no route to it.


## [v0.3.0] - 2026-08-31

### Added

- **The WarpDemuX panel can be demultiplexed by `escpod` on the sidecar path.** A WDX-barcoded run no longer has to route every read into a second full copy of the POD5: point `ldx.model` at the newly vendored `barcode_crf_wdx4_rna004@v0.2.0` and the run is demuxed by the same `escapepod_demux` rule the LDX panel uses, writing a few-MB `.p5s` sidecar beside the raw signal and basecalling the run once. **No new rules were needed** — `get_ldx_model` asks only for a CRF bundle directory and the rule passes neither `--barcodes` nor `--method`, because the bundle carries its own references, so selecting the panel is a path in config. The bundle is byte-for-byte upstream, verified against the hash `escapepod-models`' `MANIFEST.json` declares and paired with `adapter_rna004@v1.1.0`, which is the boundary model its provenance requires: serving a CRF one it was not trained against silently changes the reading frame.

  **Provisional on two counts, both documented beside the bundle.** The suggested `min_crf_margin: 2.0` is measured on a run that is *in* this model's training corpus — the v0.2.0 split is a random read-level 10%, not leave-one-run-out, and every WDX4 run available to this project is in it — so it is an upper bound to be re-derived against a genuinely held-out run. And the labels it learned are ungated: upstream sets `no_gate: true` and records the cost, WarpDemuX declining ~15% of reads at ~0.68 balanced accuracy on those. Expect a gated retrain to supersede it. **`barcode11` is not in the panel** and cannot be added by configuration; a WDX4b run using it has no escpod equivalent today.

- **Tooling for the measurement that gates that switch.** `workflow/scripts/demux_concordance.py` compares two backends' routing read by read, and the decomposition is the point: a read the CRF never decoded (`refused` — its adapter ended before the model's window, a `boundary_margin` question) is not one it got wrong (`differed`), and a read WarpDemuX rejected that the CRF called (`claimed`) is yield rather than disagreement. On a 1.7M-read run the same data reads 87.63% pooled and **99.39%** decomposed, so this is not a presentational nicety. Agreement is also reported stratified by WarpDemuX's own confidence, because the CRF's labels are WarpDemuX calls at `conf>=0.9` and scoring it over every call that merely clears WarpDemuX's per-class reject threshold — as low as 0.173 for bc07 — measures the student on reads the teacher never taught it. `scripts/crf_margin_sweep.py` derives a `min_crf_margin` operating point from a run rather than inheriting the LDX panel's, and `workflow/scripts/barcode_names.py` holds the name crosswalk both need.

- **`reference.max_mismatch`** merges near-identical tRNAs when building the reference, not only exact duplicates. Genome-wide sets are heavily redundant — danRer11's GtRNAdb mature-tRNA FASTA is 8,879 records but 3,315 distinct sequences, and many of the rest differ by one or two bases, close enough that reads cannot be confidently assigned between them. At 2 mismatches the reference goes 3,315 → 1,246 entries. Grouping is greedy leader clustering, so every member is within the threshold of its representative rather than chained to it through intermediates; distance is measured before CCA is appended (CCA is added to 8,008 of 8,879 danRer11 records, so clustering after would silently regroup the equal-length buckets); and only equal-length sequences are compared, so an indel never merges. **Lossy and off by default** (`0`): merged tRNAs share one reference name, so charging and CPM cannot be resolved between members of a group. Every merge is recorded in `reference/build_report.txt`.

- **Per-read base-calling error calls and per-site charging-error counts**, behind `mismatch_calls.enabled` (off by default; the BAM walk is expensive). A charging odds ratio needs to know, per read, both whether it carries an error at a site and whether it is charged — `bcerror` gives the first only in aggregate and carries no read identity, so it can *select* sites but never supply the evidence to test them. `get_mismatch_calls.py` reports both, its primary output being a per-site error × charging 2×2 (488 KB against 52 MB for the per-read form on a zebrafish sample); `select_bcerror_sites.py` picks the sites once across all samples so every sample is tested at the same positions.

### Changed

- **`escpod` bumped 0.17.1 -> 0.18.1, a correctness floor for `signal classify` on any reference whose tRNA bodies are not sacCer3's.** The CCA|adapter junction was located by requiring `CCAGGC` to occur exactly once per reference record, and that check ran *before* the common-arm check that actually identifies the junction. The motif is only CCA plus the adapter's opening bases, so it collides over ~75 nt of tRNA body in ~1.5% of records — 4 of 282 in hg38, 51 of 3315 in danRer11 — and one such record aborted the entire run for that sample. Dropping the offenders would not have been neutral: **48 of danRer11's 51 are Glu**, so it would have removed most of one amino-acid family from a charging analysis. Candidates are now filtered by the arm first and uniqueness required of what survives. Two more fixes come with it: 0.18.0 stopped the reads table being written as a single Arrow record batch however large it was (`demux` shards its reader threads by batch index, so a single-batch file is decoded by one thread), and 0.17.2 made every POD5-taking command accept a directory and refuse a path that does not exist — before it, a directory input logged one WARN, wrote a header-only CSV and **exited 0**, which nothing downstream can tell from a run where no read passed.

- **The barcode vocabulary is panel-aware.** Every escpod CRF bundle emits whatever its `metadata.json` calls its references, and the two panels differ in *which* of a sample's two names is configured: LDX samples are written `nbc01`, upstream's name, and renamed to `ldx01` for the BAM, while WDX4 samples are written `barcode03`, ours, and the bundle emits `bc03`. `get_sample_barcode_label` now renames on the strength of the name rather than of `is_ldx_enabled()` — a guard that was safe only while `nbc` could mean one thing — and the two sites matched against `classifications.csv` key on the emitted form. Without that, a WDX4 sample configured `barcode03` matches nothing in a file full of `bc03` and *every* sample on the run fails with "No reads were assigned", after the demux pass and the whole-run basecall have already been paid for.

- **`rebasecall_ldx_run` is sized to finish, and gets every GPU** (`gpu:4`, `285m`, 128 GB, 32 cpu, from `gpu:1`/`24h`/48 GB/8 cpu); `escapepod_demux` goes 16 -> 32 cpu. The rule has **no resume** — it is a plain `dorado basecaller ... >{output}`, so a job killed on wall clock loses the whole basecall rather than a tail of it, and the profile comment claiming dorado would resume was false comfort. That matters because Slurm will not *start* a job whose walltime runs into a maintenance reservation, so an over-large request does not buy margin, it means the job silently never runs — which is what happened to the 2026-08-21 flowcell, after `escapepod_demux` had already spent 4h43 on the GPU. The walltime is now a worked number: 1,632 reads/s wall measured on 4× A30, ~25.5M assigned reads ≈ 4.3 h, plus 25 min. It encodes a date rather than a property of the data, and should be re-derived against the next reservation.

- CI: `prefix-dev/setup-pixi` 0.10.1 -> 0.10.2.

### Fixed

- **The v6 dorado basecalling model was not ignored.** `resources/models/.gitignore` matched `rna004_130bps_sup*`, but dorado v6 dropped the `130bps` segment, so `rna004_sup@v6.0.0` — 365 MB of tensors — showed up as untracked and was one `git add -A` from being committed.

## [v0.2.1] - 2026-08-26

### Fixed
- **`get_align_stats` and `get_charging_summary` closed the interpreter's stdout.** Both called `fout.close()` unconditionally, but `fout` is `sys.stdout` whenever `--out` is not given. Harmless in the pipeline, which always passes `--out` and where the process was exiting anyway, but wrong for any interactive or piped use. Both now borrow stdout through a context manager and only close a handle they opened.
- **`get_trna_charging_cpm` never closed its output file.** It relied on CPython refcounting to close the handle at function exit, which happens to work and would not under a different interpreter, or if an exception kept the frame alive. On the gzip branch that is a truncated file, since the trailer is only written on close. Now a `with` block.
- **`scripts/install-escpod-gpu.sh` verified its own build with a check that 0.17.1 made impossible.** It grepped `escpod demux --help` for `--gpu` and failed the build if absent — but `--help` is now byte-for-byte identical between the musl release and the GPU build, since clap renders `--device <auto|cpu|gpu>` whether or not the gpu features are compiled in. No flag grep can distinguish them, so a correct GPU build would have been rejected as "the features did not take". It now greps the binary for the `CUDAExecutionProvider` string, which is linked in only when the features actually took (0 occurrences in the musl build, present in the GPU one).

### Changed
- **`escpod` bumped 0.12.0 -> 0.17.1, and `ldx.gpu` no longer needs a source build.** 0.17.1 is the first release to publish a GPU artifact (`escpod-v<ver>-x86_64-unknown-linux-gnu-gpu.tar.gz`), so `pixi run setup` now downloads both it and the portable musl binary, each against a pinned SHA256. `scripts/install-escpod-gpu.sh` survives only for building a ref that has no release. The GPU artifact is the one dynamically linked build (glibc >= 2.28, because the CUDA runtimes are dlopened and a static musl binary cannot) and is x86_64 Linux only, so setup skips it elsewhere rather than failing; `escpod_version` stays shared with the musl default, since pointing the global pin at `-gpu` would hand every other rule a single-platform binary just to run `escpod merge`. **Do not pin 0.17.0**: its GPU artifact failed to build, which skipped the release job, so it exists on PyPI with no GitHub release and no binaries — upstream superseded rather than corrected it, and `pixi run setup` would 404.
- **`escpod demux --gpu` migrated to `--device <auto|cpu|gpu>`.** The old spelling is deprecated but still works (it warns and continues as `--device gpu`), so this is a migration rather than a break — but the meaning changed underneath it: `--gpu` fell back to the CPU when the device was unusable, where `--device gpu` fails. `escapepod_demux` now passes `--device` explicitly in **both** directions rather than relying on the `auto` default, because the failure modes are not symmetric — on this path a silent CPU fallback is a 20x slowdown (2.6 h against 458.9 s on our own flowcell) that still produces correct output, so it looks like success and is caught only by noticing the wall clock. Under `ldx.gpu: false` the explicit `--device cpu` likewise stops a GPU binary from opportunistically using a device the config said not to.
- **Verified on an A30** (driver 580.126.09, our `gpu_cuda13` onnxruntime 1.27.1 against the CUDA 12 runtime upstream documents): the encoder reports `CRF encoder: 2 worker(s) on GPU [0]`, and GPU/CPU agree on **414 of 415** reads of the LDX fixture. The single disagreement is the lowest-margin read in the fixture (`crf_margin` 0.2239 on GPU against -1.2088 on CPU, where the median is 9.55) and both paths pick the same `crf_best`, i.e. float nondeterminism between encoder paths rather than a behaviour change. It is exactly the read `ldx.min_crf_margin` exists to gate.
- **`escpod` 0.13.0 is a new correctness floor for the LDX path**, and the reason this bump is worth taking on its own. A `.p5s` sidecar's `(batch_idx, row_idx)` locators were dereferenced unchecked, so an index that passed the file-level identity guard but held wrong offsets returned *a different real read, correctly self-labelled* — undetectable downstream. The signal paths were worse: they projected `read_id` and never read it, stamping the queried UUID onto whatever row they landed on. Every indexed lookup now confirms the row's `read_id` first. Our `demux --annotate` -> sidecar -> `signal classify` flow is precisely that path.
- Lint gates for real. All three jobs in the `Lint` workflow (`python-lint`, `yaml-lint`, `snakemake-lint`) now fail the build; the first two were `continue-on-error: true` and had been reporting success over 7 unformatted files and a genuine `empty-lines` error in `ci.yml`. Ruff's rule set is declared explicitly in `ruff.toml` (`target-version = "py310"`, plus `I`/`B`/`SIM`/`FURB`/`C4`/`BLE`/`EXE`/`RUF`) and yamllint's in `.yamllint.yml`, because tool defaults are not a stable contract — ruff 0.16 widened its own and took `ruff check` from clean to 49 findings on unchanged code. `ruff` is pinned `0.16.*` (it was `"*"`, resolving to three different versions across environments in `pixi.lock`), and `.pre-commit-config.yaml` now runs the same tools at the same versions CI does: black and flake8 are gone, having been a third and fourth opinion about `workflow/scripts/` that no installed hook had ever actually run.

## [v0.2.0] - 2026-08-23

### Changed
- **BREAKING: `ldx.gpu` now defaults to `true`.** GPU demux is 20.3x on our own flowcell (458.9 s against ~2.6 h for 1,001,307 reads, identical 92.22% yield, on a quarter of the cores), so the CPU path is no longer a sensible default. It requires a source build — the released escpod carries no GPU code and rejects `--gpu` as an unknown argument — plus a CUDA onnxruntime and cuDNN. `escapepod_demux` resolves the locally built `<escpod_version>-gpu` binary itself and exits naming the build command if it is absent; `escpod_version` deliberately stays on the release, since demux is the only rule with a GPU path and moving the global pin would force a source build on every rule and break `pixi run setup`'s download URL. Set `ldx.gpu: false` for a CPU-only host — `config-ldx-test.yml` does, because the value is resolved while the DAG is built and CI has no GPU. `cluster/lsf/config.yaml` gains the `escapepod_demux` GPU entries it never had; the SLURM profile already had them.
- **BREAKING: amino-acid identity classification and the `leech` dependency are removed.** The `classify_aa` and `classify_aa_identity` rules, `workflow/rules/aatrnaseq-classify-aa.smk`, both config blocks, their `pipeline_outputs()` branches, the `bam/aa_classified/` output, and the GPU entries for them in the LSF and SLURM profiles are all gone, along with `scripts/install-leech.sh`, the leech section of `setup-tools.sh`, the `install-leech` pixi task and the `resources/leech` submodule (the repo now has no submodules). Both features were `enabled: false` by default, so the default pipeline is unaffected; charging classification never used leech — it runs `escpod signal classify` against the vendored ONNX bundle. Recovering the feature means restoring it from git history and installing leech and torch independently.
- `scripts/install-escpod-gpu.sh` no longer builds from `resources/leech/escapepod-rs`. escapepod-rs was a submodule *nested inside* leech, which coupled the GPU build to an unrelated dependency; it now clones standalone into `resources/escapepod-rs`, or builds from whatever `ESCPOD_SRC` points at.
- **BREAKING: LDX demux adopts `barcode_crf_ldx16_rna004@v0.1.0`**, the first bundle whose class names are `ldx01`..`ldx16` rather than upstream's older `nbc` vocabulary; the `barcode_crf_nbc16_rna004` family is closed at v0.3.1. It is a **retrain**, not a rename — corrected geometry (`state_len=4`, so the full 27-nt code is emitted rather than discriminating on 23 nt) and the first trained without bonito or ont-koi — so **demux results are not comparable across this switch**. Measured on the full 1,001,307-read run at identical settings, ~8% of calls change; yield is identical (92.22%) in all comparisons, so this is purely *which* barcode and never *whether*. That churn is not a property of ldx16: an nbc-only upgrade (v0.2.0 -> v0.3.1) moves 8.24%, and the disagreements are not cumulative across releases, which is the signature of retrain variance.
- **Vendored demux bundles are byte-for-byte upstream again.** `boundary.margin: 0` and `boundary.clamp_max_shift: 300` had been patched into nbc16's `metadata.json` without regenerating `SHA256SUMS.txt`, so that bundle failed its own integrity check for eleven days — and the patched copy then read as evidence that a later release had *removed* those keys. It had not: no upstream CRF bundle has ever declared either, at any version, because `build_crf_bundle.py` cannot write them. Both values now live in `config-base.yml` where they are visible and diffable, with `--clamp-max-shift` wired in `demux.smk` to make the move behaviour-neutral. Unset does not mean "the bundle decides" — it means escpod's fallback of margin 200 and no clamp, which costs 24 of 415 fixture reads.
- **`ldx.min_margin` is retired to 0, superseded by the lattice gate.** It briefly defaulted to 12 — the sub-plateau tail is almost entirely wrong calls, so 12 dropped 8,817 reads (0.95%) of which 85.2% were cross-bundle disagreements, while 13 was a trap that ate 46.76% of the run at 8.5% precision (99.0% of reads sit on a 12/13/14 margin plateau, the 16 references being >=12 apart by design). But its ceiling was always ~1% of reads. Re-measured against `min_crf_margin: 1.0` on the same flowcell, `min_margin 12` removes only 728 further reads, **90.5% of which are correct**, for +0.09 points of total error removed — so it is not complementary to the lattice gate, it just lowers the average quality of what gets discarded. Raise it again only for a run with no lattice gate, and leave it 0 when the measurement *is* the crosstalk. `config-ldx-test.yml` no longer pins `min_margin` or `ldx.model`, so the tests exercise the shipped defaults rather than bypassing them.
- **The charging classifier is now `escpod signal classify` against a vendored ONNX model bundle, and Remora is retired.** The pipeline no longer installs, imports or invokes Remora anywhere.
  - The model is `charging_feature_nn_rna004@v0.1.0`, vendored at `resources/models/charging/` (upstream `rnabioco/escapepod-models` is private and compute nodes have no route to GitHub). It is a bundle **directory** and is self-describing: the anchor, feature recipe, k-mer table (pinned by sha256), abstain rule and recommended operating point all come from its `metadata.json`, not from flags. Published test AUROC 0.9906, balanced accuracy 0.9625 on a held-out flowcell.
  - The bundle's `9mer_levels_v1.txt` is a symlink to `resources/kmers/9mer_levels_v1.txt`, which is byte-identical to the released copy — the file the model was trained against. Replacing that file invalidates the model; the sha256 in `metadata.json` gates it.
  - New config block `charging` with `model`, `min_mapq` and `ml_threshold`. `remora_cca_classifier`, `remora_kmer_table`, `classifier` and `opts.remora` are gone.
  - `charging.min_mapq` defaults to **0**, not escpod's own default of 1. tRNA references are highly redundant, so a read mapping equally well to two isodecoders gets MAPQ 0 from bwa and is still a good read; `--min-mapq 1` drops 118 of 209 records (56%) on the test data, while 0 classifies exactly the 187 reads Remora did.
  - `charging.ml_threshold` replaces the value hardcoded in `get_cca_trna_cpm`. It is 200, matching the bundle's declared operating point.
  - `classify_charging` is CPU-only; `escpod signal classify` has no GPU path. GPU entries removed from the cluster profiles.
- **`transfer_bam_tags` is removed, and there is no `cm` tag.** Remora emitted its score into `MM`/`ML` — the standard modbase tags — clobbering the calls modkit needs, so a rename step existed purely to move them to `cl`/`cm`. `escpod signal classify` writes `cl` (uint8, `round(P(charged) * 255)`) directly onto the input records, in the same order, leaving dorado's `MM`/`ML` intact. Nothing downstream consumed `cm`.
- **The `classifier: remora|leech` switch and the `classify_charging_leech` rule are removed.** They selected between two implementations of the same retired `cca_classifier.pt`, which is also deleted. leech's amino-acid identity rules (`classify_aa`, `aa_identity`) are unaffected.
- **`escpod` bumped 0.8.1 -> 0.12.0**, pinned alongside the model: the bundle is the per-base-feature ONNX variant, which escpod reads only from 0.10.0 onward (older binaries refuse it with ``missing field `gbm` ``). 0.11.0 additionally applies the bundle's abstain rule (previously parsed and dropped, so ~0.85% of scoreable reads got a confident call they should not have) and adds the `--tsv` `reason` column; 0.12.0 is required for the demux CRF scores below.
- **`escpod` is now checked unconditionally at startup**, not only for LDX runs — `merge_pods` and `classify_charging` use it on every run. The `pod5` CLI check moved under the warpdemux conditional, which is the only remaining consumer of that package.
- `pixi run setup` no longer installs PyTorch. Nothing in the default pipeline needs it; leech's `.pt` bundles still do, so enabling `classify_aa`/`aa_identity` now means installing torch yourself.
- Updated dorado from 1.4.0 to 2.1.1. The basecalling model is unchanged: `rna004_130bps_sup@v5.3.0` remains the newest RNA004 sup model in dorado 2.1.1, and all four modified-base models the pipeline uses (`m5C_2OmeC`, `inosine_m6A_2OmeA`, `pseU_2OmeU`, `2OmeG`) are still available at v5.3.0.
- Updated dorado from 1.3.1 to 1.4.0, basecalling model from v5.1.0 to v5.3.0, and modkit from 0.6.0 to >=0.6.1.
- Shell scripts (`setup-env.sh`, `setup-tools.sh`) now read `dorado_version` and `dorado_model` from `config/config-base.yml` instead of hardcoding defaults.
- Updated stale documentation references for dorado and modkit paths/versions.
- Pipeline now warns on startup if installed dorado version does not match `dorado_version` in config.

### Added
- **LDX demultiplexing: a second signal-level barcode backend, `escpod demux`** (`ldx.enabled`). Where WarpDemuX classifies boundary-gated DTW fingerprints, this reads the barcode out of the raw adapter signal with a CTC-CRF and matches the decode to references by edit distance, in one fused pass that detects, basecalls, matches and routes each read straight into its barcode's POD5 — so the per-sample split already exists when the command returns, with no separate read-ID extraction or `pod5 filter`. Barcodes are assigned with an `ldx:` key in the samples YAML. Exactly one of `warpdemux.enabled` / `ldx.enabled` may be set; both populate the same per-sample barcode field, so enabling both is rejected at parse time rather than producing an ambiguous DAG. No barcode kit is configured — the model is a self-describing bundle *directory* carrying its own references and pinning the boundary detector it was calibrated against.
- **A barcoded test fixture, and an LDX path that is actually tested end to end** (closes #120). `.tests/fixtures/ldx-demux` is 415 reads (4.6 MB) cut from the 2026-08-06 LDX flowcell, selected from that run's own completed outputs so each read's barcode, 3' adapter and alignment status were known before it was picked — which is what lets it be this small without relying on chance survival. Five populations, each covering a branch: two EDX-filtered barcodes, one unfiltered pooled barcode, one barcode no sample claims, and reads the CRF could not call; the filtered barcodes deliberately carry wrong-adapter reads, or `filter_{fastq,pod5}_by_edx` would be a no-op that still passes. Committed rather than added to the S3 tarball, so it cannot drift from the configs referencing it. `tests/integration/test_ldx_demux.py` is two-tiered — 11 fixture/config checks run anywhere including CI, 22 output checks skip without a GPU run — and CI gains `pixi run dry-run-ldx`.
- **`BC:Z:` on every read of a demultiplexed BAM.** The barcode previously existed only in the output path: given `bam/final/ldx04/ldx04.bam`, nothing in the file said which barcode produced it, and the per-read record that could recover it is deleted by `clean` (and is `temp()` under `demux_scratch` on the WarpDemuX path, so it can vanish mid-run). Stamped at `inject_ubam_tags`, the first point both backends have converged. Unbarcoded samples get no tag at all — absence means "not demultiplexed", not "unknown". LDX values are normalised to this project's naming (`nbc04` -> `ldx04`), with upstream's name preserved in an `@CO` line so no reader has to guess the vocabulary.
- New `transfer_tags.py` flags: `--set-tag TAG:TYPE:VALUE` (a constant tag, which the script previously could not write at all — every tag had to exist on a source read), `--rg-sample` / `--rg-library` / `--rg-barcode`, and `--comment`.
- `pixi run verify-demux-model` checks demux bundles against the hashes **upstream** declares (`metadata.boundary.sha256`, `provenance.sha256`) plus any `SHA256SUMS` the release itself shipped — never a checksum file written here — and names `boundary.margin` / `clamp_max_shift` explicitly if a bundle declares them, since that can only mean the copy was edited.
- **LDX demux records the CRF's own confidence** (`ldx.ref_scores`, on by default; needs escpod >= 0.12.0). `classifications.csv` gains `crf_logp`, `crf_margin`, `crf_best` and `mean_logpost` — `log P(barcode | signal)` read out of the CRF lattice by constraining its forward recursion to the paths emitting each reference. This is the per-read score `confidence` never was: on a designed panel the edit-distance margin measures how far apart the references are, not how sure the model is, which is why `min_margin` plateaus (99% of reads on three values) and why 90% of cross-bundle disagreements are exact reference matches no edit-distance gate can see. Costs +3.6% on the fused command, more under `ldx.gpu`. Columns append after `confidence`, so `summarize_demux.awk` and the QC report are unaffected.
- **`ldx.min_crf_margin: 1.0`** — the pipeline's false-positive control, replacing `min_margin`. Below it a read becomes `unclassified` and keeps its scores, so the row still says which gate dropped it and by how much. Swept on the 2026-08-06 donor flowcell (1,001,307 reads, 923,431 calls), scoring precision against a second independently trained bundle — the same floor-on-error methodology as the old `min_margin` table:

  ```
  gate                 recall%   of all errors removed   discard precision%
  min_margin 12         91.34            9.52                  83.78
  min_crf_margin 0.5    89.53           29.87                  86.12
  min_crf_margin 1.0    87.85           47.69                  84.56   <- set
  min_crf_margin 2.3    85.59           64.48                  75.38
  min_crf_margin 4.6    80.18           77.71                  50.02
  ```

  1.0 discards reads that are 84.56% real errors — the same quality as `min_margin 12`'s 83.78% — while removing 47.69% of all error against that gate's 9.52%. Five times the error caught at equal discard precision, for 3.5 points of recall. Past ~1.3 discard precision decays; by 4.6 half of what is dropped is correct. Corroboration that the score means what it claims: 9,617 of 9,618 reads where the lattice prefers a different reference carry a negative margin, and 87.9% of those are genuine cross-bundle disagreements. Measured on the GPU encoder, unlike the CPU-measured `min_margin` table; calls agree with the CPU encoder on 99.76% of fixture reads at a median `crf_margin` difference of 0.0016 nats. `min_crf_prob` is left unset — a different cut, and not the one swept.
- **`summary/tables/{sample}/{sample}.charging_calls.tsv.gz`** — one row per read the charging model saw, with `p_charged`, `cl`, and a `reason` naming why every unscored read was unscored (`no_aligned_arm`, `no_signal`, `ns_mismatch`). Not optional: the model **abstains** on reads whose common arm did not align rather than guessing, abstention is charging-correlated, and so a charging fraction over called reads alone is an **underestimate**. Report the no-call rate beside it.
- `read_attrition` folds those reasons in, so the `aligned -> charge-called` gate names its loss instead of inferring it. `anchor_coverage` still covers the complementary half — reads that never reached the model at all.
- `pixi run verify-charging-model` verifies the vendored bundle against its pinned checksums.
- `trim_reference` rule produces a tRNA-only FASTA (`trna_only.fa`) by stripping 5'/3' adapter sequences from the adapted reference. This FASTA is used by clover for MODOMICS annotation and structure visualization.
- `build_trna_reference.py --mode trim` for generating adapter-stripped tRNA-only FASTA files.
- `get_bcerror_freqs.py` and `compute_odds_ratios.py` accept `--offset-5p` and `--offset-3p` to filter adapter positions and convert to tRNA-only coordinates.

### Fixed
- **The QC report ignored `charging.ml_threshold` entirely.** `aatrnaseq-report.smk` read `config.get("ml-threshold", 200)` — a top-level key that has never existed — so the lookup always fell through to the literal 200. Nothing broke while the threshold *was* 200, which is why it survived: the two agree by coincidence at the shipped value. Move the threshold and the CPM tables shift while the report keeps reporting the old cut, with no error and no warning. `get_cca_trna_cpm` takes it from `config["charging"]["ml_threshold"]`; now so does the report. `CLAUDE.md` and `docs/troubleshooting/faq.md`, which both still described the threshold as hardcoded in the rule, now also say what 200 *is* — the bundle's declared operating point (FPR 1.74% / TPR 0.939 against ligation chemistry) — and that its precision depends on the sample's own charged fraction, so a low-charging sample needs a higher cut.
- **Every final BAM was invalid SAM.** dorado emits a proper `@RG` carrying runid, basecall model, flowcell and device; `bwa_align` builds the aligned header fresh from the reference and drops it; `inject_ubam_tags` then copied dorado's per-read `RG:Z:` straight back on while writing with `template=target` — the bwa header, which declares no read groups. So every read in every final BAM referenced an `@RG` that did not exist. `samtools quickcheck` passes on this, which is why it went unnoticed, but Picard `ValidateSamFile` and GATK reject it, and the basecall-model provenance was destroyed at the header while a dangling pointer stayed on each read. `transfer_tags.py` now builds the output header explicitly, splicing the source's `@RG` in; dorado's `ID` is deliberately preserved (it is what `RG:Z:` points at) along with `PU`/`PM`/`DT`/`PL`/`DS`, and only `SM`/`LB`/`BC` are overwritten with the pipeline's sample, run id and barcode.
- **An LDX run's QC report claimed demultiplexing was not enabled**, and showed an empty Barcode column, while `escapepod_demux` was writing a perfectly good `demux_summary.tsv.gz` that nothing read. Three separate bugs, any one of which alone leaves the section empty: the gate was `isTRUE(config$warpdemux$enabled)`, false on an LDX run; the samples parser read `val$wdx` and `val$edx` but never `val$ldx`, so `barcode` was all NA; and `run_ids` was derived by filtering `!is.na(barcode_kit)`, which finds nothing on an LDX run because LDX deliberately configures no kit. Both backends already write the same `predicted_barcode / n_reads / pct` schema, so no plotting or table code changed — `wdx_summary_all` / `has_wdx_summary` are renamed to `demux_summary_all` / `has_demux_summary`, and a `demux_backend` variable names which one is in play.
- **`align_stats`'s `classified` row counted every read in the final BAM, not the reads the model actually called.** Under Remora the two were the same, because Remora only emitted reads it classified; `escpod signal classify` passes unscored records through untouched, so the row silently became equal to `aligned` and the charge-calling gate in `read_attrition` reported zero loss. `get_align_stats.py` gains `--require-tag`, and the rule passes `cl` for that BAM. On the test data the gate goes from a reported 0 lost to the true 68 of 624.
- **Bash strict mode was off for every rule.** The Snakefile's `onstart` calls `shell.prefix()` to put the pinned tool bins on PATH, and that REPLACES snakemake's default prefix — which is `set -euo pipefail`. Every `a | b | c` in the workflow was therefore reporting only `c`'s exit status, so a failing `bwa mem | samtools view | samtools sort` would look like a success with an empty BAM. The prefix now restores strict mode explicitly.
- Removed an unused `OrderedDict` import in `get_align_stats.py` that was failing `pixi run lint-python`.
- `get_charging_table.py` no longer drops reads whose charging ML tag is exactly 0. The write gate used `if tag_value` (falsy at 0), silently discarding maximally-confident *uncharged* reads (ML score range is 0-255, >=200 = charged). This biased charging fraction upward and shrank the CPM denominator in `get_trna_charging_cpm.py`. Gate now checks tag presence (`is not None`). Added a regression test covering ML==0.
- `get_charging_table.py` now warns on stderr (instead of dropping silently) when a tag is a multi-element array that cannot be reduced to a single charging score — the same denominator-biasing failure mode as the ML==0 bug. This is reachable via the default `--tag ML` on a dorado mod-base BAM, where every read would otherwise be discarded with no signal. Also warns when the output table ends up empty.
- Pipeline summary outputs (bcerror, odds_ratios) now report positions in tRNA-only coordinates (1-indexed) instead of full-reference coordinates that included adapter sequences. This fixes incorrect nucleotide positions in downstream tools like clover's `plot_tRNA_structure()`.

## [v0.1.1] - 2026-02-11

### Fixed
- `bwa_align` OOM at 48 GB: decoupled dorado tags from alignment by stripping tags from FASTQ, dropping `bwa mem -C`, and injecting tags from the unaligned BAM afterward via new `inject_ubam_tags` rule

### Added
- Quarto QC report with per-sample tabs (#81)
- Per-tRNA pairwise modification odds ratios (#85)
- Reference sequence similarity QC (#84)
- Squiggy session JSON export for Positron IDE
- Utility to collapse redundant GtRNAdb FASTA sequences
- Multiple 3' adapter support for PT tag detection
- Skip mode for reference validation
- Pre-download dorado mod base models rule (avoids race conditions)
- nvitop GPU monitoring dependency

### Changed
- `classify_charging` switched from GPU to CPU with parallel workers (8 threads)
- WarpDemuX workflow simplified: eliminated `merge_pods_for_demux`, passes raw POD5 dirs directly
- `bwa_align` filtering changed from `-F 4` to `-F 20` (also excludes reverse-strand reads)
- Removed redundant awk position filter from `bwa_align`
- Removed `protected()` directive from `rebasecall` output

### Fixed
- Race condition when parallel GPU jobs download dorado modification models simultaneously
- Reverse-strand reads not filtered at alignment step
- Redundant awk position filter in `bwa_align` superseded by adapter-based filtering
- Graceful fallback for `get_pipeline_commit` when git unavailable
- Various snakefmt formatting and test corrections

## [v0.1.0] - 2025-01-16

### Added
- Run manifest generation for reproducibility tracking
- Native SLURM cluster support with GPU configuration
- Pytest unit tests and expanded CI coverage (#80)
- Adapter position tagging with parasail alignment - adds PT tags (#79)
- MkDocs documentation website (#77)
- Publication citation (White et al. 2025 Nat Commun)
- Reference validation and building step to ensure tRNA sequences have proper CCA endings and adapter structure required for charging classification
- Optional WarpDemuX barcode demultiplexing support for pooled/multiplexed sequencing runs (#74)
- Optimized modkit thresholds from ModkitOpt
- Pixi package manager support as primary environment manager
- Mermaid diagram for workflow visualization

### Changed
- Standardized output directory structure to nested sample paths
- Updated dorado version from 0.9.1 to 1.3.1
- Migrated GitHub Actions CI from conda to pixi (#75)
- Separated tool installation from environment activation (`pixi run setup`)
- Updated README to use pixi instead of conda (#72)
- Normalized file permissions across repository

### Fixed
- **Major fix**: Improved 5' adapter detection from 0.04% to 82% (#82)
- Resolved GLIBCXX version errors with configurable CUDA support
- Shell glob pattern errors in tests

## 2025-11-07

### Added
- Claude Code session start hook for automated development setup (#69)
- Comprehensive CI/CD build and test checks (#68)
- New project initialization structure (#65)

### Changed
- Prefer pandas over polars for stability on some cluster nodes
- Reduced Remora logging level

## 2025-06-22

### Added
- LSF-specific cluster configuration (#63)

### Changed
- Renamed test directory to recommended `.tests/` location
- Reduced and made optional dorado verbosity
- Downgraded numpy to fix remora stats compatibility

## 2025-03-16

### Added
- Modkit integration for RNA modification analysis (#59)
  - `modkit_pileup` rule for modification pileups
  - `modkit_summary` rule for modification summaries
  - `modkit_extract` and `modkit_extract_full` rules for detailed modification data
- Automatic dorado and model download/installation (#56)
- Modified base calling support (pseU, m5C, inosine_m6A)
- Full modkit outputs with optimized memory allocation

### Changed
- Eliminated support for FAST5 files - pipeline now POD5-only (#43)
- Reorganized output directory structure
- Renamed charging tags during transfer (ML→CL, MM→CM)
- Updated model download strategy
- Increased memory allocation for modkit rules

### Fixed
- Restored `-v` option in dorado for proper verbosity control

## 2025-01-08

### Added
- Rule for calculating CPM of charged/uncharged tRNAs (#28)
- Remora CCA classifier for charging state classification (#18)
- GPU pipeline support for `cca_classify` rule (#23)
- Charging probability extraction and analysis (#46)

### Changed
- Implemented ML threshold (≥200 = charged, <200 = uncharged)
- Compressed output files for storage efficiency
- Various tweaks to file handling (#27)

### Fixed
- Actually use the threshold value in classification

## 2024-08-13

### Added
- Alignment filtering capabilities with configurable parameters (#13)
- Optional Remora signal metrics extraction
- Kmer models included in pipeline
- Logging and optional failed BAM outputs
- New BAM tag indicating why reads are filtered
- Support for processing reads from both pass and fail directories

### Changed
- Cleanup filtering approach for full-length tRNA reads
- Updated test data
- Ignore supplementary and secondary alignments
- Snakemake v8 compatibility (#10)

### Fixed
- Insertion double-counting bug (#15)
- Dropped redundant summary align stats (#14)
- Added 'pod5' to list of possible pod5 directories

## 2024-05-19

### Added
- Support for merging multiple sequencing runs per sample
- Support for unmapped BAM as input (#5)
- Pipeline commit and config recording for reproducibility (#9)
- Bedgraph output generation
- Alignment statistics calculations
- Base calling error frequency calculations

### Changed
- Use v5.0.0 dorado models with modification calling
- Expose dorado and bwa command-line options
- Reworked alignment stats output (#8)
- Set rebasecalled outputs as read-only

### Fixed
- Keep additional BAM flags (e.g., pi) during processing (#1)
- Use -T → -C options to preserve all BAM tags from dorado

## 2024-02-07

### Added
- Initial pipeline release
- Core workflow: POD5 merge → rebasecall → align → filter
- BWA MEM alignment to tRNA + adapter reference
- Post-alignment filtering for full-length tRNAs
- Basic summary statistics generation
- Snakemake workflow with modular rule structure
- Conda environment specification
- Sample configuration via TSV files
