# Changelog

All notable changes to the aa-tRNA-seq pipeline are documented in this file.

## [Unreleased]

### Changed
- **The charging classifier is now `escpod signal classify` against a vendored ONNX model bundle, and Remora is retired.** The pipeline no longer installs, imports or invokes Remora anywhere.
  - The model is `charging_feature_nn_rna004@v0.1.0`, vendored at `resources/models/charging/` (upstream `rnabioco/escapepod-models` is private and compute nodes have no route to GitHub). It is a bundle **directory** and is self-describing: the anchor, feature recipe, k-mer table (pinned by sha256), abstain rule and recommended operating point all come from its `metadata.json`, not from flags. Published test AUROC 0.9906, balanced accuracy 0.9625 on a held-out flowcell.
  - The bundle's `9mer_levels_v1.txt` is a symlink to `resources/kmers/9mer_levels_v1.txt`, which is byte-identical to the released copy — the file the model was trained against. Replacing that file invalidates the model; the sha256 in `metadata.json` gates it.
  - New config block `charging` with `model`, `min_mapq` and `ml_threshold`. `remora_cca_classifier`, `remora_kmer_table`, `classifier` and `opts.remora` are gone.
  - `charging.min_mapq` defaults to **0**, not escpod's own default of 1. tRNA references are highly redundant, so a read mapping equally well to two isodecoders gets MAPQ 0 from bwa and is still a good read; `--min-mapq 1` drops 118 of 209 records (56%) on the test data, while 0 classifies exactly the 187 reads Remora did.
  - `charging.ml_threshold` replaces the value hardcoded in `get_cca_trna_cpm`. It is 200, matching the bundle's declared operating point.
  - `classify_charging` is CPU-only; `escpod signal classify` has no GPU path. GPU entries removed from the cluster profiles.
- **`transfer_bam_tags` is removed, and there is no `cm` tag.** Remora emitted its score into `MM`/`ML` — the standard modbase tags — clobbering the calls modkit needs, so a rename step existed purely to move them to `cl`/`cm`. `escpod signal classify` writes `cl` (uint8, `round(P(charged) * 255)`) directly onto the input records, in the same order, leaving dorado's `MM`/`ML` intact. Nothing downstream consumed `cm`.
- **The `classifier: remora|leech` switch and the `classify_charging_leech` rule are removed.** They selected between two implementations of the same retired `cca_classifier.pt`, which is also deleted. leech's amino-acid identity rules (`classify_aa`, `aa_identity`) are unaffected.
- **`escpod` bumped 0.8.1 -> 0.11.0**, pinned alongside the model: the bundle is the per-base-feature ONNX variant, which escpod reads only from 0.10.0 onward (older binaries refuse it with ``missing field `gbm` ``). 0.11.0 additionally applies the bundle's abstain rule (previously parsed and dropped, so ~0.85% of scoreable reads got a confident call they should not have) and adds the `--tsv` `reason` column.
- **`escpod` is now checked unconditionally at startup**, not only for LDX runs — `merge_pods` and `classify_charging` use it on every run. The `pod5` CLI check moved under the warpdemux conditional, which is the only remaining consumer of that package.
- `pixi run setup` no longer installs PyTorch. Nothing in the default pipeline needs it; leech's `.pt` bundles still do, so enabling `classify_aa`/`aa_identity` now means installing torch yourself.
- Updated dorado from 1.4.0 to 2.1.1. The basecalling model is unchanged: `rna004_130bps_sup@v5.3.0` remains the newest RNA004 sup model in dorado 2.1.1, and all four modified-base models the pipeline uses (`m5C_2OmeC`, `inosine_m6A_2OmeA`, `pseU_2OmeU`, `2OmeG`) are still available at v5.3.0.
- Updated dorado from 1.3.1 to 1.4.0, basecalling model from v5.1.0 to v5.3.0, and modkit from 0.6.0 to >=0.6.1.
- Shell scripts (`setup-env.sh`, `setup-tools.sh`) now read `dorado_version` and `dorado_model` from `config/config-base.yml` instead of hardcoding defaults.
- Updated stale documentation references for dorado and modkit paths/versions.
- Pipeline now warns on startup if installed dorado version does not match `dorado_version` in config.

### Added
- **`summary/tables/{sample}/{sample}.charging_calls.tsv.gz`** — one row per read the charging model saw, with `p_charged`, `cl`, and a `reason` naming why every unscored read was unscored (`no_aligned_arm`, `no_signal`, `ns_mismatch`). Not optional: the model **abstains** on reads whose common arm did not align rather than guessing, abstention is charging-correlated, and so a charging fraction over called reads alone is an **underestimate**. Report the no-call rate beside it.
- `read_attrition` folds those reasons in, so the `aligned -> charge-called` gate names its loss instead of inferring it. `anchor_coverage` still covers the complementary half — reads that never reached the model at all.
- `pixi run verify-charging-model` verifies the vendored bundle against its pinned checksums.
- `trim_reference` rule produces a tRNA-only FASTA (`trna_only.fa`) by stripping 5'/3' adapter sequences from the adapted reference. This FASTA is used by clover for MODOMICS annotation and structure visualization.
- `build_trna_reference.py --mode trim` for generating adapter-stripped tRNA-only FASTA files.
- `get_bcerror_freqs.py` and `compute_odds_ratios.py` accept `--offset-5p` and `--offset-3p` to filter adapter positions and convert to tRNA-only coordinates.

### Fixed
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
