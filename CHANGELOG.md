# Changelog

All notable changes to the aa-tRNA-seq pipeline are documented in this file.

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
