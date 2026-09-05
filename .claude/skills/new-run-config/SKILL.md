---
name: new-run-config
description: Set up the config pair (config-<name>.yml + samples file) for a new aa-tRNA-seq sequencing run, validate it against the vendored demux bundles and the run directory, and dry-run the DAG. Use when a user has a new run to process, names a POD5 directory, lists barcodes or libraries, or asks how to write a config or samples file.
---

# Configure a new run

A run needs two files under `config/`: `config-<name>.yml`, holding only what
differs from `config/config-base.yml`, and `samples-<name>.yml` (or `.tsv` for
unbarcoded runs), one entry per sample naming the barcode(s) its library
carries. `workflow/scripts/run_config.py` writes both from a description of the
libraries and refuses a plan the pipeline would refuse. Prefer it over writing
YAML by hand; the rules it enforces are the ones that otherwise surface hours
later as "no reads were assigned" (#120 is the canonical case).

## 1. Gather

Ask for, or find in the run's own files, before writing anything:

| what | where it usually is | why it matters |
|---|---|---|
| run directory | the MinKNOW output dir holding `pod5_pass/`, `pod5_fail/` or `pod5/` | `find_raw_inputs` globs exactly those three |
| library design | the wet-lab sheet: which sample is which barcode | one sample per barcode tuple; replicates are separate samples |
| barcode axes | 3' **LDX** code (`ldx01..ldx16`), 5' **FDX** code (`fdx01..fdx04`, dual-index libraries only), 3' adapter **EDX** identity (`edx01..`) | each axis is a different demux mechanism; FDX needs LDX |
| 3' adapter(s) actually ligated | the library prep; recover from soft-clipped tails if undocumented (see escapepod-models `dev-notes/fdx-three-prime-chemistry.md`) | must match `adapters.three_prime` by name AND sequence; the base config's `edx01/edx02` are the sacCer3 dual-adapter set, not necessarily yours |
| reference | raw mature tRNA FASTA for the organism (build mode) or an already-adapted one | `resources/ref/` and escapepod-models `resources/ref/` hold the ones used so far |
| GPU for demux | is a GPU build of escpod installed (`pixi run setup`)? | `ldx.gpu: true` is a 20x speedup and fails loudly without one |

Sample spec grammar for `--sample NAME=SPEC`:

```
ldx01               3' LDX code
ldx01+fdx01         dual index: 3' LDX + 5' FDX
ldx01/edx01         LDX code, keep only reads with 3' adapter edx01
ldx01+fdx01/edx07   all three
wdx:barcode03       WarpDemuX barcode (retired backend)
-                   unbarcoded: the whole run is this sample
```

## 2. Write

```bash
pixi run new-run-config -- --name <name> \
    --output-dir results/<name> \
    --run /path/to/run_dir \
    --sample wt_rep1=ldx01+fdx01 --sample wt_rep2=ldx02+fdx01 \
    --sample ko_rep1=ldx04+fdx02 \
    --reference-raw resources/ref/<organism>-mature-tRNAs-collapsed.fa \
    --three-prime plain=GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA \
    --dry-run
```

`--sample` entries belong to the `--run` before them; give several runs in
sequence. Omit `--three-prime` only when the libraries really used the base
config's adapters. `--json` gives a machine-readable report. Existing files are
never overwritten without `--force`.

The checks, and what a failure means:

- **code not emitted by the bundle**: a typo, or the wrong panel (`nbc` names
  are retired; the shipped bundles emit `ldx`/`fdx` names)
- **fdx without ldx / mixed axis sets / identical tuples**: the join in
  `select_demux_reads.py` would swallow or duplicate reads
- **edx not a declared adapter**: the `edx:` name must be in
  `adapters.three_prime`
- **3' adapters differ in length, or do not start with GGC**: the reference
  builder needs one length; the charging model anchors on the `CCA|GGC` junction
- **run directory has no pod5 dirs / reference missing**: the DAG would be
  empty or fail at `build_reference`

## 3. Check and dry-run

```bash
pixi run check-run-config -- config/config-<name>.yml --dry-run
```

re-applies every check to an existing pair and prints the job table. A dry-run
proves the DAG, not the data: it never opens a POD5, so an adapter set that
matches nothing still passes. When the run is barcoded and small enough, the
real proof is the first sample's `demux/read_ids/<run>/demux_summary.tsv.gz`.

## 4. Run

The controller is itself a Slurm job (see the `hpc-compute` skill): submit it
with `sbatch` on a CPU partition with a long wall, and let it drive the rules
through `--profile=cluster/slurm`. Watch `summary/read_attrition.tsv.gz` when
it finishes; every gate's loss is one row there.

## Do not

- edit `config/config-base.yml` for a project: that file is the shared
  default, and a change there moves every run
- edit a vendored model bundle under `resources/models/`; overrides belong in
  config (`resources/models/demux/README.md` records why)
- copy POD5 into the repository or a worktree; configs point at the canonical
  run directory
