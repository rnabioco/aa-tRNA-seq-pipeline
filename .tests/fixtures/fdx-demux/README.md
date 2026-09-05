# Dual-index (LDX + FDX) demultiplexing test fixture

365 reads of a real dual-indexed nanopore run: every molecule carries a 3' LDX
code and a 5' FDX code. The second barcoded input in this repository, and the
one that lets the FDX axis be tested end to end.

## Provenance

Donor: Run2 of the 20260828 FDX pilot (`20260828_1246_P2S-01617-A_PBK76651_a6de399f`,
flowcell PBK76651) — acylated *E. coli* MG1655 tRNA, one SQK-RNA004 flowcell,
2,471,694 reads, five libraries each carrying one 5' FDX code and three 3' LDX
codes (`config/fdx4_libraries.yaml` in escapepod-models):

| library | fdx | ldx codes | 3' adapter |
|---|---|---|---|
| A | fdx01 | ldx01 ldx02 ldx03 | plain |
| B | fdx02 | ldx04 ldx05 ldx06 | plain |
| C | fdx03 | ldx07 ldx08 ldx09 | plain |
| D | fdx04 | ldx10 ldx11 ldx12 | plain |
| E | fdx01 | ldx13 ldx14 ldx15 | edx07 |

Two properties make this run the right donor:

- **The fdx label is independent of the fdx signal.** It comes from the LDX
  call (3' end, `barcode_crf_nbc16_rna004@v0.3.0`) mapped through the pool
  design, not from any 5' model. So the FDX axis can be scored against it.
- **Run2 is held out.** The shipped `barcode_crf_fdx4_rna004@v0.2.0` was trained
  on Run1 only (`manifest_runs: [Run1]`), so what the FDX call does on these
  reads is a cross-flowcell number, not a fit to training data.

Reads were selected from the donor's own analysis (`results/fdx/demux/
Run2_libraries.parquet`, and `scratch/fdx/Run2.bam` for alignment status), so
each read's codes and alignment status were known before it was chosen. POD5
records are copied byte-for-byte out of the instrument POD5 by `escpod filter`.

Rebuild with `.tests/make_fdx_test_data.sh` (it documents the exact source paths;
run it in its own allocation, it scans 53 GB of POD5).

## Composition

| reads | population | selected as | covers |
|------:|---|---|---|
| 100 | ldx01 + fdx01 (library A), aligned | `ldx01:fdx01:aligned` | sample `fdx01_ldx01` |
| 100 | ldx04 + fdx02 (library B), aligned | `ldx04:fdx02:aligned` | sample `fdx02_ldx04` |
| 100 | ldx07 + fdx03 (library C), aligned | `ldx07:fdx03:aligned` | sample `fdx03_ldx07` |
| 40 | ldx10 + fdx04 (library D) | `ldx10:fdx04:unclaimed` | a pair no sample claims; demux must route it nowhere |
| 25 | `unclassified` on the LDX axis | `unclassified` | reads the CRF could not call |

`fixture_manifest.tsv` records which population every read came from. Only
plain-3'-adapter libraries are included, so `config-fdx-test.yml` configures the
plain adapter as a single string and no `edx:` keys; the EDX filter is covered by
the LDX fixture.

## What the outputs measure

`demux/read_ids/run/demux_summary.tsv.gz` is the LDX axis, `demux/read_ids/run/
fdx/demux_summary.tsv.gz` the FDX axis, and `demux/read_ids/run/assigned_summary.tsv`
the join — how many of each sample's 100 selected reads survived both calls.
The LDX call reproduces the fixture's own routing up to retrain variance
(the donor used nbc16 v0.3.0, the pipeline ships ldx16 v0.1.0); the FDX call is
the held-out measurement described above. `tests/integration/test_fdx_demux.py`
asserts the shape of both, not the decimals.

## Reference

`collapsed.fa`: the 47 collapsed *E. coli* K-12 MG1655 mature tRNA sequences
(`eschColi_K_12_MG1655-mature-tRNAs-collapsed.fa` in escapepod-models), raw,
every one ending in CCA. The pipeline builds the adapted reference from it
(`reference.mode: build`) with the plain 3' adapter,
`GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA` — the EDX scaffold truncated before the
barcode, recovered from the donor run's soft-clipped tails (escapepod-models
`dev-notes/fdx-three-prime-chemistry.md`) because no design document recorded it.
