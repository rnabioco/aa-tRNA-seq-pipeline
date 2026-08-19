# LDX demultiplexing test fixture

415 reads of a real pooled, barcoded nanopore run — the only barcoded input in
this repository, and the reason the demux path can be tested end to end at all.

## Why this exists

`config/config-demux-test.yml` could never complete (issue #120). It layered
fabricated WDX barcode assignments and `edx:` keys on top of `.tests/sample1`,
which is unbarcoded sacCer3 data: every read detected as adapter `none`, so
`extract_edx_read_ids` wrote an empty list and `escpod filter` failed — correctly,
since a barcode matching zero reads should fail loudly. Dropping the `edx:` keys
did not rescue it either, because bwa maps **none** of the reads WarpDemuX happens
to route to those barcodes.

No amount of config editing fixes that. The data has to actually be barcoded.

## Provenance

Donor: `20260806_1113_P2S-00519-A_PBG58575_3669d73d`
(`20260806_EDC_activation_test_EDXs/EDC_act_test`) — one SQK-RNA004 flowcell
carrying all 16 LDX barcodes, 1,001,307 classified reads. Reads were selected
from that run's own completed pipeline outputs, so each read's barcode, 3'
adapter and alignment status were all known *before* it was chosen. POD5 records
are copied byte-for-byte out of the instrument POD5 by `pod5 filter`; nothing is
synthesised or re-encoded.

Rebuild with `.tests/make_ldx_test_data.sh` (it documents the exact source paths).

## Composition

| reads | population | selected as | covers |
|------:|---|---|---|
| 100 | `nbc01`, adapter `edx01`, aligned | `nbc01:edx01:on_target` | survives the EDX filter, aligns, gets classified |
| 25 | `nbc01`, other adapters, aligned | `nbc01:other_adapter` | **removed** by the EDX filter |
| 100 | `nbc02`, adapter `edx02`, aligned | `nbc02:edx02:on_target` | as above, second barcode |
| 25 | `nbc02`, other adapters, aligned | `nbc02:other_adapter` | **removed** by the EDX filter |
| 100 | `nbc08`, mixed adapters, aligned | `nbc08:pool_aligned` | the *unfiltered* path — no `edx:` key |
| 40 | `nbc05` | `nbc05:unclaimed` | a barcode no sample claims; demux must route it nowhere |
| 25 | `unclassified` | `unclassified` | reads the CRF could not call |

`fixture_manifest.tsv` records which population every read id came from, so a
test can assert on the expected outcome per read rather than on totals alone.

The off-target reads matter: without them `filter_{fastq,pod5}_by_edx` would be
a no-op that still passes, which is a test that cannot fail.

## Reference

`collapsed.fa` — 47 human tRNA sequences, raw (no adapters). The donor run is
human, not sacCer3, so this fixture ships its own reference and
`config-ldx-test.yml` sets `reference.mode: build` to adapt it. 4.6 KB.

## Why it is committed rather than downloaded

The rest of the test data lives in an S3 tarball fetched by
`pixi run dl-test-data`. This fixture is 4.6 MB, versioned alongside the configs
that reference it, so a checkout can run the demux path without a download step
and without the fixture drifting out of sync with the code. `.gitignore` carries
two narrow negations for it — the repo otherwise ignores `.tests` and `*.pod5`.
