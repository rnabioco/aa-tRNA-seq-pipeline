# Vendored demux models

Barcode/boundary models for `escpod demux`, committed here rather than fetched.

## Why vendored and not `escpod demux models fetch`

`escpod` ships a pinned manifest with a verified cache
(`escpod demux models fetch <bundle>`), and that is the right way to obtain a
model when the manifest knows about it. As of escapepod-rs v0.7.0 the manifest
carries exactly one bundle, `wdx4_rna004`. The CRF NBC16 bundle below is not in
it, so there is no `fetch` name that resolves to it.

The upstream repository (`rnabioco/escapepod-models`) is also private, so
`fetch` needs a `$GITHUB_TOKEN` even for the bundles it does know — which
compute nodes on this cluster generally cannot use, since they have no route to
GitHub. Committing the bundle removes the token, the network, and the cache
from the run-time path entirely.

## Vendored bundles are byte-for-byte upstream

**Do not edit a bundle, and do not add files to one.** Not `metadata.json`, not
a regenerated `SHA256SUMS.txt`. Verify with:

```bash
pixi run verify-demux-model
```

which checks each bundle against the hashes UPSTREAM declares
(`metadata.boundary.sha256`, `provenance.sha256`) plus any `SHA256SUMS.txt` the
release itself shipped — never a checksum file written here.

This rule exists because it was broken. PRs #107 and #108 patched
`boundary.margin: 0` and `boundary.clamp_max_shift: 300` into the nbc16
`metadata.json`. `SHA256SUMS.txt` was not regenerated, so that bundle failed its
own integrity check for eleven days, and the patched copy then read as evidence
that a later upstream release had *removed* those keys — it had not; **no
upstream CRF bundle has ever declared either key**, because
`build_crf_bundle.py` cannot write them. Both values now live in
`config-base.yml` under `ldx:`, where they are visible and diffable, and the
measurements that justify them are recorded below.

## `barcode_crf_ldx16_rna004@v0.1.0`  (default)

16-plex CTC-CRF barcode basecaller, and **the first bundle whose class names are
`ldx01`..`ldx16`** rather than upstream's older `nbc` vocabulary. Released
2026-08-19 from
`https://github.com/rnabioco/escapepod-models/releases/tag/barcode_crf_ldx16_rna004%40v0.1.0`.

Successor to the `barcode_crf_nbc16_rna004` family, which is **closed at
v0.3.1**. It is a retrain, not just a rename: corrected geometry
(escapepod-models#36, `state_len=4` so the full 27-nt code is emitted rather
than discriminating on 23 nt) and the first trained without bonito or ont-koi,
which upstream established as equivalent over 3 paired seeds.

Because it is a retrain, **its calls are not identical to nbc16's**: on the
415-read LDX fixture the two disagree on 34 reads (8.2%), some into barcodes the
fixture contains none of. Neither is ground truth there — the fixture was built
from nbc16's routing — but demux results are not comparable across the switch.

```bash
pixi run escpod-model-info      # geometry, references, published metrics
```

The zip ships no `SHA256SUMS.txt`; both ONNX graphs are instead pinned by hashes
inside the bundle's own metadata, which `verify-demux-model` checks.

## `barcode_crf_wdx4_rna004@v0.2.0`  (WDX4 panel, provisional)

4-plex CTC-CRF over the **WarpDemuX** panel — `bc03`, `bc04`, `bc05`, `bc07`,
which this project calls `barcode03`..`barcode07`. It exists so a WDX-barcoded
run can be demuxed by escpod on the sidecar path instead of by WarpDemuX
routing reads into a second full copy of the POD5.

Select it per run; it is not the default:

```yaml
ldx:
  enabled: true
  model: resources/models/demux/barcode_crf_wdx4_rna004@v0.2.0
  min_crf_margin: 2.0   # see below -- provisional
```

Nothing else changes. The bundle carries its own references, so no
`--barcodes` and no `--method` are passed, and the barcode-name crosswalk
(`workflow/scripts/barcode_names.py`) maps a sample configured `barcode03` onto
the `bc03` the bundle emits.

**The panel is four codes. `barcode11` is not in it** and cannot be added by
configuration — upstream's note is that bc11 and the other seven WarpDemuX
codes "need new sequencing, not a config change". A WDX4b run using bc11 has no
escpod equivalent today.

### Provisional, on two counts

**`min_crf_margin` is not yet a measured operating point.** Swept on a
1.7M-read run, `>= 2.0` gives 99.87% agreement with WarpDemuX's confident calls
at 0.113% worst-barcode cross-contamination, while keeping 1,236,235 calls
against WarpDemuX's own 1,107,936 — better on both axes at once. But **that run
is in this model's training corpus**: the v0.2.0 split is a random read-level
10%, not leave-one-run-out, and every WDX4 run in `2026-aars-in-vitro` is in it.
So the figure is an upper bound. Re-derive it with
`scripts/crf_margin_sweep.py` against a genuinely held-out run before treating
it as settled.

**The labels it learned are ungated.** Upstream's config sets `no_gate: true`
and says so plainly — "THE GATE IS OFF, DELIBERATELY, AND THIS IS A REAL COST".
WarpDemuX declines ~15% of reads and its balanced accuracy on those is ~0.68,
so v0.2.0 inherits noise that v0.1.0's confidence-gated corpus did not. A gated
retrain, held out by run, is the fix; expect a v0.3.0 to supersede this bundle.

Measured behaviour worth knowing before reading a concordance report: **every
read this model leaves `unclassified` is one it never decoded** — the adapter
ended before the 2000-sample window — not one it was unsure about. That is
120,158 reads (7%) on a 1.7M-read run, and it is a `boundary_margin` /
`clamp_max_shift` question rather than a classifier one.

## `barcode_crf_nbc16_rna004@v0.2.0`  (retained)

Kept so runs pinned to it stay reproducible. Not the default.


16-plex CTC-CRF barcode basecaller for the LDX adapters (upstream calls these
barcodes `nbc01`..`nbc16`; **LDX is the name we use for them**). Released
2026-08-01 from
`https://github.com/rnabioco/escapepod-models/releases/tag/barcode_crf_nbc16_rna004%40v0.2.0`.

The bundle is self-describing: `metadata.json` carries the 16 barcode
references, the signal geometry, the standardisation constants, and the
boundary detector it is pinned to. That is why the run-time invocation needs no
`--barcodes` and no `--method`:

```bash
escpod demux <pod5>... --model resources/models/demux/barcode_crf_nbc16_rna004@v0.2.0 -d out/
escpod demux --model resources/models/demux/barcode_crf_nbc16_rna004@v0.2.0 --info
```

### Contents

| File | sha256 | Notes |
|---|---|---|
| `barcode_crf_nbc16_rna004.onnx` | `5b626e5e…3134c6` | CRF encoder. Decode is *not* in the graph — standard ONNX ops cannot express it, so `escpod` runs the Viterbi over 256 states × 5 transitions itself. |
| `adapter_rna004.onnx` | `b59f8667…3d26b5` | Boundary CNN, `adapter_rna004@v1.1.0`. |
| `metadata.json` | — | Runtime sidecar: references, geometry, standardisation, boundary pin. |
| `provenance.json` | — | Training provenance and published metrics. |

Both checksums were verified against two independent sources: the release's own
`provenance.json`/release notes, and — for `adapter_rna004` — the pinned member
hash in escapepod-rs's `demux/models.rs` manifest. Re-check with
`sha256sum -c SHA256SUMS.txt` from inside the bundle directory.

### The boundary model is not interchangeable

`adapter_rna004.onnx` ships *inside* the bundle because the CRF's training
window is defined by that detector's `adapter_end`. Substituting the LLR
detector costs 17.2 points of balanced recall, and the failure is silent — it
runs and produces plausible output. `escpod` refuses `--method llr` against a
bundle pinned to `cnn`, so do not pass `--method` at all.

### Why `ldx.boundary_margin: 0` and `ldx.clamp_max_shift: 300`

Both are set in `config-base.yml` and passed to `escpod` as flags. No bundle
declares them, so leaving them unset does not mean "the model decides" — it
means escpod's fallback of margin 200 and no clamp. On the 415-read fixture that
fallback costs 24 reads (390/415 -> 366/415).

The measurements below are what justify the two values. They were made on the
nbc16 model; the geometry they describe (`chunk` 3000) is unchanged in ldx16.

`margin` is the samples of `adapter_end` a read needs *beyond* `signal.chunk`
before the CRF will decode it. Absent, escpod falls back to 200 — the filter
`extract_chunks.py` applied when building the training corpus, so reads below it
were never represented. But that describes how the corpus was *selected*, not
what the encoder requires, which is a full `chunk` of history and no more. Reads
in `[chunk, chunk + 200)` were being routed to `unclassified` undecoded, with
confidence 0.

Measured on the 2026-08-06 nbc16 run (1,001,307 reads, 145,775 unclassified):
every unclassified read had a detected adapter, and 98.2% failed only this gate.
Decoding the affected band at margin 0 returns **36,921 reads, 100% of the band,
at median edit distance 0 with 98.3% within 2 edits** — against references whose
minimum pairwise distance is 12, and cleaner than the reads that already passed
(96.4%). Calls spread across all 16 barcodes, and 84.5% align to tRNA. Demux
yield 85.44% -> 89.13%. Reads below `chunk` still decode 0%, confirming the
window genuinely does not exist there.

`escpod demux --boundary-margin` documents the bundle as the intended long-term
home for this value (rnabioco/escapepod-rs#193), and that remains right — but it
is an UPSTREAM fix, in `build_crf_bundle.py`, not something to patch in after
the fact. Until a release declares it, the config is the honest place: it
survives re-fetching a bundle, and it shows up in a diff.

#### `clamp_max_shift: 300`

`margin` cannot reach a read whose adapter ends before `chunk` (3000): its window
would start before sample 0, so there is nothing to relax. `clamp_max_shift`
instead keeps the window width and anchors it at the read start — `[0, chunk]` —
sliding `chunk - adapter_end` samples of downstream signal into the tail, for
reads within the bound.

The model tolerates that slide well. Known-good reads deliberately slid forward
still call the same barcode 98.6% of the time at shift 0 and **93.5% at shift
500**, so the ceiling is set by the data, not the decoder. Applied to the real
`adapter_end` 2,500-2,999 band, every read decodes at median edit distance 0, but
two things decay together across it:

| shift | within 2 edits | still aligns to a tRNA |
|---|---|---|
| 0-99 | 97.4% | 78.5% |
| 100-199 | 96.3% | 72.6% |
| 200-299 | 95.6% | 66.0% |
| 300-399 | 94.9% | 60.4% |
| 400-499 | 92.9% | 50.0% |

**300 is a judgement, not a measurement**: it keeps agreement above ~95% and the
aligning fraction near two thirds, and gives up the ~5,800 reads past it where
half no longer align. The alignment decay is a property of these reads — a bigger
shift means more of the adapter was truncated to begin with — not of the model,
so raise `ldx.clamp_max_shift` per run if an analysis wants the tail.

### Published accuracy

Exact match to the emitted reference is 0.9736; balanced precision/recall at a
0.97 recovery threshold is 0.9858/0.9565. Note the training pilot
(`20260728_LDX_Demux_Pilot`) is **confounded** — each LDX sat on exactly one of
two flowcells — so upstream reports two 8×8 matrices rather than one 16×16.
That caveat is about the *evaluation*, not about applying the model.

### Matching against references

Edit distances are computed against the sequences the model **emits**
(`target[state_len:]`, i.e. 44 nt), not the 48-nt training targets. The bundle
already carries the emitted form, which is the reason to prefer it over a
hand-written `--barcodes` CSV: full-length targets still call the same barcode
but inflate every distance and compress the confidence margin that
`--min-margin` gates on.
