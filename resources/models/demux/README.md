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

## `barcode_crf_nbc16_rna004@v0.2.0`

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

### Local amendment: `boundary.margin: 0`

**`metadata.json` in this copy differs from the released v0.2.0 by one key.** The
`boundary` block gains `"margin": 0`. Nothing else is touched — the ONNX graphs,
references, standardisation and geometry are byte-identical to upstream, and the
sidecar is not covered by any checksum (the pinned adapter is, via
`boundary.sha256`, which v0.2.0 does not declare).

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

Declaring it here rather than in the run config keeps the value with the model,
which is what `escpod demux --boundary-margin` documents as the intended home
(rnabioco/escapepod-rs#193). **Drop this amendment** when escapepod-models
re-exports the bundle with the key set upstream; until then, re-fetching the
released v0.2.0 silently reverts to 200.

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
