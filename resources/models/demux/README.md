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

Both ONNX checksums were verified against two independent sources: the release's
own `provenance.json`/release notes, and — for `adapter_rna004` — the pinned
member hash in escapepod-rs's `demux/models.rs` manifest. Re-check with
`sha256sum -c SHA256SUMS.txt` from inside the bundle directory, or
`pixi run verify-demux-models` for every bundle at once.

Note that the `metadata.json` hash here is of the **locally amended** file (see
below), so it attests to our copy, not to upstream's. It went stale once
already: the `clamp_max_shift` amendment edited `metadata.json` without
regenerating this file, and the mismatch went unnoticed until
`verify-demux-models` existed. **Regenerate `SHA256SUMS.txt` whenever a sidecar
is amended.**

### The boundary model is not interchangeable

`adapter_rna004.onnx` ships *inside* the bundle because the CRF's training
window is defined by that detector's `adapter_end`. Substituting the LLR
detector costs 17.2 points of balanced recall, and the failure is silent — it
runs and produces plausible output. `escpod` refuses `--method llr` against a
bundle pinned to `cnn`, so do not pass `--method` at all.

### Local amendment: `boundary.margin: 0` and `boundary.clamp_max_shift: 300`

**`metadata.json` in this copy differs from the released v0.2.0 by two keys**, both
in the `boundary` block. Nothing else is touched — the ONNX graphs,
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
so raise the bound per run (`ldx.boundary_margin`'s sibling, or
`--clamp-max-shift`) if an analysis wants the tail.

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

## `barcode_crf_wdx4_rna004@v0.2.0`

4-plex CTC-CRF barcode basecaller for the **WarpDemuX** panel — the kit this
pipeline documents as `WDX4_tRNA_rna004_v1_0`. Upstream names these codes
`bc03`/`bc04`/`bc05`/`bc07`; the pipeline's project-facing names are
`barcode03`/`barcode04`/`barcode05`/`barcode07`, and `barcode_names.py`
translates between them (the same job it does for `nbc`->`ldx`, in the opposite
direction — here the *configured* name is ours and the *emitted* name is
upstream's).

It exists to replace WarpDemuX itself. `escpod demux` already runs the LDX
bundle above through one fused pass; pointing it at this bundle puts the WDX
panel on the same path and removes WarpDemuX from the dependency graph.

### Provenance is a worktree commit, not a release

Unlike the nbc16 bundle, **there is no release to fetch**. `models/*.onnx` is
gitignored in escapepod-models, no zip exists in `dist/`, and escpod's built-in
manifest has no name that resolves to it. The only complete copy upstream is:

```
escapepod-models/.claude/worktrees/wdx-crf/models/barcode_crf_wdx4_rna004@v0.2.0/
  (branch wdx-crf, commit b9b9fd8
   "feat(registry): ship barcode_crf_wdx4_rna004@v0.2.0, and separate the training seed (#68)")
```

Do not assume `escpod demux models fetch` reproduces these bytes — it cannot.
Re-vendor from that worktree, or from a release once one is tagged (see the
licence note below for why one has not been).

### Contents

| File | sha256 | Verifiable against |
|---|---|---|
| `barcode_crf_wdx4_rna004.onnx` | `a33459d1…62b93` | the bundle's own `provenance.json` **and** escapepod-models `MANIFEST.json` — two independent sources |
| `adapter_rna004.onnx` | `b59f8667…3d26b5` | this bundle's `boundary.sha256`, the pinned member hash in escapepod-rs's `demux/models.rs`, and byte-for-byte equality with the nbc16 bundle's copy |
| `metadata.json` | `4986ceb5…82af2` | **local integrity only** — nothing upstream publishes a hash for the sidecar |
| `provenance.json` | `b46e6ef0…e9820` | **local integrity only** |

Re-check with `sha256sum -c SHA256SUMS.txt` from inside the bundle directory.

Note what the last two rows are worth. A `SHA256SUMS.txt` entry for a sidecar
proves only that nobody edited it since vendoring — and if a local amendment is
ever applied (as one was for nbc16), the recorded hash becomes *ours*, not
upstream's. The same caveat applies to the nbc16 `SHA256SUMS.txt` above, whose
`metadata.json` hash is of the amended file.

`adapter_rna004.onnx` is byte-identical to the copy in the nbc16 bundle, and it
is **deliberately duplicated rather than shared or symlinked**: escpod resolves
`boundary.onnx` relative to the bundle directory, and a bundle that reaches
outside itself stops being self-contained and independently hash-verifiable.
512 KB is the price of that property. Please do not "fix" it.

### Geometry differs from the nbc16 bundle

|  | nbc16 | **wdx4** |
|---|---|---|
| `signal.chunk` | 3000 | **2000** |
| timesteps (chunk/stride) | 300 | **200** |
| emitted reference length | 44 nt | **28 nt** |
| references | 16 | **4** |
| min pairwise edit distance | 12 | **8** |

The last row is the one to keep in mind. With half the redundancy, a decode
error here is likelier to land on another *valid* barcode than to fall through
to `unclassified` — the panel is closer to its discrimination limit than the
16-plex is.

### No local amendment

**This bundle is vendored byte-identical to upstream**, and that is deliberate
rather than an oversight. It declares no `boundary.margin` and no
`boundary.clamp_max_shift`, so escpod falls back to margin 200 (`--info` reports
`min adapter_end 2200`) and refuses to clamp.

The *structural* argument from the nbc16 amendment above does carry over: 200
records how `extract_chunks.py` selected the training corpus, not what the
encoder requires, which is a full `chunk` of history and no more. The
*measurement* does not carry over, on three counts:

- **Different chunk.** The affected band here is `[2000, 2200)`, a different
  population of reads from nbc16's `[3000, 3200)`.
- **Four references at min pairwise 8**, against sixteen at 12. Relaxing the
  window is strictly riskier on this panel — see above.
- **Different corpus**, and one whose labels came from a teacher.

Copying `margin: 0` / `clamp_max_shift: 300` across would be asserting a result
nobody measured. Use the run-level `demux.boundary_margin` and
`demux.clamp_max_shift` config keys to evaluate them instead, and amend here
only once the numbers earn it — recording them in this section as the nbc16
section does.

### Published accuracy, and the one caveat that matters most

Exact match to the emitted reference is 0.96075; balanced precision/recall at a
0.97 recovery threshold is 0.9667/0.9393; per-base error 0.0117.

Held out **by run across 17 runs** (`grouped_by: run`), with every one of the 17
groups at median edit distance 0 and median margin 8 — i.e. margin equal to the
panel's minimum pairwise distance, which is what a clean call looks like. This
is a genuine cross-run holdout, and on that count it is *better evidence* than
the nbc16 bundle's confounded pilot.

**But the labels are WarpDemuX's own calls.** From `provenance.json`:

> `label_source`: "WarpDemuX WDX4_tRNA_rna004_v1_0 calls (fpt_boost teacher)
> … **UNGATED** … Additionally restricted to reads that aligned to a tRNA".

So every number above is **student↔teacher agreement, not an independent
accuracy measurement**, and WarpDemuX is a hard ceiling on all of them. Nothing
here says the CRF is more *correct* than WarpDemuX. What it is, is higher
*yield*: upstream's build config records that WarpDemuX declines ~21% of reads
(~6% hard `unclassified`, ~15% below its 0.9 confidence gate) — ~1.9M reads on
one run alone. Yield is the reason to adopt it; correctness has to be argued
from arbiters that do not depend on the teacher (3' adapter purity, tRNA
alignment rate, charging-fraction invariance).

### Four of twelve codes

The panel covers `bc03`, `bc04`, `bc05`, `bc07` only. Upstream:

> "The panel is only 4 codes because that is all the data covers. bc11 and the
> other seven WarpDemuX codes need new sequencing, not a config change."

Those four are exactly the `WDX4_tRNA_rna004_v1_0` kit, so runs using that kit
are fully covered. A run using any other WarpDemuX code has **no path** through
this model. The pipeline validates every sample's barcode against the bundle's
own reference list at DAG-construction time, so this surfaces as an immediate
error naming the four supported codes rather than as a silent misroute.

### Licence review is still open

`provenance.json` carries an explicit hold:

> "**Do not tag a release until this is signed off.** … the twelve barcode
> sequences are WarpDemuX's, transcribed from `ext/WarpDemuX/README.md`, and the
> labels are distilled from WarpDemuX's own calls. Neither the sequences nor the
> teacher are ours." (escapepod-models#40)

`rnabioco/aa-tRNA-seq-pipeline` is a **private** repository, so vendoring the
bundle here is not a public distribution and the hold does not block it.
**If this repository is ever made public, that sign-off is a prerequisite.**

### Upstream provenance is internally inconsistent

Three fields in the upstream `provenance.json` disagree with the rest of it, and
are reproduced here unedited rather than silently corrected:

- `training.dataset` names **3 runs**, while `training.label_source` and
  `metrics.per_group` both describe **17**. The 17-run figure is the one the
  metrics were computed on.
- `runtime` still points at `@v0.1.0`.

Worth an upstream issue; nothing in this pipeline reads those fields.
