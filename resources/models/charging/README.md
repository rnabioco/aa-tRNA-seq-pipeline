# Vendored charging model

The tRNA aminoacylation (charged vs uncharged) classifier read by
`escpod classify`, committed here rather than fetched — for the same
reasons as [the demux bundles](../demux/README.md): the upstream repository
`rnabioco/escapepod-models` is private, so a fetch needs a `$GITHUB_TOKEN`, and
compute nodes on this cluster have no route to GitHub. Committing the bundle
takes the token, the network and the cache out of the run-time path entirely.

Two bundles are vendored, one per basecalling model. They are **not
interchangeable**, and which one is correct is decided by how the data was
basecalled:

| basecalling model | charging bundle | status |
|---|---|---|
| `rna004_sup@v6.0.0` | `charging_feature_nn_sup6_rna004@v0.1.0` | **default** |
| `rna004_130bps_sup@v5.3.0` | `charging_feature_nn_rna004@v0.1.1` | for data already basecalled with v5.3.0 |

Pointing `charging.model` at the wrong one for your `base_calling_model` is an
error at DAG construction (`charging.basecaller_check`), because it is not a
detail — see below.

## Why the pairing is a rule, not a preference

The charging feature set is `mean + z-scored k-mer residual`, and the EXPECTED
level is predicted from the read's own basecall. So a charging model
substantially detects *how the basecaller fails* at the aminoacyl adduct, and
swapping basecaller changes what its dominant feature means.

Measured upstream on the same reads called two ways: ~0.0097 AUROC and
3.0–3.2 pp of TPR lost, ~0.4 pp of FPR gained (so no threshold recovers it), and
**3.9% of per-read calls flip** — while the aggregate charged fraction moves
0.04 pp. The one number anyone would check when changing basecaller reads "no
change" while one read in 26 answers differently.

Since escapepod-models#106 each bundle **declares** the basecaller it was
trained against, and escpod states it at load:

```
INFO bundle was trained on basecalls from rna004_sup@v6.0.0
     (dorado 2.1.1+d66c17c); scoring reads called with another model is a
     domain shift on the k-mer residual — escpod does not check this
```

escpod carries the declaration without enforcing it. This pipeline enforces it:
`workflow/scripts/basecaller_compat.py` compares model identity at DAG
construction, and `pixi run verify-basecaller` proves byte identity against the
declared `model_sha256`.

**Both bundles require escpod >= 0.19.0.** A charging bundle's schema is
`deny_unknown_fields`, so every older escpod refuses a bundle carrying a
`basecaller` block outright with ``unknown field `basecaller` ``. That is the
schema working as designed, and it is why `escpod_version` moved alongside.

## `charging_feature_nn_sup6_rna004@v0.1.0` (default)

Released 2026-09-01 from
`https://github.com/rnabioco/escapepod-models/releases/tag/charging_feature_nn_sup6_rna004%40v0.1.0`.
Trained on `rna004_sup@v6.0.0` basecalls produced by dorado `2.1.1+d66c17c` —
which is exactly the binary `pixi run setup` installs and exactly the model
`base_calling_model` points at, verified byte-for-byte (190 files,
`9cab42f3…bbadb1`).

Its operating point is **measured on the v6 corpus with its own checkpoint**,
not carried over: `operating_point.cl` is 200, the same value
`charging.ml_threshold` already used, but arrived at independently.

Against the v5.3.0 pairing it trades a little discrimination for substantially
less bias:

| | v5.3.0 pairing | v6 pairing |
|---|---|---|
| test AUROC | 0.9906 | 0.9872 |
| balanced accuracy | 0.9625 | 0.9530 |
| charged reads excluded as unreadable | 4.03% | **1.52%** |

The third row is why this is the default. Reads whose common arm did not align
are abstained on, the exclusion rate is far higher for charged reads, and both
bundles say so in their own `coverage_note` — so a charging fraction over called
reads alone is an UNDERESTIMATE. Cutting the charged-class exclusion rate by
2.6x reduces that bias directly, which matters more here than 0.0034 of AUROC.
Report the no-call rate beside the fraction either way.

## `charging_feature_nn_rna004@v0.1.1`

Released 2026-09-01 from
`https://github.com/rnabioco/escapepod-models/releases/tag/charging_feature_nn_rna004%40v0.1.1`.
This family replaced the Remora `cca_classifier.pt` in pipeline v0.2.0 and was
the default until v6.

v0.1.1 is a **sidecar-only reissue of v0.1.0**: the ONNX and the k-mer table are
byte-identical (`cmp` clean, and the graph's sha256 is unchanged at
`6cfebc2d…820637`), and `metadata.json` differs only by the version string and
the `basecaller` block. Calls do not move — the fixture scores identically, read
for read, under either.

Kept vendored because data basecalled with `rna004_130bps_sup@v5.3.0` is still
correctly scored by it, and re-basecalling an existing run to move to v6 is a
GPU cost rather than a correctness fix. To use it, set BOTH:

```yaml
base_calling_model: "resources/models/rna004_130bps_sup@v5.3.0"
dorado_model: rna004_130bps_sup@v5.3.0
charging:
  model: "resources/models/charging/charging_feature_nn_rna004@v0.1.1"
```

Note its declared dorado is `1.4.0+ba44a013` against the pinned 2.1.1, so this
pairing warns about the major-version difference; the basecalling model, which
is the rule that governs charging calls, matches exactly.

## Both bundles

The bundle is self-describing: `metadata.json` carries the anchor definition,
the feature recipe (offsets, stat layout, standardisation constants), the k-mer
table it is defined against pinned by sha256, the abstain rule, and the
recommended operating point. That is why the run-time invocation passes no
motif, no offsets and no threshold — a caller computing the features
differently gets a wrong answer rather than an error.

```bash
escpod classify reads.pod5 -b aln.bam -r ref.fa \
    -m resources/models/charging/charging_feature_nn_sup6_rna004@v0.1.0 \
    -o out.bam --tsv calls.tsv --min-mapq 0
```

### Contents

| File | sha256 | Notes |
|---|---|---|
| `charging_feature_nn_rna004.onnx` | `6cfebc2d…820637` | LSTM over the offset axis, opset 17. Input `float32 [B,4,33]`, output `[B,2]` logits (uncharged, charged). |
| `9mer_levels_v1.txt` | `1d366c9e…f13e63` | **A symlink**, see below. The k-mer level table the `resid` feature is defined against. |
| `metadata.json` | — | Runtime sidecar: anchor, features, abstain rule, standardisation, operating point. |
| `provenance.json` | — | Training provenance and published metrics. |

Re-check with `sha256sum -c SHA256SUMS.txt` from inside the bundle directory.

### `9mer_levels_v1.txt` is a symlink, not a copy

The released bundle ships its own copy of the ONT RNA004 9-mer level table.
That copy is **byte-identical** to `resources/kmers/9mer_levels_v1.txt`, which
this repository already tracked — the bundle's `kmer_table.source_path` records
that exact path, because it is the file the model was trained against. So the
vendored bundle symlinks it rather than committing a second 7.7 MB copy:

```
9mer_levels_v1.txt -> ../../../kmers/9mer_levels_v1.txt
```

`escpod` opens the table by path and follows the link, and the sha256 pinned in
`metadata.json` still gates it, so a divergence between the two would fail
loudly at load rather than silently changing the feature space. **If
`resources/kmers/9mer_levels_v1.txt` is ever replaced with a different table,
this bundle must get its own copy back** — the residual is defined relative to
this table and a different one makes the model invalid.

### Charging fractions are biased low, and that is declared

The bundle **abstains** on reads where `aligner_arm_depth == 0` (bwa placed no
base of the common arm at all). Those reads get no `cl` tag, rather than a
default class: on that population the model scores balanced accuracy 0.4993 and
calls 100% of the uncharged library charged, so a call there would be noise
shaped like a result.

Arm resolvability is itself charging-correlated — the adduct is what stops the
aligner — so **any charging fraction computed over called reads alone is an
underestimate**, and the no-call rate has to be reported beside it. The pipeline
does this for you: `escpod classify --tsv` emits a row per unscored read
with a `reason` column (`no_aligned_arm`, `no_signal`, `ns_mismatch`), which
lands in `summary/tables/{sample}/{sample}.charging_calls.tsv.gz` and is folded
into `summary/read_attrition.tsv.gz`.

The Remora `cca_classifier` had the same bias. It differed only in not saying
so — its uncallable reads simply never appeared in the output, and the loss was
visible only as a gap between two rows of `align_stats`.

### Operating point

`cl >= 200` (P(charged) >= 0.7824). This is a **recommendation measured on
held-out data, not a property of the model**, and it is what
`charging.ml_threshold` in `config/config-base.yml` is set to. Measured against
library chemistry rather than labels: FPR 1.74% on an enzymatic-only (0%
charged) library, TPR 0.939 on a chemical-only (100% charged) library.

Precision depends on the *sample's* charged fraction, not on the model: 95%
precision needs `cl >= 205` at f=0.25, `cl >= 248` at f=0.10 and `cl >= 254` at
f=0.05. 99% is unreachable below f≈0.25. A caller wanting a different operating
point should move `charging.ml_threshold` and say so.

### Out of distribution

Trained on *S. cerevisiae* tRNA with the edx01/edx02 adapter family. The anchor
needs `...CCA` followed by the common arm `GGCTTCTTCTTGCTCTT`; the 13 nt after
the arm may differ between library adapters and are never read. Reads whose
common arm differs from that sequence are out of distribution.

### The runtime version is pinned alongside the model

This is the per-base-feature ONNX variant. `escpod classify` reads it
from escapepod-rs#231 onward (**v0.10.0+**); an older `escpod` refuses the
bundle with ``missing field `gbm` ``. `escpod_version` in
`config/config-base.yml` is pinned at or above that for this reason.
