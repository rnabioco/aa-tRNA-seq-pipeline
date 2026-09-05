#! /usr/bin/env bash
#
# Rebuild the dual-index (LDX + FDX) demultiplexing test fixture
# (.tests/fixtures/fdx-demux).
#
# Like make_ldx_test_data.sh, the output is COMMITTED rather than uploaded: it
# is ~4 MB, and versioning it beside the configs that reference it stops the two
# drifting apart. Run this only to change the fixture's composition.
#
# Requires access to the donor run and the escapepod-models analysis of it,
# both outside this repository (paths below). Heavy enough for its own
# allocation -- it scans a 53 GB POD5 directory and an 830 MB BAM:
#
#   srun -p rna -c 8 --mem 16G -t 00:30:00 -J fdx-fixture --comment=fdx-fixture \
#       -- .tests/make_fdx_test_data.sh
#
# See .tests/fixtures/fdx-demux/README.md for what the fixture contains and why.

set -euo pipefail

# --- Sources --------------------------------------------------------------
# The donor: Run2 of the 20260828 FDX pilot (PBK76651), acylated E. coli
# MG1655, five libraries each carrying one 5' FDX code and three 3' LDX codes.
# Run2 is the flowcell the shipped fdx4 model was NOT trained on, which is what
# makes anything measured on this fixture a held-out number.
MODELS=/beevol/home/jhessel/devel/rnabioco/escapepod-models
RUN=$MODELS/data/20260828_FDX1-4_Demux/FDX1-4_Run2/20260828_1246_P2S-01617-A_PBK76651_a6de399f
LIBS=$MODELS/results/fdx/demux/Run2_libraries.parquet   # read -> ldx code -> library -> fdx
CALLS=$MODELS/results/fdx/demux/Run2_calls.csv          # every read's ldx call, `unclassified` included
BAM=/beevol/home/jhessel/scratch/fdx/Run2.bam           # Run2 basecalled and aligned to the raw tRNA bodies
REF=$MODELS/resources/ref/eschColi_K_12_MG1655-mature-tRNAs-collapsed.fa

OUT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/fixtures/fdx-demux
# The pipeline's own pixi env and escpod. A worktree has neither (both are
# gitignored), so point PIPELINE at the main checkout when running from one.
PIPELINE=${PIPELINE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
PY=$PIPELINE/.pixi/envs/default/bin/python
ESCPOD_VERSION=$(awk '/^escpod_version:/ {print $2}' "$PIPELINE/config/config-base.yml")
export PATH=$PIPELINE/resources/tools/escpod/$ESCPOD_VERSION/bin:$PIPELINE/.pixi/envs/default/bin:$PATH

# Reads per population. ~11 KB of POD5 per read, so 365 reads is ~4 MB.
KEEP=100          # aligned reads per claimed (ldx, fdx) pair
DECOY=40          # reads of a pair no sample claims
UNCLASSIFIED=25   # reads the LDX CRF could not call

for p in "$RUN" "$LIBS" "$CALLS" "$BAM" "$REF"; do
  [ -e "$p" ] || { echo "missing source: $p" >&2; exit 1; }
done

mkdir -p "$OUT/run/pod5"

# --- 1. Choose the reads --------------------------------------------------
# Selection reads the donor run's OWN analysis, so every read's ldx code, fdx
# library and alignment status is known before it is picked. The fdx label is
# the LIBRARY's, derived from the ldx call through the pool design -- an
# independent channel from the 5' signal the fdx model reads, which is why it
# can stand as truth for it.
"$PY" - "$LIBS" "$CALLS" "$BAM" "$OUT" "$KEEP" "$DECOY" "$UNCLASSIFIED" <<'PYEOF'
import csv, subprocess, sys
from pathlib import Path

import pandas as pd

libs, calls, bam, out = (Path(a) for a in sys.argv[1:5])
keep, decoy, unclassified = (int(a) for a in sys.argv[5:8])

# Reads that survived alignment (primary, forward) in the donor run. Streamed:
# the BAM is 830 MB and its SAM text several GB, which buffered whole is what
# an earlier version of this step was OOM-killed on at 16 GB.
aligned = set()
with subprocess.Popen(["samtools", "view", str(bam)], stdout=subprocess.PIPE,
                      text=True) as proc:
    for line in proc.stdout:
        aligned.add(line.split("\t", 1)[0])
if proc.returncode:
    sys.exit(f"samtools view {bam} failed with {proc.returncode}")
print(f"donor run: {len(aligned)} aligned reads", file=sys.stderr)

df = pd.read_parquet(libs, columns=["read_id", "code", "fdx"])
print(f"donor run: {len(df)} ldx-classified reads", file=sys.stderr)

picked, manifest = [], []


def take(ids, n, why):
    # sorted + head, never sampled, so a rebuild selects exactly the same
    # reads. (The POD5 will not be byte-identical -- escpod stamps a fresh file
    # UUID on every write -- but its read set will be.)
    chosen = sorted(ids)[:n]
    if len(chosen) < n:
        sys.exit(f"only {len(chosen)} reads available for {why}, wanted {n}")
    picked.extend(chosen)
    manifest.extend((r, why) for r in chosen)
    print(f"  {why}: {len(chosen)}", file=sys.stderr)


for ldx, fdx, why in (("ldx01", "fdx01", "claimed"), ("ldx04", "fdx02", "claimed"),
                      ("ldx07", "fdx03", "claimed")):
    sel = df[(df.code == ldx) & (df.fdx == fdx)].read_id
    take([r for r in sel if r in aligned], keep, f"{ldx}:{fdx}:aligned")

sel = df[(df.code == "ldx10") & (df.fdx == "fdx04")].read_id
take(list(sel), decoy, "ldx10:fdx04:unclaimed")

# The donor run was demultiplexed with barcode_crf_nbc16, so its calls CSV
# carries upstream's `nbc` names; only `unclassified` is read from it.
with open(calls, newline="") as fh:
    uncl = [row["read_id"] for row in csv.DictReader(fh) if row["barcode"] == "unclassified"]
take(uncl, unclassified, "unclassified")

(out / "read_ids.txt").write_text("\n".join(picked) + "\n")
with open(out / "fixture_manifest.tsv", "w") as fh:
    fh.write("read_id\tselected_as\n")
    for r, why in manifest:
        fh.write(f"{r}\t{why}\n")
print(f"selected {len(picked)} reads (unique {len(set(picked))})", file=sys.stderr)
PYEOF

# --- 2. Cut the POD5 ------------------------------------------------------
# Straight out of the instrument POD5: the fixture must be raw input, not
# another tool's output. `escpod filter` takes the run's pod5 directory and
# writes the selected reads to one file, records copied byte-for-byte.
escpod filter "$RUN/pod5" -i "$OUT/read_ids.txt" -o "$OUT/run/pod5/fdx_demux_fixture.pod5"
rm -f "$OUT/read_ids.txt"

# --- 3. Reference ---------------------------------------------------------
cp "$REF" "$OUT/collapsed.fa"

echo
echo "fixture rebuilt at $OUT"
du -sh "$OUT"
