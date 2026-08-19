#! /usr/bin/env bash
#
# Rebuild the LDX demultiplexing test fixture (.tests/fixtures/ldx-demux).
#
# Unlike make_test_data.sh, the output of this script is COMMITTED to the repo
# rather than uploaded to S3 — it is 4.6 MB, and keeping it versioned alongside
# the configs that reference it stops the two drifting apart. You only need to
# run this to change the fixture's composition.
#
# Requires access to the donor run and its completed pipeline outputs, both of
# which live outside this repository (paths below).
#
# See .tests/fixtures/ldx-demux/README.md for what the fixture contains and why.

set -euo pipefail

# --- Sources --------------------------------------------------------------
# The donor: one SQK-RNA004 flowcell carrying all 16 LDX barcodes, plus all
# seven EDX 3' adapters mixed within each barcode by design. That mixture is
# what makes EDX filtering testable — a barcode holding a single adapter would
# give filter_{fastq,pod5}_by_edx nothing to remove.
PROJ=/beevol/home/jhessel/devel/rnabioco/2026-aars-in-vitro
RUN=$PROJ/data/20260806_EDC_activation_test_EDXs/EDC_act_test/20260806_1113_P2S-00519-A_PBG58575_3669d73d
RESULTS=$PROJ/results/ldx-ligation          # that run, through this pipeline
REF=$PROJ/resources/fa/collapsed.fa         # 47 human tRNAs, raw

OUT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/fixtures/ldx-demux

# Reads per population. Total is ~415 reads at ~11.4 KB/read of POD5, so these
# knobs are what trade fixture size against how much each sample exercises.
KEEP=100          # on-target adapter, per EDX-filtered barcode
OFFTARGET=25      # wrong adapter, per EDX-filtered barcode (must be dropped)
DECOY=40          # reads of a barcode no sample claims
UNCLASSIFIED=25   # reads the CRF could not call

for p in "$RUN" "$RESULTS" "$REF"; do
  [ -e "$p" ] || { echo "missing source: $p" >&2; exit 1; }
done

mkdir -p "$OUT/run/pod5"

# --- 1. Choose the reads --------------------------------------------------
# Selection reads the donor run's OWN pipeline outputs, so every read's
# barcode, 3' adapter and alignment status is known before it is picked. That
# is what lets the fixture be small: no reliance on chance survival.
python3 - "$RESULTS" "$OUT" "$KEEP" "$OFFTARGET" "$DECOY" "$UNCLASSIFIED" <<'PYEOF'
import csv, gzip, subprocess, sys
from collections import defaultdict
from pathlib import Path

results, out = Path(sys.argv[1]), Path(sys.argv[2])
keep, offtarget, decoy, unclassified = (int(a) for a in sys.argv[3:7])

RUN_ID = "20260806_1113_P2S-00519-A_PBG58575_3669d73d"
SAMPLE_OF = {"nbc01": "ldx01_fresh_edx01", "nbc02": "ldx02_fresh_edx02",
             "nbc08": "ldx08_fresh_pool", "nbc05": "ldx05_fresh_edx05"}


def aligned_ids(sample):
    """Read ids that survived alignment and filtering in the donor run."""
    bam = results / "bam" / "final" / sample / f"{sample}.bam"
    p = subprocess.run(["samtools", "view", str(bam)],
                       capture_output=True, text=True, check=True)
    return {ln.split("\t", 1)[0] for ln in p.stdout.splitlines() if ln}


def adapters(sample):
    """read id -> detected 3' adapter, from the donor run's EDX detection."""
    f = results / "demux" / "edx" / sample / f"{sample}.edx_adapters.tsv.gz"
    with gzip.open(f, "rt") as fh:
        return {r["read_id"]: r["adapter_3p"] for r in csv.DictReader(fh, delimiter="\t")}


bc = {}
with open(results / "demux" / "read_ids" / RUN_ID / "classifications.csv") as fh:
    rd = csv.reader(fh); next(rd)
    for rid, barcode, _conf in rd:
        bc[rid] = barcode
print(f"donor run: {len(bc)} classified reads", file=sys.stderr)

by_bc = defaultdict(list)
for rid, b in bc.items():
    by_bc[b].append(rid)

picked, manifest = [], []


def take(ids, n, why):
    # sorted + head, never sampled, so a rebuild selects exactly the same
    # reads. (The POD5 itself will not be byte-identical -- pod5 stamps a
    # fresh file UUID on every write -- but its read set will be.)
    chosen = sorted(ids)[:n]
    if len(chosen) < n:
        sys.exit(f"only {len(chosen)} reads available for {why}, wanted {n}")
    picked.extend(chosen)
    manifest.extend((r, why) for r in chosen)
    print(f"  {why}: {len(chosen)}", file=sys.stderr)


for barcode, want in (("nbc01", "edx01"), ("nbc02", "edx02")):
    sample = SAMPLE_OF[barcode]
    ali, adp = aligned_ids(sample), adapters(sample)
    take([r for r in by_bc[barcode] if adp.get(r) == want and r in ali],
         keep, f"{barcode}:{want}:on_target")
    take([r for r in by_bc[barcode] if adp.get(r) not in (want, None) and r in ali],
         offtarget, f"{barcode}:other_adapter")

ali = aligned_ids(SAMPLE_OF["nbc08"])
take([r for r in by_bc["nbc08"] if r in ali], keep, "nbc08:pool_aligned")
take(by_bc["nbc05"], decoy, "nbc05:unclaimed")
take(by_bc["unclassified"], unclassified, "unclassified")

(out / "read_ids.txt").write_text("\n".join(picked) + "\n")
with open(out / "fixture_manifest.tsv", "w") as fh:
    fh.write("read_id\tselected_as\n")
    for r, why in manifest:
        fh.write(f"{r}\t{why}\n")
print(f"selected {len(picked)} reads (unique {len(set(picked))})", file=sys.stderr)
PYEOF

# --- 2. Cut the POD5 ------------------------------------------------------
# Straight out of the instrument POD5, not out of the demux-split POD5: the
# fixture must be raw input, not another tool's output, or the test would be
# partly asserting that escpod reproduces its own routing.
pod5 filter "$RUN"/pod5/*.pod5 \
  --ids "$OUT/read_ids.txt" \
  --force-overwrite --missing-ok \
  -o "$OUT/run/pod5/ldx_demux_fixture.pod5"
rm -f "$OUT/read_ids.txt"

# --- 3. Reference ---------------------------------------------------------
cp "$REF" "$OUT/collapsed.fa"

echo
echo "fixture rebuilt at $OUT"
du -sh "$OUT"
