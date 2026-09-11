// aa-tRNA-seq pipeline — what changed. Three slides, tool-first.
// Build: typst compile pipeline-changes.typ
// (no typst on PATH here — grab a static binary from
// https://github.com/typst/typst/releases if rebuilding on a node with no
// route out, or `pixi global install typst` where one exists)
//
// Hand-rolled rather than polylux: no package registry needed, so this
// compiles on a compute node with no route to the outside.
//
// Last updated 2026-09-11 for main@078df00 (v0.9.3): symlink staging replaces
// the old merge step, calmd was added ahead of charging, and the headline
// numbers on slide 3 predate the v6.0.0 basecall-model migration — see the
// note on that slide.

#let ink      = rgb("#17191c")
#let muted    = rgb("#7a8189")
#let accent   = rgb("#0f766e")   // teal  — escpod / the new thing
#let danger   = rgb("#b4231f")   // red   — retired
#let good     = rgb("#15803d")   // green — verified
#let hairline = rgb("#dcdfe3")
#let wash     = rgb("#f4f6f6")
#let blush    = rgb("#fdf3f2")

#set page(paper: "presentation-16-9", margin: (x: 54pt, y: 40pt), fill: white)
#set text(font: "Nimbus Sans", size: 14pt, fill: ink, lang: "en")
#set par(leading: 0.58em)
#show raw: set text(font: "Source Code Pro", size: 0.87em)

// ─── components ──────────────────────────────────────────────────────────

#let slide(title, kicker: none, body) = {
  if kicker != none {
    text(size: 11pt, weight: 700, fill: accent, tracking: 1.4pt)[#upper(kicker)]
    v(-6pt)
  }
  text(size: 26pt, weight: 800)[#title]
  line(length: 100%, stroke: 0.8pt + hairline)
  v(8pt)
  body
  pagebreak(weak: true)
}

// 5pt, not 3: descenders (the p in "16-plex") collide with the caption below.
#let stat(value, label, tint: accent, size: 32pt) = stack(
  spacing: 5pt,
  text(size: size, weight: 800, fill: tint)[#value],
  text(size: 10pt, fill: muted)[#label],
)

// One box in a pipeline chain; `tool` prints underneath in mono.
#let step(label, tool: none, tint: ink, bg: white) = block(
  fill: bg,
  stroke: 0.9pt + (if bg == white { hairline } else { tint }),
  radius: 3pt,
  inset: (x: 5pt, y: 6pt),
  width: 100%,
)[
  #set align(center)
  #text(size: 10pt, weight: 600, fill: tint)[#label]
  #v(-4pt)
  #if tool != none {
    text(size: 8pt, font: "Source Code Pro",
         fill: if bg == white { muted } else { tint })[#tool]
  }
]

#let arr = align(horizon + center)[#text(size: 12pt, fill: muted)[→]]

// Lay a chain out evenly: boxes on even indices, arrows on odd.
#let chain(..cells) = {
  let items = cells.pos()
  let cols = ()
  for (i, _) in items.enumerate() {
    cols.push(if calc.odd(i) { 12pt } else { 1fr })
  }
  grid(columns: cols, align: horizon, ..items)
}

#let note(body, tint: accent) = block(
  width: 100%, fill: wash, stroke: (left: 3pt + tint),
  inset: (x: 13pt, y: 9pt), radius: 2pt,
)[#text(size: 12pt)[#body]]

// ─── 1. the streamlining ─────────────────────────────────────────────────

#slide(kicker: "before / after", [One binary replaced three tools])[

  #text(size: 10pt, weight: 700, fill: danger, tracking: 1.2pt)[BEFORE]
  #v(3pt)
  #chain(
    step("merge", tool: "pod5", tint: danger, bg: blush), arr,
    step("demux", tool: "warpdemux"), arr,
    step("split POD5", tool: "pod5 filter", tint: danger, bg: blush), arr,
    step("basecall ×N", tool: "dorado"), arr,
    step("align", tool: "bwa"), arr,
    step("charge", tool: "remora", tint: danger, bg: blush), arr,
    step("retag ML→cl", tool: "transfer_tags", tint: danger, bg: blush), arr,
    step("mods", tool: "modkit"),
  )

  #v(7pt)
  #text(size: 10pt, weight: 700, fill: accent, tracking: 1.2pt)[AFTER]
  #v(3pt)
  #chain(
    step("stage", tool: "symlinks", tint: accent, bg: wash), arr,
    step("demux → sidecar", tool: "escpod demux", tint: accent, bg: wash), arr,
    step("basecall ×1", tool: "dorado"), arr,
    step("split uBAM", tool: "samtools"), arr,
    step("align", tool: "bwa"), arr,
    step("calmd", tool: "samtools", tint: accent, bg: wash), arr,
    step("charge", tool: "escpod classify", tint: accent, bg: wash), arr,
    step("mods", tool: "modkit"),
  )

  #v(8pt)
  #grid(columns: (1.3fr, 1fr), gutter: 24pt,
    text(size: 12.5pt)[
      *#text(fill: accent)[escpod] absorbed the POD5 CLI, Remora, and the leech
      charging backend* — one static Rust binary where there were three Python
      tools and a PyTorch stack. #raw("pixi run setup") no longer installs torch
      at all.
      #v(3pt)
      The #raw("ML")→#raw("cl") retag step went with them: escpod writes its call
      directly onto the aligned records, so dorado's modbase tags are never
      clobbered and never need restoring.
    ],
    grid(columns: (1fr, 1fr), gutter: 12pt, row-gutter: 9pt,
      stat("3 → 1", "tools doing signal work"),
      stat("0", "torch in the default install", tint: good),
      stat("×N → ×1", "basecall passes per run"),
      stat("12 GB", "no longer duplicated", tint: good),
    ),
  )
]

// ─── 2. the models ───────────────────────────────────────────────────────

#slide(kicker: "the two new models", [What they actually do])[
  #grid(columns: (1fr, 1fr), gutter: 26pt,
    [
      #text(size: 10pt, weight: 700, fill: accent, tracking: 1.2pt)[LDX DEMUX]
      #v(1pt)
      #text(size: 16pt, weight: 700)[A CTC-CRF that #text(fill: accent)[basecalls the barcode]]
      #v(5pt)
      #text(size: 11pt)[
        Reads the DNA barcode out of the raw adapter signal and matches the
        decode to references by *edit distance*, rather than classifying
        boundary-gated fingerprints as WarpDemuX does.
        #v(4pt)
        Those fingerprints barely separate — escpod's default segment-mean scores
        *0.03* discriminability, "no barcode signal", and the whole ladder tops
        out near 0.49. A dead end, not a tuning problem.
        #v(4pt)
        The bundle is *self-describing* — its own references, geometry and
        pinned boundary detector — so no #raw("--barcodes") or #raw("--method")
        is ever passed.
      ]
      #v(6pt)
      #grid(columns: (1fr, 1fr), gutter: 10pt,
        stat("16-plex", "against 4–5 for WDX", size: 24pt),
        stat("12", "min edit distance between refs", size: 24pt),
      )
    ],
    [
      #text(size: 10pt, weight: 700, fill: accent, tracking: 1.2pt)[CHARGING]
      #v(1pt)
      #text(size: 16pt, weight: 700)[A recurrent net over the #text(fill: accent)[offset axis]]
      #v(5pt)
      #text(size: 11pt)[
        Per-base features — signal mean plus a z-scored k-mer residual — over
        offsets −8..+24 around the CCA–adapter junction, read as a *sequence*
        rather than a flat vector.
        #v(4pt)
        An MLP control with identical inputs but no offset axis lands below every
        architecture that has one: the gain belongs to *respecting that axis*,
        not to "neural beats trees".
        #v(4pt)
        Labels are *library chemistry, not a model*: chemical ligation cannot
        attach to a deacylated tRNA, enzymatic needs a free 3′-OH.
      ]
      #v(6pt)
      #note(tint: danger)[
        *It abstains* where the common arm did not align — it scores 0.4993
        there. Charging fractions over called reads alone run low.
      ]
    ],
  )
]

// ─── 3. what it bought ───────────────────────────────────────────────────

#slide(kicker: "measured on real data", [What it bought])[
  #grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 16pt,
    stat("0.9906", "charging AUROC, held-out flowcell"),
    stat("99.951%", "demux agreement vs the prior run"),
    stat("395 KB", "sidecar, replacing a full run copy", tint: good),
    stat("+0", "reads unaccounted for in the split", tint: good),
  )

  #v(11pt)
  #text(size: 12pt)[
    Before the arm features were fixed, the charging model was *learning an
    artifact*: with arm availability derived from the aligner, and bwa stopping
    at the adduct, *20.9–25.3%* of charged-library reads lost the entire arm
    against *0.9%* of uncharged — so the model read that missing pattern as a
    class. Counting along the query instead removed it.
  ]

  #v(11pt)
  #text(size: 10pt, weight: 700, fill: danger, tracking: 1.2pt)[
    AND THREE THINGS THAT HAD BEEN FAILING QUIETLY
  ]
  #v(4pt)
  #grid(columns: (1fr, 1fr, 1fr), gutter: 16pt,
    [
      #stat("0", "reads reported lost at charge-calling", tint: danger, size: 28pt)
      #v(1pt)
      #text(size: 11pt)[…while *68 of 624* actually were.]
    ],
    [
      #stat("100%", "of final BAMs were invalid SAM", tint: danger, size: 28pt)
      #v(1pt)
      #text(size: 11pt)[Every read cited an #raw("@RG") that did not exist.]
    ],
    [
      #stat("every", "pipe ran without pipefail", tint: danger, size: 28pt)
      #v(1pt)
      #text(size: 11pt)[A failed align exited 0 with an empty BAM.]
    ],
  )

  #v(10pt)
  #note[
    All figures from 20,400 reads of the 2026-08-06 flowcell, all 16 LDX
    barcodes — not fixtures. #text(fill: muted)[The demux test path has no
    working fixture (issue \#120), which is why. Measured pre-v6 basecall
    model; directional, not current.]
  ]
]
