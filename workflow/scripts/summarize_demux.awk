# Tabulate escapepod's per-read classifications into a per-barcode summary.
#
# Input:  an `escpod demux --classifications` CSV, header row first. Either the
#         single-model shape, read_id,barcode,confidence,..., or the fused
#         multi-model shape, read_id,ldx,ldx_confidence,...,fdx,fdx_confidence,...
# Output: predicted_barcode<TAB>n_reads<TAB>pct
#
# Which column is tallied: `-v col=NAME` names it (an axis of a fused CSV, e.g.
# `fdx`); without it the second column is used, which is `barcode` in the
# single-model shape and the first axis in the fused one.
#
# Column names match parse_warpdemux's output so the QC report reads one format
# regardless of which demux backend produced it.
#
# This lives in a file rather than inline in the rule on purpose. Snakemake
# `shell:` blocks are ordinary Python strings, so backslash escapes in them are
# consumed by Python before the shell ever sees them: an inline "\n" arrives as
# a real newline and breaks the awk string literal it sits in. Keeping the
# program here means awk reads it verbatim.

NR == 1 {
    c = 0
    for (i = 1; i <= NF; i++) if (col != "" && $i == col) c = i
    if (c == 0) {
        if (col != "") {
            print "summarize_demux.awk: no column named " col " in " FILENAME > "/dev/stderr"
            exit 1
        }
        c = 2
    }
    next
}

{ n[$c]++; total++ }

END {
    printf "predicted_barcode\tn_reads\tpct\n"
    for (bc in n) {
        printf "%s\t%d\t%.2f\n", bc, n[bc], 100 * n[bc] / total
    }
}
