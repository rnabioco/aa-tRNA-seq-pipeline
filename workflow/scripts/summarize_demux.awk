# Tabulate escapepod's per-read classifications into a per-barcode summary.
#
# Input:  read_id,barcode,confidence  (CSV, header row)
# Output: predicted_barcode<TAB>n_reads<TAB>pct
#
# Column names match parse_warpdemux's output so the QC report reads one format
# regardless of which demux backend produced it.
#
# This lives in a file rather than inline in the rule on purpose. Snakemake
# `shell:` blocks are ordinary Python strings, so backslash escapes in them are
# consumed by Python before the shell ever sees them: an inline "\n" arrives as
# a real newline and breaks the awk string literal it sits in. Keeping the
# program here means awk reads it verbatim.

NR > 1 { n[$2]++; total++ }

END {
    printf "predicted_barcode\tn_reads\tpct\n"
    for (bc in n) {
        printf "%s\t%d\t%.2f\n", bc, n[bc], 100 * n[bc] / total
    }
}
