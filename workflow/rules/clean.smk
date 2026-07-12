"""
On-demand cleanup of large, regenerable intermediate artifacts.

Run explicitly with:

    snakemake clean --configfile=<your-config.yml>

This deletes the big space-consuming intermediates (merged/split POD5,
rebasecalled uBAMs, FASTQs, intermediate/charging BAMs, demux scratch) from
`output_directory` while KEEPING everything downstream you actually want to
retain: the `summary/` tables (bg/tsv/modkit pileups), the final BAM
(`bam/final/`), any AA-classified BAMs, the built reference, and logs.

Everything removed here is reproducible: re-running the pipeline regenerates
whatever is needed to (re)build a requested output. Because the kept summary
outputs already exist and are up to date, a plain re-run is a no-op — force a
specific target (e.g. `--forcerun merge_pods`) if you actually want the raw
intermediates back.

Note: this is a superset of the auto-`temp()` cleanup controlled by
`cleanup_intermediates` (a bool or list of tier names; see common.smk
`maybe_temp` / `_enabled_cleanup_tiers` and config/README.md). It is meant for
runs completed with intermediates retained (e.g. `cleanup_intermediates: false`,
or a partial tier list), or for reclaiming space on an older output directory.

Unlike the tiered auto-cleanup, this rule also removes `demux/pod5` and (via the
`demux` target) the split/EDX POD5 that the tiers may intentionally keep — so run
it only when you no longer need to re-run classification from the EDX POD5.
"""

# Directories under `output_directory` that are large and cheaply regenerable.
# Paths are relative to outdir. Order does not matter.
CLEAN_TARGETS = [
    "pod5",  # merged per-sample POD5 (merge_pods)
    "fq",  # per-sample FASTQ (ubam_to_fastq)
    "bam/rebasecall",  # rebasecalled uBAM (dorado) — very large
    "bam/aln",  # bwa_align output
    "bam/tagged",  # inject_ubam_tags output
    "bam/charging",  # classify_charging output
    "bam/classified",  # transfer_bam_tags output
    "bam/adapter_tagged",  # add_adapter_tags output (hardlinked into bam/final)
    "demux",  # warpdemux_output, split/edx POD5+FASTQ, read_ids
]

# Directories intentionally preserved (listed for documentation / safety review):
#   summary/          - bg, tsv, modkit pileups (the outputs you keep)
#   bam/final/        - final tagged BAM (key output; hardlink survives adapter_tagged removal)
#   bam/aa_classified/- AA-identity BAMs (key output when aa_identity enabled)
#   reference/        - built + indexed reference
#   logs/             - run logs
#   *.json            - manifest.json, squiggy-session.json


rule clean:
    """Delete large regenerable intermediates from output_directory (keeps summaries + final BAM)."""
    params:
        outdir=outdir,
        targets=" ".join(CLEAN_TARGETS),
    shell:
        r"""
        freed_kb=0
        for rel in {params.targets}; do
            target="{params.outdir}/$rel"
            if [ -L "$target" ]; then
                # e.g. a reuse_outputs_from symlink to another run's data:
                # remove only our link, never the shared source it points to.
                rm -f "$target"
                echo "clean: unlinked symlink $rel (source preserved)"
            elif [ -d "$target" ]; then
                sz_kb=$(du -sk "$target" | cut -f1)
                sz_h=$(du -sh "$target" | cut -f1)
                rm -rf "$target"
                freed_kb=$((freed_kb + sz_kb))
                echo "clean: removed $rel ($sz_h)"
            else
                echo "clean: skip $rel (absent)"
            fi
        done
        awk -v kb="$freed_kb" 'BEGIN {{
            split("KB MB GB TB", u); v=kb; i=1;
            while (v>=1024 && i<4) {{ v/=1024; i++ }}
            printf "clean: reclaimed %.1f %s from {params.outdir}\n", v, u[i]
        }}'
        """
