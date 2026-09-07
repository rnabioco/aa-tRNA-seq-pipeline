"""
On-demand cleanup of large, regenerable intermediate artifacts.

Run explicitly with:

    snakemake clean --configfile=<your-config.yml>

This deletes the big space-consuming intermediates (split POD5, rebasecalled
uBAMs, intermediate/charging BAMs, demux scratch -- and, from runs before v0.7,
merged POD5 and FASTQ) from `output_directory` while KEEPING everything
downstream you actually want to retain: the `summary/` tables (bg/tsv/modkit
pileups), the final BAM (`bam/final/`), any AA-classified BAMs, the built
reference, and logs.

Everything removed here is reproducible: re-running the pipeline regenerates
whatever is needed to (re)build a requested output. Because the kept summary
outputs already exist and are up to date, a plain re-run is a no-op — force a
specific target (e.g. `--forcerun rebasecall`) if you actually want the raw
intermediates back.

Note: this is a superset of the auto-`temp()` cleanup controlled by
`cleanup_intermediates` (a bool or list of tier names; see common.smk
`maybe_temp` / `_enabled_cleanup_tiers` and config/README.md). It is meant for
runs completed with intermediates retained (e.g. `cleanup_intermediates: false`,
or a partial tier list), or for reclaiming space on an older output directory.

Unlike the tiered auto-cleanup, this rule also removes `demux/pod5` and (via the
`demux` target) the split/EDX POD5 that the tiers may intentionally keep — so run
it only when you no longer need to re-run classification from the EDX POD5.

Nothing here touches the `.p5s` demux sidecars: those live beside the raw POD5,
outside `output_directory` entirely, and they are the record of which reads
belong to which sample on an LDX run. Deleting one costs a full re-demux, and
this rule is for things that are cheap to rebuild.
"""

# Directories under `output_directory` that are large and cheaply regenerable.
# Paths are relative to outdir. Order does not matter.
CLEAN_TARGETS = [
    "pod5",  # per-sample POD5 symlink dirs (stage_pod5); merged copies on pre-v0.7 runs
    "fq",  # retired ubam_to_fastq output (pre-v0.7 runs; alignment streams from the uBAM now)
    "bam/rebasecall",  # rebasecalled uBAM (dorado) — very large
    "bam/rebasecall_run",  # LDX run-level uBAM (dorado), split into the above
    "bam/aln",  # bwa_align output (carries dorado's tags since v0.7)
    "bam/calmd",  # calmd output (MD/NM added for the TCN charging bundle)
    "bam/tagged",  # retired inject_ubam_tags output (pre-v0.7 runs)
    "bam/charging",  # classify_charging output
    "bam/classified",  # retired transfer_bam_tags output (pre-v0.4.0 runs)
    "bam/adapter_tagged",  # add_adapter_tags output (hardlinked into bam/final)
    "demux",  # warpdemux_output, split/edx POD5+FASTQ, read_ids
]

# Directories intentionally preserved (listed for documentation / safety review):
#   summary/          - bg, tsv, modkit pileups (the outputs you keep)
#   bam/final/        - final tagged BAM (key output; hardlink survives adapter_tagged removal)
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
