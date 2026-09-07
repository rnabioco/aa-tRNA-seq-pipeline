# Rules for quality control metrics and statistics
# Base calling errors, alignment stats, and read attrition


rule compute_reference_similarity:
    """
    Compute pairwise sequence similarity matrix for reference FASTA.

    This QC step identifies potential cross-mapping issues by calculating
    all-vs-all sequence similarities using global alignment.

    Alignment count is quadratic in the number of reference sequences, so the
    script collapses identical sequences before aligning (lossless) and, when
    `qc.reference_similarity_max_mismatch` is set, additionally collapses
    near-identical ones by Hamming distance. See `config/config-base.yml`.
    """
    input:
        fasta=get_raw_reference(),
    output:
        matrix=os.path.join(outdir, "summary", "qc", "reference_similarity.tsv"),
        clusters=os.path.join(
            outdir, "summary", "qc", "reference_similarity.clusters.tsv"
        ),
    log:
        os.path.join(outdir, "logs", "qc", "reference_similarity.log"),
    # Only ~1.5x from threading even at 16 cores: parasail is a ctypes binding
    # that holds the GIL during marshalling, and the TSV write is serial. The
    # real win is collapsing duplicates, so do not over-request cores here.
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    params:
        src=SCRIPT_DIR,
        max_mismatch=get_similarity_max_mismatch(),
    shell:
        """
        python {params.src}/compute_seq_similarity.py \
            {input.fasta} \
            {output.matrix} \
            --threads {threads} \
            --max-mismatch {params.max_mismatch} \
            --clusters {output.clusters} \
            2>&1 | tee {log}
        """


rule base_calling_error:
    """
    extract base calling error metrics to tsv file
    """
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.bcerror.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "bcerror", "{sample}.bwa"),
    params:
        src=SCRIPT_DIR,
        fa=get_validated_reference(),
        offset_5p=get_5p_offset(),
        offset_3p=get_3p_offset(),
    shell:
        """
        python {params.src}/get_bcerror_freqs.py \
            {input.bam} \
            {params.fa} \
            {output.tsv} \
            --offset-5p {params.offset_5p} \
            --offset-3p {params.offset_3p}
        """


rule bcerror_sites:
    """
    Select candidate modification sites from base-calling error rates.

    Sites are chosen once across all samples so that the per-read calls below
    cover the same positions in every sample and can be compared directly.
    """
    input:
        tsv=expand(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.bcerror.tsv.gz"
            ),
            sample=samples.keys(),
        ),
    output:
        tsv=os.path.join(outdir, "summary", "tables", "bcerror_sites.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "mismatch_calls", "bcerror_sites"),
    params:
        src=SCRIPT_DIR,
        min_error=config.get("mismatch_calls", {}).get("min_error", 0.1),
        min_cov=config.get("mismatch_calls", {}).get("min_coverage", 20),
        min_samples=config.get("mismatch_calls", {}).get("min_samples", 1),
    shell:
        """
        python {params.src}/select_bcerror_sites.py \
            {input.tsv} \
            --output {output.tsv} \
            --min-error {params.min_error} \
            --min-cov {params.min_cov} \
            --min-samples {params.min_samples} \
            2>&1 | tee {log}
        """


rule mismatch_calls:
    """
    Emit per-read base-calling error calls at the selected sites.

    This is the read-level counterpart of the base_calling_error rule, and gives
    downstream co-occurrence analysis (for example modification against charging)
    a per-read signal that does not depend on the modification caller.
    """
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
        sites=rules.bcerror_sites.output.tsv,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.mismatch_calls.tsv.gz"
        ),
        counts=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging_error.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "mismatch_calls", "{sample}"),
    params:
        src=SCRIPT_DIR,
        fa=get_validated_reference(),
        offset_5p=get_5p_offset(),
        offset_3p=get_3p_offset(),
        ml_thresh=config.get("mismatch_calls", {}).get("ml_threshold", 200),
    shell:
        """
        python {params.src}/get_mismatch_calls.py \
            {input.bam} \
            {params.fa} \
            {output.tsv} \
            --counts {output.counts} \
            --sites {input.sites} \
            --offset-5p {params.offset_5p} \
            --offset-3p {params.offset_3p} \
            --ml-threshold {params.ml_thresh} \
            2>&1 | tee {log}
        """


rule align_stats:
    """
    extract alignment stats

    The `classified` row counts reads carrying a `cl` tag, not every read in
    the final BAM: `escpod classify` passes unscored records through
    untouched, so without --require-tag that row would equal `aligned` and the
    charge-calling gate in read_attrition would report zero loss.
    """
    input:
        unmapped=rules.rebasecall.output,
        aligned=rules.calmd.output.bam,
        classified=rules.finalize_bam.output.bam,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.align_stats.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "stats", "{sample}.align_stats"),
    params:
        src=SCRIPT_DIR,
    shell:
        """
        python {params.src}/get_align_stats.py \
            -o {output.tsv} \
            -a unmapped aligned classified \
            --require-tag - - cl \
            -i {wildcards.sample} \
            -b {input.unmapped} \
            {input.aligned} \
            {input.classified}
        """


rule anchor_coverage:
    """
    Does each aligned read span the CCA anchor the charging model reads?

    Runs on the ALIGNED bam, which is temp() under cleanup_intermediates — after a
    run finishes this cannot be recomputed without re-basecalling from POD5, which
    is exactly what recovering it for the 2026-08-06 LDX run required.
    """
    input:
        aligned=rules.calmd.output.bam,
        reference=get_validated_reference(),
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.anchor_coverage.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "stats", "{sample}.anchor_coverage"),
    params:
        src=SCRIPT_DIR,
        adapter_prefix=lambda wildcards: get_adapter_3p()[:25],
    shell:
        """
        python {params.src}/anchor_coverage.py \
            --bam {input.aligned} \
            --reference {input.reference} \
            --sample {wildcards.sample} \
            --adapter-prefix {params.adapter_prefix} \
            --output {output.tsv} 2>&1 | tee {log}
        """


rule read_attrition:
    """
    Where this run's reads were lost, as one table.

    Always produced. Each gate's loss was already derivable, but only by differencing
    rows across files, so nobody did — a 12.37% drop at charge-calling survived every
    run until it was reconstructed by hand (issue #110).

    The charge-calling gate no longer has to be inferred: `escpod classify`
    emits a row per unscored read with a `reason`, so that loss is now broken out
    by cause rather than reported as one unexplained difference.
    """
    input:
        align_stats=expand(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.align_stats.tsv.gz"
            ),
            sample=samples.keys(),
        ),
        anchor=expand(
            os.path.join(
                outdir,
                "summary",
                "tables",
                "{sample}",
                "{sample}.anchor_coverage.tsv.gz",
            ),
            sample=samples.keys(),
        ),
        charging_calls=expand(
            os.path.join(
                outdir,
                "summary",
                "tables",
                "{sample}",
                "{sample}.charging_calls.tsv.gz",
            ),
            sample=samples.keys(),
        ),
        demux=get_demux_summaries,
        assigned=get_assigned_summaries,
    output:
        tsv=os.path.join(outdir, "summary", "read_attrition.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "stats", "read_attrition.log"),
    params:
        src=SCRIPT_DIR,
        demux_arg=lambda wildcards, input: (
            f"--demux-summary {' '.join(input.demux)}" if input.demux else ""
        ),
        # The per-sample join on the escpod path: separates "a code no sample
        # claims, or a dual-index pair the axes did not agree on" from what the
        # basecaller could not read, which the barcode-assigned -> basecalled
        # row otherwise lumps together.
        assigned_arg=lambda wildcards, input: (
            f"--assigned-summary {' '.join(input.assigned)}" if input.assigned else ""
        ),
    shell:
        """
        python {params.src}/read_attrition.py \
            --align-stats {input.align_stats} \
            --anchor-coverage {input.anchor} \
            --charging-calls {input.charging_calls} \
            {params.demux_arg} \
            {params.assigned_arg} \
            --output {output.tsv} 2>&1 | tee {log}
        """
