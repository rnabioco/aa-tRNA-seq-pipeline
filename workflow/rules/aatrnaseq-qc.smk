# Rules for quality control metrics and statistics
# Base calling errors, alignment stats, and signal metrics


rule compute_reference_similarity:
    """
Compute pairwise sequence similarity matrix for reference FASTA.

This QC step identifies potential cross-mapping issues by calculating
all-vs-all sequence similarities using global alignment.
"""
    input:
        fasta=get_raw_reference(),
    output:
        matrix=os.path.join(outdir, "summary", "qc", "reference_similarity.tsv"),
    log:
        os.path.join(outdir, "logs", "qc", "reference_similarity.log"),
    params:
        src=SCRIPT_DIR,
    shell:
        """
        python {params.src}/compute_seq_similarity.py \
            {input.fasta} \
            {output.matrix} \
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


rule align_stats:
    """
extract alignment stats
"""
    input:
        unmapped=rules.rebasecall.output,
        aligned=rules.bwa_align.output.bam,
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
            -i {wildcards.sample} \
            -b {input.unmapped} \
            {input.aligned} \
            {input.classified}
        """


rule remora_signal_stats:
    """
run remora to get signal stats
"""
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
        pod5=get_classification_pod5,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.remora.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "remora", "{sample}"),
    params:
        src=SCRIPT_DIR,
        kmer=config["remora_kmer_table"],
        opts=config["opts"]["remora"],
    shell:
        """
        python {params.src}/extract_signal_metrics.py \
            --pod5_dir {input.pod5} \
            --bam {input.bam} \
            --kmer {params.kmer} \
            --sample_name {wildcards.sample} \
            {params.opts} \
            | gzip -c \
                >{output.tsv}
        """
