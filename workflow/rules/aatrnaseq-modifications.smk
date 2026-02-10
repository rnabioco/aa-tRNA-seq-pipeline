# Rules for modification calling and coverage analysis
# Uses modkit for RNA modification detection and deeptools for coverage


rule bam_to_coverage:
    input:
        bam=rules.add_adapter_tags.output.bam,
        bai=rules.add_adapter_tags.output.bai,
    output:
        counts_tmp=temp(
            os.path.join(outdir, "summary", "tables", "{sample}", "{sample}.counts.bg")
        ),
        cpm_tmp=temp(
            os.path.join(outdir, "summary", "tables", "{sample}", "{sample}.cpm.bg")
        ),
        counts=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.counts.bg.gz"
        ),
        cpm=os.path.join(outdir, "summary", "tables", "{sample}", "{sample}.cpm.bg.gz"),
    params:
        bg_opts=config["opts"]["coverage"],
    log:
        os.path.join(outdir, "logs", "bg", "{sample}.txt"),
    threads: 4
    shell:
        """
    bamCoverage \
      -b {input.bam} \
      -o {output.cpm_tmp} \
      --normalizeUsing CPM \
      --outFileFormat bedgraph \
      -bs 1 \
      -p {threads} \
      {params.bg_opts}

    bamCoverage \
      -b {input.bam} \
      -o {output.counts_tmp} \
      --outFileFormat bedgraph \
      -bs 1 \
      -p {threads} \
      {params.bg_opts}

    gzip -c {output.counts_tmp} > {output.counts}
    gzip -c {output.cpm_tmp} > {output.cpm}
    """


rule modkit_pileup:
    """
    """
    input:
        bam=rules.add_adapter_tags.output.bam,
        bai=rules.add_adapter_tags.output.bai,
    output:
        bed=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.pileup.bed.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "pileup", "{sample}"),
    params:
        fa=get_validated_reference(),
        threshold_opts=get_modkit_threshold_opts(),
    shell:
        """
    modkit pileup \
        --log-filepath {log} \
        --ref {params.fa} \
        {params.threshold_opts} \
        {input.bam} - \
        | gzip -9 -c > {output.bed}
    """


rule modkit_extract_calls:
    """
    Extract per-read modification calls with optimized thresholds.
    """
    input:
        bam=rules.add_adapter_tags.output.bam,
        bai=rules.add_adapter_tags.output.bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "extract_calls", "{sample}"),
    params:
        fa=get_validated_reference(),
        threshold_opts=get_modkit_threshold_opts(),
    shell:
        """
    modkit extract calls \
        --bgzf \
        --reference {params.fa} \
        --log-filepath {log} \
        --edge-filter 10 \
        --mapped --pass \
        {params.threshold_opts} \
        {input.bam} {output.tsv}
    """


rule modkit_extract_full:
    """
    Extract full modification information.
    """
    input:
        bam=rules.add_adapter_tags.output.bam,
        bai=rules.add_adapter_tags.output.bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_full.tsv.gz"
        ),
    threads: 12
    log:
        os.path.join(outdir, "logs", "modkit", "extract_full", "{sample}"),
    params:
        fa=get_validated_reference(),
    shell:
        """
    modkit extract full \
        --bgzf \
        --threads {threads} \
        --reference {params.fa} \
        --log-filepath {log} \
        --edge-filter 10 \
        --mapped-only \
        {input.bam} {output.tsv}
    """
