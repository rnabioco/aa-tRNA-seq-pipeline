# Rules for modification calling and coverage analysis
# Uses modkit for RNA modification detection and deeptools for coverage
# All outputs are converted to tRNA-only coordinates (adapter positions removed)


rule bam_to_coverage:
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
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
    log:
        os.path.join(outdir, "logs", "bg", "{sample}.txt"),
    threads: 4
    params:
        bg_opts=config["opts"]["coverage"],
        convert_script=os.path.join(SCRIPT_DIR, "convert_to_trna_coords.py"),
        fa=get_validated_reference(),
        offset_5p=get_5p_offset(),
        offset_3p=get_3p_offset(),
    shell:
        """
        bamCoverage \
            -b {input.bam} \
            -o {output.cpm_tmp} \
            --normalizeUsing CPM \
            --outFileFormat bedgraph \
            -bs 1 \
            -p {threads} \
            {params.bg_opts} 2>>{log}

        bamCoverage \
            -b {input.bam} \
            -o {output.counts_tmp} \
            --outFileFormat bedgraph \
            -bs 1 \
            -p {threads} \
            {params.bg_opts} 2>>{log}

        python {params.convert_script} \
            --input {output.counts_tmp} \
            --output {output.counts} \
            --format bedgraph \
            --reference {params.fa} \
            --offset-5p {params.offset_5p} \
            --offset-3p {params.offset_3p} 2>>{log}

        python {params.convert_script} \
            --input {output.cpm_tmp} \
            --output {output.cpm} \
            --format bedgraph \
            --reference {params.fa} \
            --offset-5p {params.offset_5p} \
            --offset-3p {params.offset_3p} 2>>{log}
        """


rule modkit_pileup:
    """ """
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
    output:
        bed=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.pileup.bed.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "pileup", "{sample}"),
    params:
        fa=get_validated_reference(),
        threshold_opts=get_modkit_threshold_opts(),
        convert_script=os.path.join(SCRIPT_DIR, "convert_to_trna_coords.py"),
        offset_5p=get_5p_offset(),
        offset_3p=get_3p_offset(),
    shell:
        """
        modkit pileup \
            --log-filepath {log} \
            --ref {params.fa} \
            {params.threshold_opts} \
            {input.bam} - \
            | python {params.convert_script} \
                --input - \
                --output {output.bed} \
                --format bedmethyl \
                --reference {params.fa} \
                --offset-5p {params.offset_5p} \
                --offset-3p {params.offset_3p}
        """


rule modkit_extract_calls:
    """
    Extract per-read modification calls with optimized thresholds.
    Positions are converted to 1-indexed tRNA-only coordinates.
    """
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "extract_calls", "{sample}"),
    params:
        fa=get_validated_reference(),
        threshold_opts=get_modkit_threshold_opts(),
        convert_script=os.path.join(SCRIPT_DIR, "convert_to_trna_coords.py"),
        offset_5p=get_5p_offset(),
        offset_3p=get_3p_offset(),
    shell:
        """
        modkit extract calls \
            --reference {params.fa} \
            --log-filepath {log} \
            --edge-filter 10 \
            --mapped --pass \
            {params.threshold_opts} \
            {input.bam} - \
            | python {params.convert_script} \
                --input - \
                --output {output.tsv} \
                --format modkit_calls \
                --reference {params.fa} \
                --offset-5p {params.offset_5p} \
                --offset-3p {params.offset_3p}
        """


rule modkit_extract_full:
    """
    Extract full modification information.
    Positions are converted to 1-indexed tRNA-only coordinates.
    """
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_full.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "extract_full", "{sample}"),
    threads: 4
    params:
        fa=get_validated_reference(),
        convert_script=os.path.join(SCRIPT_DIR, "convert_to_trna_coords.py"),
        offset_5p=get_5p_offset(),
        offset_3p=get_3p_offset(),
    shell:
        """
        modkit extract full \
            --threads {threads} \
            --reference {params.fa} \
            --log-filepath {log} \
            --edge-filter 10 \
            --mapped-only \
            {input.bam} - \
            | python {params.convert_script} \
                --input - \
                --output {output.tsv} \
                --format modkit_full \
                --reference {params.fa} \
                --offset-5p {params.offset_5p} \
                --offset-3p {params.offset_3p}
        """
