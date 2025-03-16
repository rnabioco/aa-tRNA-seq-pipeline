rule get_cca_trna:
    """
    extract and report charing probability (ML tag) per read
    """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
    output:
        charging_tab=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging_prob.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "get_cca_trna", "{sample}"),
    params:
        src=SCRIPT_DIR,
    shell:
        """
    python {params.src}/get_charging_table.py \
      {input.bam} \
      {output.charging_tab}
    """


rule get_cca_trna_cpm:
    """
    calculate cpm for cca classified trnas
    """
    input:
        charging_tab=rules.get_cca_trna.output.charging_tab,
    output:
        cpm=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging.cpm.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "cca_trna_cpm", "{sample}"),
    params:
        src=SCRIPT_DIR,
        # XXX move `ml_thresh` to config file
        ml_thresh=200,
    shell:
        """
    python {params.src}/get_trna_charging_cpm.py \
      --input {input.charging_tab} \
      --output {output.cpm} \
      --ml-threshold {params.ml_thresh}
    """


rule base_calling_error:
    """
  extract base calling error metrics to tsv file
  """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.bcerror.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "bcerror", "{sample}.bwa"),
    params:
        src=SCRIPT_DIR,
        fa=config["fasta"],
    shell:
        """
    python {params.src}/get_bcerror_freqs.py \
      {input.bam} \
      {params.fa} \
      {output.tsv}
    """


rule align_stats:
    """
  extract alignment stats
  """
    input:
        unmapped=rules.rebasecall.output,
        aligned=rules.bwa_align.output.bam,
        classified=rules.transfer_bam_tags.output.classified_bam,
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


rule bam_to_coverage:
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        counts_tmp=temp(
            os.path.join(outdir, "summary", "tables", "{sample}", "{sample}.counts.bg")
        ),
        cpm_tmp=temp(
            os.path.join(outdir, "summary", "tables", "{sample}", "{sample}.cpm.bg")
        ),
        counts=protected(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.counts.bg.gz"
            )
        ),
        cpm=protected(
            os.path.join(outdir, "summary", "tables", "{sample}", "{sample}.cpm.bg.gz")
        ),
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


rule remora_signal_stats:
    """
  run remora to get signal stats
  """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
        pod5=rules.merge_pods.output,
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
      > {output.tsv}
    """


rule modkit_pileup:
    """
    """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        bed=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.pileup.bed.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "pileup", "{sample}"),
    shell:
        """
    modkit pileup \
        --bgzf \
        --log-filepath {log}
        {input.bam} {output.bed}
    """


rule modkit_extract_calls:
    """
    """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        tsv=os.path.join(
            outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "modkit", "extract_calls", "{sample}"),
    shell:
        """
    modkit extract calls \
        --bgzf \
        --log-filepath {log} \
        --edge-filter 20 \
        --mapped --pass \
        {input.bam} {output.tsv}
    """
