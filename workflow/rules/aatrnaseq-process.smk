"""
Rules for processing raw data from aa-tRNA-seq experiments
"""


rule merge_pods:
    """
  merge pod5s into a single pod5
  """
    input:
        get_raw_inputs,
    output:
        os.path.join(outdir, "pod5", "{sample}", "{sample}.pod5"),
    log:
        os.path.join(outdir, "logs", "merge_pods", "{sample}"),
    threads: 12
    shell:
        """
      pod5 merge -t {threads} -f -o {output} {input}
    """


rule rebasecall:
    """
  rebasecall using different accuracy model

  TODO: remove `-v` to reduce log file size. Removing it cases the call to fail.
  """
    input:
        get_sample_pod5,
    output:
        protected(
            os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam")
        ),
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}"),
    params:
        model=config["base_calling_model"],
        raw_data_dir=get_basecalling_dir,
        temp_pod5=os.path.join(outdir, "{sample}", "{sample}.pod5"),
        dorado_opts=config["opts"]["dorado"],
    shell:
        """
    if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
      echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
      export CUDA_VISIBLE_DEVICES
    fi

    dorado basecaller {params.dorado_opts} {params.model} {input} > {output}
    """


rule ubam_to_fastq:
    """
  extract reads from bam into FASTQ format for alignment
  """
    input:
        rules.rebasecall.output,
    output:
        os.path.join(outdir, "fq", "{sample}", "{sample}.fq.gz"),
    log:
        os.path.join(outdir, "logs", "ubam_to_fastq", "{sample}"),
    shell:
        """
    samtools fastq -T "*" {input} | gzip > {output}
    """


rule bwa_idx:
    """
    Build BWA index for the validated/built reference.
    Depends on reference validation/building completing first.
    """
    input:
        get_validated_reference(),
    output:
        multiext(get_validated_reference(), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        os.path.join(outdir, "logs", "bwa_idx", "log"),
    shell:
        """
    bwa index {input}
    """


rule bwa_align:
    """
    Align reads to tRNA references with bwa mem.
    Uses the validated/built reference.
    """
    input:
        reads=rules.ubam_to_fastq.output,
        idx=rules.bwa_idx.output,
    output:
        bam=os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam"),
        bai=os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam.bai"),
    params:
        index=get_validated_reference(),
        bwa_opts=config["opts"]["bwa"],
    log:
        os.path.join(outdir, "logs", "bwa_align", "{sample}"),
    threads: 12
    shell:
        """
    bwa mem -C -t {threads} {params.bwa_opts} {params.index} {input.reads} \
        | samtools view -F 4 -h \
        | awk '($1 ~ /^@/ || $4 <= 25)' \
        | samtools view -Sb - \
        | samtools sort -o {output.bam}

    samtools index {output.bam}
    """


rule classify_charging:
    """
  run remora trained model to classify charged and uncharged reads
  """
    input:
        pod5=get_sample_pod5,
        bam=rules.bwa_align.output.bam,
    output:
        charging_bam=os.path.join(
            outdir, "bam", "charging", "{sample}", "{sample}.charging.bam"
        ),
        charging_bam_bai=os.path.join(
            outdir, "bam", "charging", "{sample}", "{sample}.charging.bam.bai"
        ),
        temp_sorted_bam=temp(
            os.path.join(
                outdir, "bam", "charging", "{sample}", "{sample}.charging.bam.tmp"
            )
        ),
    log:
        os.path.join(outdir, "logs", "classify_charging", "{sample}"),
    params:
        model=config["remora_cca_classifier"],
    shell:
        """
    if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
      echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
      export CUDA_VISIBLE_DEVICES
    fi

    remora infer from_pod5_and_bam {input.pod5} {input.bam} \
      --model {params.model} \
      --out-bam {output.charging_bam} \
      --log-filename {log} \
      --reference-anchored \
      --device 0

    # sort the result
    samtools sort {output.charging_bam} > {output.temp_sorted_bam}
    cp {output.temp_sorted_bam} {output.charging_bam}

    samtools index {output.charging_bam}
    """


rule transfer_bam_tags:
    """
  creates classified bam with MM and ML tags transferred to CM/CL

  MM/ML tags from the charging classification are transferred to CM/CL so as not to interfere with
  base modifications.
  """
    input:
        source_bam=rules.classify_charging.output.charging_bam,
        target_bam=rules.bwa_align.output.bam,
    output:
        classified_bam=os.path.join(
            outdir, "bam", "classified", "{sample}", "{sample}.bam"
        ),
        classified_bam_bai=os.path.join(
            outdir, "bam", "classified", "{sample}", "{sample}.bam.bai"
        ),
    log:
        os.path.join(outdir, "logs", "transfer_bam_tags", "{sample}"),
    params:
        src=SCRIPT_DIR,
    shell:
        """
    python {params.src}/transfer_tags.py \
      --tags ML MM \
      --rename ML=CL MM=CM \
      --source {input.source_bam} \
      --target {input.target_bam} \
      --output {output.classified_bam}

    samtools index {output.classified_bam}
    """


rule add_adapter_tags:
    """
    Detect adapter positions in reads using parasail alignment
    and add PT tags (SAM-spec read annotation format) to BAM file.

    PT tag format: start;end;strand;type|start;end;strand;type
    Example: PT:Z:0;24;+;5p_adapter|118;135;+;3p_adapter

    This produces the final BAM with all tags: CM/CL (charging) and PT (adapters).
    """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        bam=os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam"),
        bai=os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam.bai"),
    log:
        os.path.join(outdir, "logs", "add_adapter_tags", "{sample}"),
    params:
        src=SCRIPT_DIR,
        adapter_5p=config["adapters"]["five_prime"],
        adapter_3p_args=lambda wc: " ".join(
            f'--adapter-3p "{name}:{seq}"' for name, seq in get_adapter_3p_list()
        ),
        min_score_5p=config["adapters"]["min_score_5p"],
        min_score_3p=config["adapters"]["min_score_3p"],
        infer_5p_flag=(
            "--infer-5p-from-alignment"
            if config["adapters"].get("infer_5p_from_alignment", False)
            else ""
        ),
        max_ref_start_for_5p=config["adapters"].get("max_ref_start_for_5p", 20),
    shell:
        """
    python {params.src}/add_adapter_tags.py \
      -i {input.bam} \
      -o {output.bam} \
      --adapter-5p "{params.adapter_5p}" \
      {params.adapter_3p_args} \
      --min-score-5p {params.min_score_5p} \
      --min-score-3p {params.min_score_3p} \
      {params.infer_5p_flag} \
      --max-ref-start-for-5p {params.max_ref_start_for_5p} \
      2> {log}

    samtools index {output.bam}
    """
