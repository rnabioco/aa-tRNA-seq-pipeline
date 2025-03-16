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
        os.path.join(rbc_outdir, "{sample}", "{sample}.pod5"),
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
  requires a GPU
  """
    input:
        rules.merge_pods.output,
    output:
        protected(os.path.join(rbc_outdir, "{sample}", "{sample}.unmapped.bam")),
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}"),
    params:
        model=config["base_calling_model"],
        raw_data_dir=get_basecalling_dir,
        temp_pod5=os.path.join(rbc_outdir, "{sample}", "{sample}.pod5"),
        dorado_opts=config["opts"]["dorado"],
    shell:
        """
    if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
      echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
      export CUDA_VISIBLE_DEVICES
    fi

    dorado basecaller {params.dorado_opts} -v {params.model} {input} > {output}
    """


rule ubam_to_fastq:
    """
  extract reads from bam into FASTQ format for alignment
  """
    input:
        rules.rebasecall.output,
    output:
        os.path.join(outdir, "fastqs", "{sample}.fastq.gz"),
    log:
        os.path.join(outdir, "logs", "ubam_to_fq", "{sample}"),
    shell:
        """
    samtools fastq -T "*" {input} | gzip > {output}
    """


rule bwa_idx:
    input:
        config["fasta"],
    output:
        multiext(config["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        os.path.join(outdir, "logs", "bwa_idx", "log"),
    shell:
        """
    bwa index {input}
    """


rule bwa_align:
    """
  align reads to tRNA references with bwa mem
  """
    input:
        reads=rules.ubam_to_fastq.output,
        idx=rules.bwa_idx.output,
    output:
        bam=os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.unfiltered.bam"),
        bai=os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.unfiltered.bam.bai"),
    params:
        index=config["fasta"],
        bwa_opts=config["opts"]["bwa"],
    log:
        os.path.join(outdir, "logs", "bwa", "{sample}"),
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
        pod5=rules.merge_pods.output,
        bam=rules.bwa_align.output.bam,
    output:
        mod_bam=os.path.join(outdir, "mod_bams", "{sample}_mod.bam"),
        mod_bam_bai=os.path.join(outdir, "mod_bams", "{sample}_mod.bam.bai"),
        txt=os.path.join(outdir, "mod_bams", "{sample}.txt"),
        temp_sorted_bam=temp(os.path.join(outdir, "mod_bams", "{sample}_mod.bam.tmp")),
    log:
        os.path.join(outdir, "logs", "classify_chargin", "{sample}"),
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
      --out-bam {output.mod_bam} \
      --log-filename {output.txt} \
      --reference-anchored \
      --device 0

    # sort the result
    samtools sort {output.mod_bam} > {output.temp_sorted_bam}
    cp {output.temp_sorted_bam} {output.mod_bam}

    samtools index {output.mod_bam}
    """


rule transfer_bam_tags:
    """
  creates final bam with classified reads MM and ML tags and table with charging probability per read
  """
    input:
        source_bam=rules.classify_charging.output.mod_bam,
        target_bam=rules.bwa_align.output.bam,
    output:
        classified_bam=os.path.join(outdir, "classified_bams", "{sample}.bam"),
        classified_bam_bai=os.path.join(outdir, "classified_bams", "{sample}.bam.bai"),
    log:
        os.path.join(outdir, "logs", "transfer_bam_tags", "{sample}"),
    params:
        src=SCRIPT_DIR,
    shell:
        """
    python {params.src}/transfer_tags.py \
      -s {input.source_bam} \
      -t {input.target_bam} \
      -o {output.classified_bam}

    samtools index {output.classified_bam}
    """
