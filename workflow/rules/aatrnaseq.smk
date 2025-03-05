
rule merge_pods:
    """
  merge all fast5/pod5s into a single pod5
  """
    input:
        get_raw_inputs,
    output:
        os.path.join(rbc_outdir, "{sample}", "{sample}.pod5"),
    log:
        os.path.join(outdir, "logs", "merge_pods", "{sample}"),
    params:
        is_fast5=config["input_format"],
    shell:
        """
    if [ "{params.is_fast5}" == "FAST5" ]; then
      pod5 convert fast5 -f --output {output} {input}
    else
      pod5 merge -f -o {output} {input}
    fi
    """


rule rebasecall:
    """
  rebasecall using using multiple chemistry-specific Dorado models for tRNA-seq or a single model for aa-tRNA-seq
  add the name of the basecalling model used to the output file
  requires a GPU
    """
    input:
        rules.merge_pods.output,
    output:
        lambda wildcards: (
            expand(
                os.path.join(rbc_outdir, "{sample}", "{sample}.{model}.unmapped.bam"),
                model=config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
            ) if samples[wildcards.sample]["sequencing_input"] == "tRNA" else
            os.path.join(rbc_outdir, "{sample}", "{sample}.unmapped.bam")
        ),
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}"),
    params:
        models=lambda wildcards: config["base_calling_models"][samples[wildcards.sample]["chemistry"]][samples[wildcards.sample]["sequencing_input"]],
        is_fast5=config["input_format"],
        raw_data_dir=get_basecalling_dir,
        temp_pod5=os.path.join(rbc_outdir, "{sample}", "{sample}.pod5"),
        dorado_opts=config["opts"]["dorado"],
    shell:
        """
        if [[ "{params.models}" == *,* ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        if [[ "{params.models}" == *,* ]]; then
            for model in {params.models}; do
                output_file="{rbc_outdir}/{wildcards.sample}/{wildcards.sample}.$(basename $model).unmapped.bam"
                dorado basecaller {params.dorado_opts} -v $model {input} > $output_file
            done
        else
            dorado basecaller {params.dorado_opts} -v {params.models} {input} > {output}
        fi
        """



def get_optional_bam_inputs(wildcards):
    sample = wildcards.sample

    if config["input_format"] == "BAM":
        return samples[sample]["raw_files"]
    else:
        return os.path.join(rbc_outdir, sample, sample + ".unmapped.bam")


rule ubam_to_fq:
    """
  extract reads from bam into FASTQ format for alignment
  """
    input:
        get_optional_bam_inputs,
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
  align reads to organism-specific tRNA references with bwa mem
  processes each BAM file produced by rebasecalling with a different model separately
  """
    input:
        reads=lambda wildcards: expand(
            os.path.join(rbc_outdir, "{sample}", "{sample}.{model}.unmapped.bam"),
            sample=wildcards.sample,
            model=config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ) if samples[wildcards.sample]["sequencing_input"] == "tRNA" else 
        os.path.join(rbc_outdir, "{sample}", "{sample}.unmapped.bam"),
        idx=lambda wildcards: config["references"][samples[wildcards.sample]["organism"]],
    output:
        bam=os.path.join(outdir, "bams", "{sample}", "{model}", "{sample}.{model}.bwa.unfiltered.bam"),
        bai=os.path.join(outdir, "bams", "{sample}", "{model}", "{sample}.{model}.bwa.unfiltered.bam.bai"),
    params:
        index=lambda wildcards: config["references"][samples[wildcards.sample]["organism"]],
        bwa_opts=config["opts"]["bwa"],
    log:
        os.path.join(outdir, "logs", "bwa", "{sample}", "{model}"),
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


rule cca_classify:
    """
  run remora trained model to classify charged and uncharged reads
  skip this step if it's just tRNA sequencing not aa-tRNA-seq 
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
        os.path.join(outdir, "logs", "cca_classify", "{sample}"),
    params:
        model=config["remora_cca_classifier"],
        sequencing_input=lambda wildcards: samples[wildcards.sample]["sequencing_input"],  # Fetch sequencing type
    shell:
        """
    if [[ "{params.sequencing_input}" == "aa-tRNA" ]]; then
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

    else
        # Skip classification for tRNA-seq by creating empty output (to avoid MissingOutputException error from Snakemake)
            touch {output.mod_bam} {output.mod_bam_bai} {output.txt}
    """


rule transfer_bam_tags:
    """
    creates final BAM with classified reads MM and ML tags
    processes each basecalling model separately for tRNA samples (which don't have charging info but may have mod calls)
    """
    input:
        source_bam=lambda wildcards: os.path.join(outdir, "bams", "{sample}", "{model}", "{sample}.{model}.bwa.unfiltered.bam"),
        target_bam=lambda wildcards: os.path.join(outdir, "bams", "{sample}", "{model}", "{sample}.{model}.bwa.unfiltered.bam"),
    output:
        classified_bam=os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam"),
        classified_bam_bai=os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam.bai"),
    log:
        os.path.join(outdir, "logs", "transfer_bam_tags", "{sample}", "{model}"),
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



rule get_cca_trna:
    """
    extract and report charing probability (ML tag) per read
    """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
    output:
        charging_tab=os.path.join(
            outdir, "tables", "{sample}", "{sample}.charging_prob.tsv.gz"
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
        cpm=os.path.join(outdir, "tables", "{sample}", "{sample}.charging.cpm.tsv.gz"),
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


rule bcerror:
    """
  extract base calling error metrics to tsv file
  """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        tsv=os.path.join(outdir, "tables", "{sample}", "{sample}.bcerror.tsv.gz"),
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
         unmapped=lambda wildcards: os.path.join(rbc_outdir, "{sample}", "{model}", "{sample}.{model}.unmapped.bam"),
        aligned=lambda wildcards: os.path.join(outdir, "bams", "{sample}", "{model}", "{sample}.{model}.bwa.unfiltered.bam"),
        classified=lambda wildcards: os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam"),
    output:
        tsv=os.path.join(outdir, "tables", "{sample}", "{model}", "{sample}.{model}.align_stats.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "stats", "{sample}", "{model}.align_stats"),
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
    """
    compute coverage, treating each basecalling model separately for tRNA-seq
    """
    input:
        bam=lambda wildcards: os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam"),
        bai=lambda wildcards: os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam.bai"),
    output:
        counts_tmp=temp(os.path.join(outdir, "tables", "{sample}", "{model}", "{sample}.{model}.counts.bg")),
        cpm_tmp=temp(os.path.join(outdir, "tables", "{sample}", "{model}", "{sample}.{model}.cpm.bg")),
        counts=protected(os.path.join(outdir, "tables", "{sample}", "{model}", "{sample}.{model}.counts.bg.gz")),
        cpm=protected(os.path.join(outdir, "tables", "{sample}", "{model}", "{sample}.{model}.cpm.bg.gz")),
    params:
        bg_opts=config["opts"]["coverage"],
    log:
        os.path.join(outdir, "logs", "bg", "{sample}", "{model}.txt"),
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
  run remora to get signal stats (separate treatment of each basecalling model where relevant)
  """
    input:
        bam=lambda wildcards: os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam"),
        bai=lambda wildcards: os.path.join(outdir, "classified_bams", "{sample}", "{model}", "{sample}.{model}.bam.bai"),
        pod5=lambda wildcards: os.path.join(rbc_outdir, "{sample}", "{model}", "{sample}.{model}.pod5"),
    output:
        tsv=os.path.join(outdir, "tables", "{sample}", "{model}", "{sample}.{model}.remora.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "remora", "{sample}", "{model}.log"),
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
