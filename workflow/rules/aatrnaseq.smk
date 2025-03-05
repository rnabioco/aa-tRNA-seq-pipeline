rule merge_pods:
    """
    Merge all fast5/pod5s into a single pod5
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
    Rebasecall using multiple chemistry-specific Dorado models for tRNA-seq or a single model for aa-tRNA-seq.
    Adds the name of the basecalling model used to the output file.
    Requires a GPU.
    """
    input:
        rules.merge_pods.output,
    params:
        models=lambda wildcards: config["base_calling_models"][samples[wildcards.sample]["chemistry"]][samples[wildcards.sample]["sequencing_input"]],
        is_fast5=config["input_format"],
        raw_data_dir=get_basecalling_dir,
        temp_pod5=os.path.join(rbc_outdir, "{sample}", "{sample}.pod5"),
        dorado_opts=config["opts"]["dorado"],
    output:
        lambda wildcards: [
            os.path.join(rbc_outdir, wildcards.sample, f"{wildcards.sample}.{model}.unmapped.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(rbc_outdir, wildcards.sample, f"{wildcards.sample}.unmapped.bam")
        ],
    log:
        lambda wildcards: [
            os.path.join(outdir, "logs", "rebasecall", wildcards.sample, f"{model}.log")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(outdir, "logs", "rebasecall", wildcards.sample, "universal.log")
        ],
    shell:
        """
        for model in {params.models}; do
            output_file="{rbc_outdir}/{wildcards.sample}/{wildcards.sample}.$(basename $model).unmapped.bam"
            log_file="{outdir}/logs/rebasecall/{wildcards.sample}/$(basename $model).log"
            dorado basecaller {params.dorado_opts} -v $model {input} > $output_file 2> $log_file
        done
        """

rule bwa_align:
    """
    Align reads to organism-specific tRNA references with bwa mem
    Processes each BAM file produced by rebasecalling with a different model separately.
    """
    input:
        reads=lambda wildcards: [
            os.path.join(rbc_outdir, wildcards.sample, f"{wildcards.sample}.{model}.unmapped.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(rbc_outdir, wildcards.sample, f"{wildcards.sample}.unmapped.bam")
        ],
        idx=config["references"][samples["{sample}"]["organism"]],
    params:
        output_files=lambda wildcards: [
            os.path.join(outdir, "bams", wildcards.sample, model, f"{wildcards.sample}.{model}.bwa.unfiltered.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(outdir, "bams", wildcards.sample, f"{wildcards.sample}.bwa.unfiltered.bam")
        ],
    output:
        protected("{params.output_files}"),  # ✅ Referencing precomputed params
    log:
        os.path.join(outdir, "logs", "bwa", "{sample}"),
    threads: 12
    shell:
        """
        bwa mem -C -t {threads} {params.bwa_opts} {params.index} {input.reads} \
            | samtools view -F 4 -h \
            | awk '($1 ~ /^@/ || $4 <= 25)' \
            | samtools view -Sb - \
            | samtools sort -o {output}
        
        samtools index {output}
        """

rule transfer_bam_tags:
    """
    Creates final BAM with classified reads MM and ML tags
    Processes each basecalling model separately for tRNA samples (which don't have charging info but may have mod calls)
    """
    input:
        source_bam=lambda wildcards: [
            os.path.join(outdir, "bams", wildcards.sample, model, f"{wildcards.sample}.{model}.bwa.unfiltered.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(outdir, "bams", wildcards.sample, f"{wildcards.sample}.bwa.unfiltered.bam")
        ],
    params:
        output_files=lambda wildcards: [
            os.path.join(outdir, "classified_bams", wildcards.sample, model, f"{wildcards.sample}.{model}.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(outdir, "classified_bams", wildcards.sample, f"{wildcards.sample}.bam")
        ],
    output:
        protected("{params.output_files}"),  # ✅ Referencing precomputed params
    log:
        os.path.join(outdir, "logs", "transfer_bam_tags", "{sample}"),
    shell:
        """
        python {params.src}/transfer_tags.py \
          -s {input.source_bam} \
          -t {input.source_bam} \
          -o {output}

        samtools index {output}
        """

rule align_stats:
    """
    Extract alignment stats
    """
    input:
        unmapped=lambda wildcards: [
            os.path.join(rbc_outdir, wildcards.sample, f"{wildcards.sample}.{model}.unmapped.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(rbc_outdir, wildcards.sample, f"{wildcards.sample}.unmapped.bam")
        ],
        aligned=lambda wildcards: [
            os.path.join(outdir, "bams", wildcards.sample, model, f"{wildcards.sample}.{model}.bwa.unfiltered.bam")
            for model in config["base_calling_models"][samples[wildcards.sample]["chemistry"]]["tRNA"]
        ] if samples[wildcards.sample]["sequencing_input"] == "tRNA" else [
            os.path.join(outdir, "bams", wildcards.sample, f"{wildcards.sample}.bwa.unfiltered.bam")
        ],
    output:
        protected("{params.output_files}"),  # ✅ Referencing precomputed params
    log:
        os.path.join(outdir, "logs", "stats", "{sample}"),
    shell:
        """
        python {params.src}/get_align_stats.py \
          -o {output} \
          -a unmapped aligned classified \
          -i {wildcards.sample} \
          -b {input.unmapped} \
             {input.aligned}
        """