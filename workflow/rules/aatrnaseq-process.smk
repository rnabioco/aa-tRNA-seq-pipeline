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
        maybe_temp(os.path.join(outdir, "pod5", "{sample}", "{sample}.pod5")),
    log:
        os.path.join(outdir, "logs", "merge_pods", "{sample}"),
    threads: 12
    shell:
        """
        rm -f {output}
        escpod merge -t {threads} -o {output} {input} 2>{log}
        """


rule download_mod_models:
    """
Download dorado modified bases models if not already present.
Runs once on the submission node before basecalling to avoid
race conditions from parallel GPU jobs downloading simultaneously.
"""
    output:
        sentinel=os.path.join(PIPELINE_DIR, "resources", "models", ".mod_models_ready"),
    params:
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
        base_model=config["dorado_model"],
        mod_bases=get_modified_bases(),
    shell:
        """
        for mod in {params.mod_bases}; do
            model="{params.base_model}_${{mod}}@v1"
            model_path="{params.models_dir}/$model"
            if [ -f "$model_path/.downloaded" ]; then
                echo "Model $model already downloaded"
            else
                echo "Downloading $model..."
                dorado download --model "$model" --models-directory {params.models_dir}
                touch "$model_path/.downloaded"
                echo "Downloaded $model"
            fi
        done
        touch {output.sentinel}
        """


rule rebasecall:
    """
rebasecall using different accuracy model

TODO: remove `-v` to reduce log file size. Removing it cases the call to fail.
"""
    input:
        pod5=get_sample_pod5,
        mod_models=rules.download_mod_models.output.sentinel,
    output:
        maybe_temp(
            os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam")
        ),
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}"),
    params:
        model=config["base_calling_model"],
        raw_data_dir=get_basecalling_dir,
        temp_pod5=os.path.join(outdir, "{sample}", "{sample}.pod5"),
        # `dorado_opts_override` lets a run swap the dorado options wholesale
        # (e.g. drop --modified-bases for a canonical basecall when comparing
        # dorado versions whose mod models differ). Falls back to opts.dorado.
        dorado_opts=config.get("dorado_opts_override", config["opts"]["dorado"]),
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        # stdout is the BAM, so send dorado's stderr (errors + progress) to the
        # rule log; otherwise a failed basecall leaves an empty log and an
        # opaque "Reason: Unknown" on the cluster.
        dorado basecaller --models-directory {params.models_dir} {params.dorado_opts} {params.model} {input.pod5} >{output} 2>{log}
        """


rule ubam_to_fastq:
    """
extract reads from bam into FASTQ format for alignment
"""
    input:
        rules.rebasecall.output,
    output:
        maybe_temp(os.path.join(outdir, "fq", "{sample}", "{sample}.fq.gz")),
    log:
        os.path.join(outdir, "logs", "ubam_to_fastq", "{sample}"),
    shell:
        """
        samtools fastq {input} 2>{log} | gzip >{output}
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
        bwa index {input} 2>{log}
        """


rule bwa_align:
    """
Align reads to tRNA references with bwa mem.
Uses the validated/built reference.

For EDX samples, input FASTQ is pre-filtered to matching reads only.
"""
    input:
        reads=get_alignment_fastq,
        idx=rules.bwa_idx.output,
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam")
        ),
        bai=maybe_temp(
            os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam.bai")
        ),
    log:
        os.path.join(outdir, "logs", "bwa_align", "{sample}"),
    threads: 16
    params:
        index=get_validated_reference(),
        bwa_opts=config["opts"]["bwa"],
    shell:
        """
        {{
            bwa mem -t {threads} {params.bwa_opts} {params.index} {input.reads} \
                | samtools view -F 20 -Sb - \
                | samtools sort -m 2G -@ 4 -o {output.bam}

            samtools index {output.bam}
        }} 2>{log}
        """


rule inject_ubam_tags:
    """Transfer all tags from unaligned BAM (dorado) to aligned BAM by read ID."""
    input:
        source_bam=rules.rebasecall.output,
        target_bam=rules.bwa_align.output.bam,
        target_bai=rules.bwa_align.output.bai,
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "tagged", "{sample}", "{sample}.tagged.bam")
        ),
        bai=maybe_temp(
            os.path.join(
                outdir, "bam", "tagged", "{sample}", "{sample}.tagged.bam.bai"
            )
        ),
    log:
        os.path.join(outdir, "logs", "inject_ubam_tags", "{sample}"),
    threads: 4
    params:
        src=SCRIPT_DIR,
    shell:
        """
        {{
            python {params.src}/transfer_tags.py \
                --all-tags \
                --threads {threads} \
                --source {input.source_bam} \
                --target {input.target_bam} \
                --output {output.bam}

            samtools index -@ {threads} {output.bam}
        }} 2>{log}
        """


rule classify_charging:
    """
run remora trained model to classify charged and uncharged reads
runs on CPU by default (no --device flag)

For EDX samples, uses the EDX-filtered POD5 to match the filtered BAM.
"""
    input:
        pod5=get_classification_pod5,
        bam=rules.inject_ubam_tags.output.bam,
    output:
        charging_bam=maybe_temp(
            os.path.join(
                outdir, "bam", "charging", "{sample}", "{sample}.charging.bam"
            )
        ),
        charging_bam_bai=maybe_temp(
            os.path.join(
                outdir, "bam", "charging", "{sample}", "{sample}.charging.bam.bai"
            )
        ),
        temp_sorted_bam=temp(
            os.path.join(
                outdir, "bam", "charging", "{sample}", "{sample}.charging.bam.tmp"
            )
        ),
    log:
        os.path.join(outdir, "logs", "classify_charging", "{sample}"),
    threads: 8
    params:
        model=config["remora_cca_classifier"],
    shell:
        """
        remora infer from_pod5_and_bam {input.pod5} {input.bam} \
            --model {params.model} \
            --out-bam {output.charging_bam} \
            --log-filename {log} \
            --reference-anchored \
            --num-extract-alignment-workers 2 \
            --num-prepare-read-workers 2 \
            --num-prepare-nn-input-workers 2 \
            --num-post-process-workers 2

        # sort the result
        samtools sort -@ {threads} {output.charging_bam} >{output.temp_sorted_bam}
        cp {output.temp_sorted_bam} {output.charging_bam}

        samtools index {output.charging_bam}
        """


rule classify_charging_leech:
    """
run leech trained model to classify charged and uncharged reads
GPU-accelerated alternative to remora (requires leech installed from resources/leech)

For EDX samples, uses the EDX-filtered POD5 to match the filtered BAM.
"""
    input:
        pod5=get_classification_pod5,
        bam=rules.inject_ubam_tags.output.bam,
        # leech reads the index (bam.mapped) for reference-anchored inference;
        # the .bai is a temp output of inject_ubam_tags, so require it here or
        # snakemake deletes it before this rule runs.
        bai=rules.inject_ubam_tags.output.bai,
        # reference-anchored mode needs the actual reference sequences; the BAM
        # @SQ header carries only names/lengths.
        reference=get_validated_reference(),
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
        os.path.join(outdir, "logs", "classify_charging_leech", "{sample}"),
    threads: 4
    params:
        model=config["remora_cca_classifier"],
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        # motif/motif-offset are intentionally omitted: leech auto-reads them
        # from the model config (cca_classifier.pt was trained with
        # motif-offset=3, and leech refuses a mismatched override).
        leech predict \
            --model {params.model} \
            --pod5 {input.pod5} \
            --bam {input.bam} \
            --output {output.charging_bam} \
            --device cuda \
            --reference-anchored \
            --reference-fasta {input.reference} \
            --workers 4 \
            --batch-size 512 \
            2>&1 | tee {log}

        # sort the result
        samtools sort -@ {threads} {output.charging_bam} >{output.temp_sorted_bam}
        cp {output.temp_sorted_bam} {output.charging_bam}

        samtools index {output.charging_bam}
        """


rule classify_aa_identity:
    """
Run leech one-vs-all bundle to predict amino acid identity per read.
Adds aa (predicted AA), ac (confidence), pn (pair names), pp (pair probs) tags.
"""
    input:
        pod5=get_classification_pod5,
        bam=rules.inject_ubam_tags.output.bam,
    output:
        bam=os.path.join(outdir, "bam", "aa_classified", "{sample}", "{sample}.bam"),
        bai=os.path.join(outdir, "bam", "aa_classified", "{sample}", "{sample}.bam.bai"),
        temp_sorted=temp(
            os.path.join(
                outdir, "bam", "aa_classified", "{sample}", "{sample}.bam.tmp"
            )
        ),
    log:
        os.path.join(outdir, "logs", "classify_aa_identity", "{sample}"),
    threads: 4
    params:
        bundle=config.get("aa_identity", {}).get("bundle", ""),
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        leech predict \
            --bundle {params.bundle} \
            --all \
            --pod5 {input.pod5} \
            --bam {input.bam} \
            --output {output.bam} \
            --device cuda \
            --workers 4 \
            --batch-size 512 \
            --raw \
            2>&1 | tee {log}

        samtools sort -@ {threads} {output.bam} >{output.temp_sorted}
        cp {output.temp_sorted} {output.bam}
        samtools index {output.bam}
        """


rule transfer_bam_tags:
    """
creates classified bam with MM and ML tags transferred to cm/cl

MM/ML tags from the charging classification are transferred to cm/cl so as not to interfere with
base modifications.
"""
    input:
        source_bam=rules.classify_charging.output.charging_bam,
        target_bam=rules.inject_ubam_tags.output.bam,
    output:
        classified_bam=maybe_temp(
            os.path.join(outdir, "bam", "classified", "{sample}", "{sample}.bam")
        ),
        classified_bam_bai=maybe_temp(
            os.path.join(outdir, "bam", "classified", "{sample}", "{sample}.bam.bai")
        ),
    log:
        os.path.join(outdir, "logs", "transfer_bam_tags", "{sample}"),
    threads: 4
    params:
        src=SCRIPT_DIR,
    shell:
        """
        python {params.src}/transfer_tags.py \
            --tags ML MM \
            --rename ML=cl MM=cm \
            --threads {threads} \
            --source {input.source_bam} \
            --target {input.target_bam} \
            --output {output.classified_bam}

        samtools index -@ {threads} {output.classified_bam}
        """


rule add_adapter_tags:
    """
Detect adapter positions in reads using parasail alignment
and add pt tags (SAM-spec read annotation format) to BAM file.

pt tag format: start;end;strand;type|start;end;strand;type
Example: pt:Z:0;24;+;5p_adapter|118;135;+;3p_adapter

This produces the final BAM with all tags: cm/cl (charging) and pt (adapters).
"""
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "adapter_tagged", "{sample}", "{sample}.bam")
        ),
        bai=maybe_temp(
            os.path.join(
                outdir, "bam", "adapter_tagged", "{sample}", "{sample}.bam.bai"
            )
        ),
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
            2>{log}

        samtools index {output.bam}
        """


rule finalize_bam:
    """
Produce the final BAM for downstream analysis.

EDX filtering now happens early in the pipeline (before alignment) via
the detect_edx_adapters / filter_fastq_by_edx / filter_pod5_by_edx rules.
This rule hardlinks the adapter-tagged BAM as the final output so that
temp() cleanup of upstream BAMs doesn't break downstream consumers.
"""
    input:
        bam=rules.add_adapter_tags.output.bam,
        bai=rules.add_adapter_tags.output.bai,
    output:
        bam=maybe_temp(os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam")),
        bai=maybe_temp(
            os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam.bai")
        ),
    log:
        os.path.join(outdir, "logs", "finalize_bam", "{sample}"),
    shell:
        """
        ln -f $(realpath {input.bam}) {output.bam}
        ln -f $(realpath {input.bai}) {output.bai}
        echo "Hardlinked adapter-tagged BAM as final" >{log}
        """
