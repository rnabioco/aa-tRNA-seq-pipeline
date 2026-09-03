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
        maybe_temp(
            os.path.join(outdir, "pod5", "{sample}", "{sample}.pod5"),
            tier="merged_pod5",
        ),
    log:
        os.path.join(outdir, "logs", "merge_pods", "{sample}"),
    threads: 12
    shell:
        """
        escpod merge -t {threads} --force -o {output} {input}
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
            os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam"),
            tier="basecall",
        ),
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}"),
    params:
        model=config["base_calling_model"],
        raw_data_dir=get_basecalling_dir,
        temp_pod5=os.path.join(outdir, "{sample}", "{sample}.pod5"),
        dorado_opts=config["opts"]["dorado"],
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        dorado basecaller --models-directory {params.models_dir} {params.dorado_opts} {params.model} {input.pod5} >{output}
        """


rule ubam_to_fastq:
    """
    extract reads from bam into FASTQ format for alignment
    """
    input:
        rules.rebasecall.output,
    output:
        maybe_temp(
            os.path.join(outdir, "fq", "{sample}", "{sample}.fq.gz"),
            tier="fastq",
        ),
    log:
        os.path.join(outdir, "logs", "ubam_to_fastq", "{sample}"),
    shell:
        """
        samtools fastq {input} | gzip >{output}
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

    For EDX samples, input FASTQ is pre-filtered to matching reads only.
    """
    input:
        reads=get_alignment_fastq,
        idx=rules.bwa_idx.output,
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam"),
            tier="cascade",
        ),
        bai=maybe_temp(
            os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam.bai"),
            tier="cascade",
        ),
    log:
        os.path.join(outdir, "logs", "bwa_align", "{sample}"),
    threads: 16
    params:
        index=get_validated_reference(),
        bwa_opts=config["opts"]["bwa"],
    shell:
        """
        bwa mem -t {threads} {params.bwa_opts} {params.index} {input.reads} \
            | samtools view -F 20 -Sb - \
            | samtools sort -m 2G -@ 4 -o {output.bam}

        samtools index {output.bam}
        """


rule inject_ubam_tags:
    """Transfer all tags from unaligned BAM (dorado) to aligned BAM by read ID.

    Also the point where the sample's identity is written INTO the BAM, because
    it is the first place both demux backends have converged (see the note at
    the top of demux.smk). Two things happen:

    - a constant `BC` tag on every read, so a read stays attributable to its
      barcode after it leaves this directory. Until now the barcode lived only
      in the output path, and the per-read record that could recover it
      (demux/read_ids/) is deleted by the `clean` rule and, on the WarpDemuX
      path, is temp() under the demux_scratch tier.
    - a repaired @RG. bwa builds the aligned header fresh and drops dorado's
      read group, but --all-tags copies the per-read RG:Z straight back, so
      every read pointed at a header line that did not exist. transfer_tags.py
      now splices the source's @RG in and stamps SM/LB/BC onto it.

    Unbarcoded samples get neither BC nor a barcode on the @RG: absence means
    "no demultiplexing", not "unknown barcode".
    """
    input:
        source_bam=rules.rebasecall.output,
        target_bam=rules.bwa_align.output.bam,
        target_bai=rules.bwa_align.output.bai,
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "tagged", "{sample}", "{sample}.tagged.bam"),
            tier="cascade",
        ),
        bai=maybe_temp(
            os.path.join(
                outdir, "bam", "tagged", "{sample}", "{sample}.tagged.bam.bai"
            ),
            tier="cascade",
        ),
    log:
        os.path.join(outdir, "logs", "inject_ubam_tags", "{sample}"),
    threads: 4
    params:
        src=SCRIPT_DIR,
        barcode_arg=lambda wildcards: (
            f"--set-tag BC:Z:{get_sample_barcode_label(wildcards.sample)} "
            f"--rg-barcode {get_sample_barcode_label(wildcards.sample)}"
            if get_sample_barcode_label(wildcards.sample)
            else ""
        ),
        # Records the upstream (escapepod-models) barcode name next to ours
        # whenever the two differ. Since the ldx16 switch neither live panel
        # renames — a sample is configured as the name its bundle emits — so
        # this emits nothing today. Kept because the emitted vocabulary belongs
        # to the bundle: the retired nbc16 panel did differ, and a future one
        # may, and then a BAM tagged `ldx04` should still say what it came
        # from.
        comment_arg=lambda wildcards: (
            f'--comment "aa-tRNA-seq:upstream_barcode='
            f'{get_sample_barcode_upstream(wildcards.sample)}"'
            if get_sample_barcode_upstream(wildcards.sample)
            != get_sample_barcode_label(wildcards.sample)
            else ""
        ),
        rg_library=lambda wildcards: samples[wildcards.sample].get("run_id") or "",
    shell:
        """
        python {params.src}/transfer_tags.py \
            --all-tags \
            --threads {threads} \
            --source {input.source_bam} \
            --target {input.target_bam} \
            --output {output.bam} \
            --rg-sample {wildcards.sample} \
            --rg-library "{params.rg_library}" \
            {params.barcode_arg} \
            {params.comment_arg}

        samtools index -@ {threads} {output.bam}
        """


rule classify_charging:
    """
    Classify charged vs uncharged reads with `escpod classify`.

    Runs on CPU. The model bundle is self-describing — it carries the anchor
    definition, the feature recipe, the k-mer table it is defined against
    (pinned by sha256) and the recommended operating point — so no motif,
    offsets or threshold are passed here. A caller computing the features
    differently gets a wrong answer rather than an error, which is why they are
    not flags. See resources/models/charging/README.md.

    The output BAM is the INPUT records with `cl` (uint8, round(P(charged)*255))
    added, in the same order: dorado's MM/ML modbase tags survive untouched, and
    the file stays coordinate-sorted, so there is no tag round-trip and no
    re-sort. This is what retired the old transfer_bam_tags step, which existed
    only because Remora emitted its score into MM/ML and clobbered the modbase
    calls modkit needs.

    Reads the bundle abstains on (`aligner_arm_depth == 0`) get NO `cl` tag
    rather than a default class, and abstention is charging-correlated — so the
    per-read TSV, which carries a `reason` for every unscored read, is a real
    output and not a debug aid. `read_attrition` folds it in.

    For EDX samples, uses the EDX-filtered POD5 to match the filtered BAM.
    """
    input:
        pod5=get_classification_pod5,
        bam=rules.inject_ubam_tags.output.bam,
        bai=rules.inject_ubam_tags.output.bai,
        reference=get_validated_reference(),
    output:
        charging_bam=maybe_temp(
            os.path.join(
                outdir, "bam", "charging", "{sample}", "{sample}.charging.bam"
            ),
            tier="cascade",
        ),
        charging_bam_bai=maybe_temp(
            os.path.join(
                outdir, "bam", "charging", "{sample}", "{sample}.charging.bam.bai"
            ),
            tier="cascade",
        ),
        calls=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging_calls.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "classify_charging", "{sample}"),
    threads: 8
    params:
        model=get_charging_model(),
        min_mapq=config["charging"]["min_mapq"],
        # LDX samples have no POD5 of their own: input.pod5 is the raw run's
        # files (for dependency tracking) but the tool takes one path — a
        # directory, which it walks recursively.
        pod5_src=get_classification_pod5_arg,
        # escpod writes plain text regardless of extension, so hand it the
        # uncompressed path and gzip afterwards.
        tsv=lambda wildcards, output: output.calls[: -len(".gz")],
    shell:
        """
        escpod classify {params.pod5_src} \
            --bam {input.bam} \
            --reference {input.reference} \
            --model {params.model} \
            --output {output.charging_bam} \
            --tsv {params.tsv} \
            --min-mapq {params.min_mapq} \
            --threads {threads} \
            >{log} 2>&1

        gzip -f {params.tsv}

        samtools index -@ {threads} {output.charging_bam}
        """


rule add_adapter_tags:
    """
    Detect adapter positions in reads using parasail alignment
    and add pt tags (SAM-spec read annotation format) to BAM file.

    pt tag format: start;end;strand;type|start;end;strand;type
    Example: pt:Z:0;24;+;5p_adapter|118;135;+;3p_adapter

    This produces the final BAM with all tags: cl (charging) and pt (adapters).
    """
    input:
        bam=rules.classify_charging.output.charging_bam,
        bai=rules.classify_charging.output.charging_bam_bai,
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "adapter_tagged", "{sample}", "{sample}.bam"),
            tier="cascade",
        ),
        bai=maybe_temp(
            os.path.join(
                outdir, "bam", "adapter_tagged", "{sample}", "{sample}.bam.bai"
            ),
            tier="cascade",
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
        bam=os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam"),
        bai=os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam.bai"),
    log:
        os.path.join(outdir, "logs", "finalize_bam", "{sample}"),
    shell:
        """
        ln -f $(realpath {input.bam}) {output.bam}
        ln -f $(realpath {input.bai}) {output.bai}
        echo "Hardlinked adapter-tagged BAM as final" >{log}
        """
