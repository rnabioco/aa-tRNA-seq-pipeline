"""
Rules for processing raw data from aa-tRNA-seq experiments
"""


rule stage_pod5:
    """
    Stage a sample's raw POD5 files as a directory of symlinks.

    Both consumers of a sample's signal take a directory: dorado basecalls one
    with --recursive, and `escpod classify` walks one recursively, looking each
    aligned read up by id. So nothing needs the signal copied. The rule this
    replaces (`merge_pods`) wrote a merged POD5 per sample -- a full second copy
    of every run, kept for the life of the analysis, which for a 500 GB flowcell
    was 500 GB of duplicate. It existed because the pod5 CLI and Remora wanted a
    single file; neither is used any more, and the LDX path has handed the raw
    run to both tools since it was written.

    The links mirror the source layout, <run>/<pod5_pass|pod5_fail|pod5>/<file>,
    so a sample that pools several runs, or a run that keeps pass and fail reads
    apart, cannot collide on a basename. Targets are canonical (realpath'd)
    paths, so the directory works from anywhere and survives a move of the
    output directory.
    """
    input:
        get_raw_inputs,
    output:
        directory(os.path.join(outdir, "pod5", "{sample}")),
    log:
        os.path.join(outdir, "logs", "stage_pod5", "{sample}"),
    run:
        stage_pod5_links(input, output[0], log[0])


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

    The POD5 input is the sample's signal store: the staged directory of its raw
    run's files (scanned with --recursive), or on the WarpDemuX path the split
    POD5 `escpod filter` wrote for it. LDX samples never reach this rule -- the
    run is basecalled whole by rebasecall_ldx_run and cut per sample afterwards.
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
        dorado_opts=config["opts"]["dorado"],
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
        resume_sh=os.path.join(SCRIPT_DIR, "dorado_basecall_resume.sh"),
        # A staged directory needs the recursive scan; a WarpDemuX split POD5
        # is one file and takes none.
        recursive=lambda wildcards: (
            "" if sample_needs_demux(wildcards.sample) else "--recursive"
        ),
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        # Via the resume wrapper rather than a bare redirect, so a job killed on
        # wall clock costs the tail of a basecall instead of all of it. A
        # per-sample basecall is smaller than the LDX run-level one, but it is
        # the same failure mode and the same fix; see the script.
        bash {params.resume_sh} {output} \
            --models-directory {params.models_dir} {params.dorado_opts} \
            {params.model} {input.pod5} {params.recursive}
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
    Align reads to the tRNA + adapter reference with bwa mem, carrying dorado's
    tags through.

    `samtools fastq -T '*'` writes every uBAM tag into the FASTQ comment, and
    `bwa mem -C` appends that comment to each aligned record. So the move table
    (`mv`, `ns`, `ts`) the charging model reads and the MM/ML modbase calls
    modkit reads arrive on the aligned BAM directly, in one streaming pass with
    no FASTQ on disk. This retired `inject_ubam_tags`, which re-read the uBAM in
    Python to put the same tags back after alignment had dropped them: a 48 GB,
    hours-long job per sample and a fourth full copy of the BAM.

    bwa builds its header from the reference and declares no read groups, while
    `-C` copies dorado's per-read `RG:Z:` through -- so `-H` inserts the uBAM's
    own @RG lines, with SM/LB/BC stamped by stamp_read_groups.py, and every read
    still resolves to a declared read group (the dangling-@RG bug #121 fixed).
    This is also where the sample's identity is written INTO the BAM, because it
    is the first place both demux backends have converged (see the note at the
    top of demux.smk): on a demux run a constant `BC:Z:` goes onto every read
    via the comment, so a read stays attributable to its barcode after it leaves
    this directory, and an @CO records upstream's barcode name whenever it
    differs from ours. Unbarcoded samples get neither.

    Memory: the comment is ~13x the read (4.8 kB against 370 B on the fixture;
    the move table dominates), and bwa holds a whole input batch plus its
    formatted output in memory. `-K` pins the batch at 100 Mbases regardless of
    thread count -- ~600k reads, so a few GB of comment text either side --
    where bwa's default is 10 Mbases PER THREAD. That payload is what PR #86
    measured as an OOM at 48 GB and answered by dropping `-C`; the footprint
    this rule is budgeted for today (160 GB, see cluster/slurm/config.yaml) is
    set by `-k 6` on a dense reference and dwarfs it. `-K` also makes the
    output independent of the thread count.

    For EDX samples only the reads carrying the sample's 3' adapter are aligned:
    `-N` on the uBAM replaces the filtered FASTQ that used to be written for
    them, and the EDX-filtered POD5 that went with it is gone too, since the
    classifier only ever touches reads the BAM names.

    -F 2324 drops unmapped (4), reverse-strand (16), secondary (256) and
    supplementary (2048) records, so the output is primary forward alignments --
    which is what the tagged BAM this replaces contained.
    """
    input:
        ubam=rules.rebasecall.output,
        read_ids=get_alignment_read_ids,
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
        # The @RG/@CO lines handed to bwa -H. Kept beside the BAM: it is tiny,
        # and it is the record of what identity was stamped on this sample.
        header=os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.rg.sam"),
    log:
        os.path.join(outdir, "logs", "bwa_align", "{sample}"),
    threads: 16
    params:
        index=get_validated_reference(),
        bwa_opts=config["opts"]["bwa"],
        src=SCRIPT_DIR,
        batch_bases=100000000,
        rg_args=get_read_group_args,
        bc_tag=lambda wildcards: (
            f"BC:Z:{get_sample_barcode_label(wildcards.sample)}"
            if get_sample_barcode_label(wildcards.sample)
            else ""
        ),
        read_filter=lambda wildcards, input: (
            f"-N {input.read_ids[0]}" if input.read_ids else ""
        ),
    shell:
        """
        python {params.src}/stamp_read_groups.py {input.ubam} {params.rg_args} \
            >{output.header}

        samtools view -u {params.read_filter} {input.ubam} \
            | samtools fastq -T '*' - \
            | awk -v bc="{params.bc_tag}" \
                'NR % 4 == 1 && bc != "" {{ $0 = $0 "\\t" bc }} {{ print }}' \
            | bwa mem -C -K {params.batch_bases} -t {threads} -H {output.header} \
                {params.bwa_opts} {params.index} - \
            | samtools view -u -F 2324 - \
            | samtools sort -m 2G -@ 4 -o {output.bam}

        samtools index {output.bam}
        """


rule calmd:
    """
    Recompute MD/NM against the reference so the aligned BAM carries an `MD`
    tag.

    `bwa mem` does not emit `MD` on its own. The `charging_tcn_sup6_rna004`
    bundle (vendored, not yet wired to any config -- see
    resources/models/charging/README.md) reconstructs its per-read reference
    from `MD` rather than by slicing the reference FASTA by coordinate, and its
    own release notes say a runtime that does the latter must refuse to score
    it: every reference in this panel carries ambiguity codes, and a single
    unresolved one blanks nine consecutive k-mers under the FASTA path.

    `samtools calmd` reads each record against its own reference span, so it
    needs no particular sort order and can run as its own pass right after
    alignment -- `-Q` keeps its per-read debug lines out of the log.

    The per-record MD computation itself is single-threaded (htslib does not
    parallelise it); `--threads` only adds workers for BGZF (de)compression on
    a whole-BAM pass, so this is worth a modest, not a large, thread count.
    """
    input:
        bam=rules.bwa_align.output.bam,
        bai=rules.bwa_align.output.bai,
        fai=get_validated_reference() + ".fai",
    output:
        bam=maybe_temp(
            os.path.join(outdir, "bam", "calmd", "{sample}", "{sample}.calmd.bam"),
            tier="cascade",
        ),
        bai=maybe_temp(
            os.path.join(outdir, "bam", "calmd", "{sample}", "{sample}.calmd.bam.bai"),
            tier="cascade",
        ),
    log:
        os.path.join(outdir, "logs", "calmd", "{sample}"),
    threads: 2
    params:
        index=get_validated_reference(),
    shell:
        """
        samtools calmd -Qb --threads {threads} {input.bam} {params.index} \
            >{output.bam} 2>{log}

        samtools index -@ {threads} {output.bam}
        """


rule classify_charging:
    """
    Classify charged vs uncharged reads with `escpod classify`.

    CPU by default; `charging.gpu: true` scores the windowed (TCN) bundle on
    the GPU instead (escpod >= 0.23.0, see get_charging_device_arg and
    get_charging_escpod_gpu_prefix in common.smk). The GBM/feature-network
    bundles this pipeline ships by default have no GPU path and are
    unaffected either way. The model bundle is self-describing — it carries
    the anchor definition, the feature recipe, the k-mer table it is defined
    against (pinned by sha256) and the recommended operating point — so no
    motif, offsets or threshold are passed here. A caller computing the
    features differently gets a wrong answer rather than an error, which is
    why they are not flags. See resources/models/charging/README.md.

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

    The POD5 argument is the sample's whole signal store -- the staged directory
    of its raw run, a WarpDemuX split POD5, or on LDX the raw run itself -- and
    never a subset of it. classify is driven by the BAM: it looks each aligned
    read's signal up by id, so signal the BAM does not name is never touched.
    That is why an EDX sample needs no EDX-filtered POD5 (it used to get one,
    a copy of the split POD5 kept forever) and why LDX samples have no POD5 of
    their own at all.
    """
    input:
        pod5=get_sample_pod5,
        bam=rules.calmd.output.bam,
        bai=rules.calmd.output.bai,
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
    # Four for CPU, sixteen for GPU -- these no longer scale together.
    #
    # CPU-mode measurement stands as it was: 2026-09-06, identical input
    # (65,821 scored reads, warm 12 GB POD5, `--threads 8`), kept 3.1-4.1
    # cores busy -- wall 55-68 s against 211-270 CPU-seconds. That was never
    # actually the GPU path (this rule is CPU-scored by default; see
    # `resources` below), and it was never re-measured for it either -- the
    # profile comment downstream said so outright ("unmeasured for the GPU
    # path... revisit once profiled").
    #
    # Profiled 2026-09-08 (rnabioco/escapepod-rs#351 and its follow-up): the
    # GPU path was superbatch-serial (CPU prep, then a GPU-only burst, fully
    # alternating -- 13.1% GPU duty cycle on `--threads 4`) before escpod's
    # own fix. With that fixed, the *remaining* bottleneck is almost entirely
    # `escapepod_signal::resquiggle::dp::DpContext::step` (the banded-DP
    # refinement, ~90% of sampled CPU cycles per read-level `perf` profiling)
    # -- a genuinely CPU-bound, embarrassingly-parallel-across-reads cost with
    # no POD5 or GPU dependency, so unlike the CPU path it keeps scaling well
    # past 4 cores. Clean (unprofiled) A/B, same real 55,446-read production
    # sample, same escpod build, `--device gpu`: `--threads 4` 203-214 s,
    # `--threads 16` 87.2 s -- ~2.4x, on top of the ~1.8x the scheduling fix
    # itself already bought (364.7 s fully-serial baseline -> 87.2 s here).
    # Output bit-identical to the fully-serial baseline at every point
    # measured; thread count does not change which reads share a GPU batch.
    # 16 was not swept upward from there -- it is one quarter of a GPU node's
    # 64 cores here, the same fair-share logic `pod5_readers` below already
    # uses for the *other* shared resource this rule contends over.
    threads: lambda wildcards: (16 if config["charging"].get("gpu", False) else 4)
    resources:
        # How many classify jobs may page a POD5 store at once; capped globally
        # in the cluster profiles. This is the throttle that matters. `jobs: 100`
        # with no per-rule cap let 36 of these run together on 2026-09-04
        # (results_v05) and 34 on 09-05 (glnrs_tc_pilot), and per-read cost
        # doubled against runs with 5-6 in flight -- 3603 and 3300 against
        # 1436-1809 us/read, on the SAME escpod version and the same bundle.
        #
        # It is not CPU contention: Slurm allocates the cores. It is 36
        # processes demand-paging one shared multi-hundred-GB POD5 set off
        # BeeGFS. So lower the cap, not the cores, when a run is I/O-starved --
        # and note that cutting `threads` WITHOUT this cap makes it worse, by
        # letting more of them fit at once.
        pod5_readers=1,
        # Conditional on charging.gpu, unlike rebasecall/escapepod_demux's
        # static GPU queue assignment in cluster/{slurm,lsf}/config.yaml: those
        # rules are UNCONDITIONALLY GPU work (dorado basecalls inside
        # escapepod_demux regardless of ldx.gpu), where this rule's GPU-ness is
        # a runtime config toggle the static per-executor profile YAML cannot
        # see. Neither profile declares slurm_partition/gres or
        # lsf_queue/lsf_extra/ngpu for classify_charging, so these win.
        # `cpus_per_task` joined them 2026-09-08, tracking `threads` above for
        # the same reason (GPU-conditional); `cluster/slurm/config.yaml` no
        # longer sets it for this rule. mem_mb/runtime stay governed by the
        # profiles, still unmeasured for the GPU path and kept at the
        # CPU-sized budget as a safe ceiling -- 40 GB and 8 h were sized
        # against CPU-mode's largest corpus (LysRS_U_all20_b4, 3.88M reads),
        # and the GPU path has not been run at anything near that scale yet.
        slurm_partition=lambda wildcards: (
            "gpu" if config["charging"].get("gpu", False) else "rna"
        ),
        slurm_account=lambda wildcards: (
            "gpu_rbi" if config["charging"].get("gpu", False) else "rbi"
        ),
        gres=lambda wildcards: "gpu:1" if config["charging"].get("gpu", False) else "",
        cpus_per_task=lambda wildcards: (
            16 if config["charging"].get("gpu", False) else 4
        ),
        lsf_queue=lambda wildcards: (
            "gpu" if config["charging"].get("gpu", False) else "rna"
        ),
        lsf_extra=lambda wildcards: (
            "-gpu num=1:j_exclusive=yes:mode=exclusive_process"
            if config["charging"].get("gpu", False)
            else ""
        ),
        ngpu=lambda wildcards: 1 if config["charging"].get("gpu", False) else 0,
    params:
        model=get_charging_model(),
        min_mapq=config["charging"]["min_mapq"],
        # Empty unless the config forces a frame, so `auto` stays escpod's own
        # default rather than something this pipeline restates.
        orientation=get_charging_orientation_arg(),
        # The frame to retry with when detection comes back underpowered; ""
        # disables the retry. See get_charging_orientation_fallback.
        orientation_fallback=get_charging_orientation_fallback(),
        fallback_sh=os.path.join(SCRIPT_DIR, "escpod_classify_fallback.sh"),
        # "" under charging.gpu: false; otherwise shadows a GPU-enabled escpod
        # onto PATH for this command only, ahead of the portable CPU build the
        # Snakefile's onstart prefix already put there. See
        # get_charging_escpod_gpu_prefix in common.smk.
        gpu_prefix=get_charging_escpod_gpu_prefix(),
        device=get_charging_device_arg(),
        # LDX samples have no POD5 of their own: input.pod5 is the raw run's
        # files (for dependency tracking) but the tool takes one path — a
        # directory, which it walks recursively.
        pod5_src=get_classification_pod5_arg,
        # escpod writes plain text regardless of extension, so hand it the
        # uncompressed path and gzip afterwards.
        tsv=lambda wildcards, output: output.calls[: -len(".gz")],
    shell:
        """
        # Through the fallback wrapper rather than a bare redirect, so a sample
        # too thin for `auto` to decide the frame is classified with the frame
        # the rest of the run measured instead of failing `rule all` for the
        # whole corpus. It retries ONLY that error; see the script.
        {params.gpu_prefix}bash {params.fallback_sh} "{params.orientation_fallback}" {log} \
            {params.pod5_src} \
            --bam {input.bam} \
            --reference {input.reference} \
            --model {params.model} \
            --output {output.charging_bam} \
            --tsv {params.tsv} \
            --min-mapq {params.min_mapq} \
            {params.orientation} \
            {params.device} \
            --threads {threads}

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

    Hardlinks the adapter-tagged BAM as the final output so that temp() cleanup
    of the upstream BAMs (the `cascade` tier) does not break downstream
    consumers. EDX filtering happens before alignment (bwa_align aligns only the
    sample's own reads), so this is a plain passthrough.
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
