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


rule escpod_align:
    """
    Align reads to the tRNA + adapter reference with `escpod align`, carrying
    dorado's tags through, and write a coordinate-sorted BAM with MD/NM.

    `escpod align` scores every read against every reference (there is no seed
    index, so nothing to build), copies every input tag through byte for byte
    -- the move table (`mv`, `ns`, `ts`) the charging model reads and the MM/ML
    modbase calls modkit reads -- and writes MD/NM matching `samtools calmd`.
    That is what retired bwa_idx, calmd, and the
    `samtools fastq -T '*' | bwa mem -C | samtools sort` pipe (issue #200).
    The `charging_tcn_sup6_rna004` bundle reconstructs its per-read reference
    from MD, so MD on the aligned records is load-bearing, not cosmetic.

    Scoring comes from `opts.escpod_align`, defaulted to bwa's own scoring and
    threshold (`--scoring 1,-1,-2,-1 --min-score 20`). References tied with the
    best score are listed in `XA` (bwa's format) with MAPQ 0 rather than hidden;
    get_charging_table.py splits a tied read's count across that tie set.

    Sample identity: escpod align copies the input's @RG/@CO through and cannot
    add a header line or a tag, so stamp_read_groups.py writes them onto the
    uBAM first, in one pass -- SM/LB/BC on dorado's own @RG (its ID untouched,
    so every per-read RG:Z: still resolves: the dangling-@RG bug #121 fixed) and
    on a demux run a constant `BC:Z:` on every record. This is where a sample's
    identity goes INTO the BAM, because it is the first place both demux
    backends have converged (see the note at the top of demux.smk); unbarcoded
    samples get no BC at all. escpod sniffs its input's format and then reopens
    the path, so it cannot read a pipe: the stamped uBAM is a transient file
    beside the output, deleted as soon as alignment finishes.

    For EDX samples only the reads carrying the sample's 3' adapter are aligned
    (`--read-ids`, like `samtools view -N`); the classifier only ever touches
    reads the BAM names, so no filtered FASTQ or POD5 exists.

    A read scoring below `--min-score` is written unmapped rather than dropped.
    `-F 2324` removes those (and would remove reverse, secondary and
    supplementary records, which the default `--strand forward` without
    `--secondary` never writes), so the output is primary forward alignments --
    what align_stats' `aligned` row, anchor_coverage and read_attrition have
    always counted.
    """
    input:
        ubam=rules.rebasecall.output,
        read_ids=get_alignment_read_ids,
        reference=get_validated_reference(),
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
        os.path.join(outdir, "logs", "escpod_align", "{sample}"),
    threads: lambda wildcards: 16 if is_alignment_gpu() else 8
    resources:
        # GPU-conditional for the same reason as classify_charging: the static
        # per-executor profiles cannot see `alignment.gpu`, so the partition,
        # account, gres and GPU queue are resolved here. mem_mb is not
        # device-dependent and lives in cluster/slurm/config.yaml.
        slurm_partition=lambda wildcards: "gpu" if is_alignment_gpu() else "rna",
        slurm_account=lambda wildcards: "gpu_rbi" if is_alignment_gpu() else "rbi",
        gres=lambda wildcards: "gpu:1" if is_alignment_gpu() else "",
        cpus_per_task=lambda wildcards: 16 if is_alignment_gpu() else 8,
        lsf_queue=lambda wildcards: "gpu" if is_alignment_gpu() else "rna",
        lsf_extra=lambda wildcards: (
            "-gpu num=1:j_exclusive=yes:mode=exclusive_process"
            if is_alignment_gpu()
            else ""
        ),
        ngpu=lambda wildcards: 1 if is_alignment_gpu() else 0,
    params:
        src=SCRIPT_DIR,
        align_opts=config["opts"]["escpod_align"],
        rg_args=get_read_group_args,
        stamped=lambda wildcards, output: output.bam[: -len(".aln.bam")]
        + ".stamped.ubam",
        tmp_dir=lambda wildcards, output: os.path.dirname(output.bam),
        read_filter=lambda wildcards, input: (
            f"--read-ids {input.read_ids[0]}" if input.read_ids else ""
        ),
        gpu_prefix=get_alignment_escpod_gpu_prefix(),
        device=get_alignment_device_arg(),
    shell:
        """
        trap 'rm -f {params.stamped}' EXIT

        python {params.src}/stamp_read_groups.py {input.ubam} {params.rg_args} \
            --compress --output {params.stamped} 2>{log}

        ({params.gpu_prefix}escpod align {params.stamped} \
            --reference {input.reference} \
            --output - \
            --sort coordinate --tmp-dir {params.tmp_dir} \
            {params.align_opts} \
            {params.read_filter} \
            {params.device} \
            --threads {threads}) 2>>{log} \
            | samtools view -b -F 2324 -@ 2 -o {output.bam} -

        samtools index {output.bam}
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
        bam=rules.escpod_align.output.bam,
        bai=rules.escpod_align.output.bai,
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
    #
    # 16 -> 32 2026-09-08 (rnabioco/escapepod-rs#354): that issue set out to
    # fix a suspected CPU/GPU scheduling stall at `--threads 16` and instead
    # found no stall to fix -- a correlated `nvidia-smi dmon` + per-second CPU
    # trace on the same node/binary/sample showed CPU continuously saturated
    # through the whole run, never idling in step with the GPU, so the
    # pipeline's own channel/chunk tuning knobs (already ineffective) were not
    # the lever. Thread count was: same sample, one wide unconfounded
    # allocation so the comparison is apples to apples,
    #   16 threads: 68.1 s, 60% GPU duty
    #   32 threads: 45.1 s, 78% GPU duty  (~1.5x over 16)
    #   48 threads: 42.6 s, 82% GPU duty  (~1.06x over 32 -- near the plateau)
    # 32 is the point past which the curve flattens hard, so it is what moves
    # here rather than 48; see escapepod-rs#354 for the full trace and every
    # number.
    #
    # 32 -> 16 2026-09-11, in two rounds (escapepod-rs#359 then #361, shipped
    # as escpod 0.24.2 then 0.24.3). #359 transposed the dwell-penalty DP's
    # inner loop for autovectorization (bit-identical by construction,
    # 5.1-5.2x on that loop in isolation) -- the same CPU-bound cost #354
    # measured above. First A/B, on a pre-release build of #359
    # (~/scratch/escpod-gpu-duty-cycle/, escpod reporting 0.24.1): isolated
    # single job still favored 32 threads by ~10-35%, but at the NODE level
    # -- a 64-core/4-GPU node fits two `--threads 32` jobs (2 of 4 GPUs used,
    # today's practice) or four `--threads 16` jobs (all 4) -- Scenario B
    # beat Scenario A 2669 vs 1213 reads/s aggregate (~2.2x) and even ran
    # each individual job faster (828 vs 730 reads/s), because using every
    # GPU on the node outweighs any remaining per-job thread advantage. Held
    # at 32 pending release: escapepod-rs#361 was already on `main`,
    # explicitly re-tuning the GPU pipeline's `groups_in_flight`/
    # `prep_chunk` defaults FOR `--threads 16` specifically (its own
    # isolated criterion benchmark: 956 -> 1289 reads/s, ~35%, at that exact
    # thread count) -- landing 16 before that shipped would have meant
    # revising this rule again days later.
    #
    # Second A/B, on the OFFICIAL escpod 0.24.3 release (both #359 and #361;
    # ~/scratch/escpod-0243-sweep/, same real 55,446-read production sample,
    # `--device gpu`): the isolated single-job gap is now GONE -- a clean
    # (cache-warm) rep ties 32 and 16 threads exactly, 53 s / 1046 reads/s
    # both. Node-packing still favors Scenario B on the metric that matters
    # for a real queue of many pending jobs: Scenario A (2x32, 2 GPUs)
    # 1960 reads/s aggregate; Scenario B (4x16, 4 GPUs) 3579 reads/s
    # aggregate (~1.83x) -- despite Scenario B running one synchronized
    # round of exactly 4 jobs ~9.5% slower makespan than Scenario A's exactly
    # 2 (61.97 s vs 56.57 s) and each of its jobs ~6.7% slower in isolation
    # (1232 vs 1320 reads/s). That makespan/per-job cost is an artifact of
    # comparing one synchronized batch of 4 against one of 2; aggregate
    # throughput is what a continuously-refilled queue actually sees, and it
    # is what moves here. `mem_mb` below moved the same day, on the same
    # evidence.
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
        # longer sets it for this rule. `mem_mb` joined 2026-09-11: measured
        # (not carried over from the CPU path unmeasured) at 7.72-8.47 GB
        # MaxRSS across both escpod 0.24.2 (~/scratch/escpod-0242-confirm/)
        # and 0.24.3 (~/scratch/escpod-0243-sweep/), flat across --threads
        # 16/32 and across 1 vs 4 concurrent jobs on the same node
        # (55,446-read real production sample) -- the windowed/TCN GPU path
        # superbatches in bounded chunks rather than holding the whole
        # corpus resident, so this is not expected to scale hard with corpus
        # size the way the CPU path's number does. 16 GB keeps ~2x headroom
        # over the observed peak; `cluster/slurm/config.yaml` no longer sets
        # mem_mb for this rule either. `runtime` stays on the profile's
        # 8h CPU-sized ceiling -- GPU wall time here is under 3 minutes, so
        # it is nowhere near binding. LSF's static `mem_mb=40` (GB, see the
        # N.B. in cluster/lsf/config.yaml) is UNCHANGED: this measurement is
        # Slurm-side only and LSF's unit handling for a rule-level override
        # here hasn't been checked, so it stays the safe unmeasured ceiling.
        slurm_partition=lambda wildcards: (
            "gpu" if config["charging"].get("gpu", False) else "rna"
        ),
        slurm_account=lambda wildcards: (
            "gpu_rbi" if config["charging"].get("gpu", False) else "rbi"
        ),
        gres=lambda wildcards: "gpu:1" if config["charging"].get("gpu", False) else "",
        mem_mb=lambda wildcards: (
            16000 if config["charging"].get("gpu", False) else 40000
        ),
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
    consumers. EDX filtering happens before alignment (escpod_align aligns only
    the sample's own reads), so this is a plain passthrough.
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
