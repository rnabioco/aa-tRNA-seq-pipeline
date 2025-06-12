# workflow/rules/demux.smk
"""
Rules for demultiplexing aa-tRNA-seq samples using warpdemux
"""

rule warpdemux:
    """
    Demultiplex merged pod5 files using warpdemux based on barcode information
    """
    input:
        pod5=os.path.join(outdir, "pod5", "{sample}", "{sample}.pod5")
    output:
        demux_dir=directory(os.path.join(outdir, "pod5", "{sample}", "demuxed")),
        demux_summary=os.path.join(outdir, "pod5", "{sample}", "demuxed", "summary.txt")
    log:
        os.path.join(outdir, "logs", "warpdemux", "{sample}")
    params:
        barcode=lambda wildcards: samples[wildcards.sample].get("barcode", "")
    conda:
        "workflow/envs/warpdemux.yaml" 
    shell:
        """
        if [[ "{params.barcode}" != "" ]]; then
            mkdir -p {output.demux_dir}
            
            # Run warpdemux with the specified barcode
            warpdemux \
                --input {input.pod5} \
                --output-dir {output.demux_dir} \
                --barcode {params.barcode} \
                --summary {output.demux_summary} \
                2> {log}
        else
            # If no barcode specified, create empty directory and summary file
            mkdir -p {output.demux_dir}
            touch {output.demux_summary}
            echo "No barcode specified for sample {wildcards.sample}. Skipping demultiplexing." > {log}
        fi
        """

rule merge_demuxed_pods:
    """
    Merge demultiplexed pod5 files by barcode
    """
    input:
        demux_dir=rules.warpdemux.output.demux_dir
    output:
        merged=os.path.join(outdir, "pod5", "{sample}", "{sample}.{barcode}.demuxed.pod5")
    log:
        os.path.join(outdir, "logs", "merge_demuxed_pods", "{sample}.{barcode}")
    threads: 12
    shell:
        """
        if [[ -d {input.demux_dir}/{wildcards.barcode} ]]; then
            pod5 merge -t {threads} -f -o {output.merged} {input.demux_dir}/{wildcards.barcode}/*.pod5 2> {log}
        else
            echo "No demultiplexed files found for barcode {wildcards.barcode}" > {log}
            touch {output.merged}
        fi
        """

# Modify the rebasecall rule to use demultiplexed pod5 files when available
rule rebasecall_demuxed:
    """
    Rebasecall using demultiplexed pod5 files
    """
    input:
        lambda wildcards: os.path.join(outdir, "pod5", wildcards.sample, f"{wildcards.sample}.{samples[wildcards.sample]['barcode']}.demuxed.pod5") 
        if "barcode" in samples[wildcards.sample] else 
        os.path.join(outdir, "pod5", wildcards.sample, f"{wildcards.sample}.pod5")
    output:
        protected(
            os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam")
        )
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}")
    params:
        model=config["base_calling_model"],
        raw_data_dir=get_basecalling_dir,
        dorado_opts=config["opts"]["dorado"]
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
          echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
          export CUDA_VISIBLE_DEVICES
        fi

        # Use the same rebasecalling command as in the original pipeline
        # but use the demultiplexed pod5 file when available
        resources/tools/dorado/bin/dorado basecaller \
            {params.model} \
            {params.dorado_opts} \
            {input} > {output} 2> {log}
        """