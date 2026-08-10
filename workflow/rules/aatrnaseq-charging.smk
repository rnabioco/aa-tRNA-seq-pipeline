# Rules for tRNA charging classification analysis
# Extracts and summarizes charged vs uncharged tRNA classification from Remora ML model


rule get_cca_trna:
    """
    extract and report charing probability (ML tag) per read
    """
    input:
        bam=rules.finalize_bam.output.bam,
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
            --tag cl \
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
