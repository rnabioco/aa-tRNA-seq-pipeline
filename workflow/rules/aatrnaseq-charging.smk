# Rules for tRNA charging classification analysis
# Extracts and summarizes the charged vs uncharged calls (`cl` tag) that
# `escpod classify` wrote onto the BAM.


rule get_cca_trna:
    """
    extract and report charging probability (`cl` tag) per read

    Reads that the model abstained on carry no `cl` tag and so do not appear
    here. Their count and cause are in {sample}.charging_calls.tsv.gz and
    read_attrition.tsv.gz — read the two together, because abstention is
    charging-correlated and this table alone understates the charged fraction.
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
        ml_thresh=config["charging"]["ml_threshold"],
    shell:
        """
        python {params.src}/get_trna_charging_cpm.py \
            --input {input.charging_tab} \
            --output {output.cpm} \
            --ml-threshold {params.ml_thresh}
        """
