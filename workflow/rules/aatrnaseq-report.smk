# Rules for generating QC reports using Quarto
# Requires the 'report' pixi environment: pixi run -e report render-qc


rule render_qc_report:
    """
    Render Quarto QC report for a sample.
    Run with: pixi run -e report snakemake render_qc_report --configfile=config/config-test.yml
    """
    input:
        align_stats=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.align_stats.tsv.gz"
        ),
        charging_prob=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging_prob.tsv.gz"
        ),
        charging_cpm=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.charging.cpm.tsv.gz"
        ),
        bcerror=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.bcerror.tsv.gz"
        ),
    output:
        html=os.path.join(outdir, "reports", "{sample}_qc_report.html"),
    log:
        os.path.join(outdir, "logs", "report", "{sample}.log"),
    params:
        template=os.path.join(SNAKEFILE_DIR, "report", "qc-report.qmd"),
        ml_threshold=config.get("ml-threshold", 200),
        outdir=outdir,
    shell:
        """
        quarto render {params.template} \
            -P sample:{wildcards.sample} \
            -P output_dir:{params.outdir} \
            -P ml_threshold:{params.ml_threshold} \
            --output-dir $(dirname {output.html}) \
            --output $(basename {output.html}) \
            2>&1 | tee {log}
        """


rule render_all_qc_reports:
    """
    Render QC reports for all samples.
    Run with: pixi run -e report snakemake render_all_qc_reports --configfile=config/config-test.yml
    """
    input:
        expand(
            os.path.join(outdir, "reports", "{sample}_qc_report.html"),
            sample=samples.keys(),
        ),
