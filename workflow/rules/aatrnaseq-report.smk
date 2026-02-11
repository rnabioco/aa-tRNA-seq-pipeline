# Rules for generating QC reports using Quarto
# Requires the 'report' pixi environment: pixi run -e report render-qc


rule render_combined_qc_report:
    """
    Render combined Quarto QC report with faceted plots for all samples.
    Run with: pixi run -e report snakemake render_combined_qc_report --configfile=config/config-test.yml
    """
    input:
        align_stats=expand(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.align_stats.tsv.gz"
            ),
            sample=samples.keys(),
        ),
        charging_prob=expand(
            os.path.join(
                outdir,
                "summary",
                "tables",
                "{sample}",
                "{sample}.charging_prob.tsv.gz",
            ),
            sample=samples.keys(),
        ),
        charging_cpm=expand(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.charging.cpm.tsv.gz"
            ),
            sample=samples.keys(),
        ),
        bcerror=expand(
            os.path.join(
                outdir, "summary", "tables", "{sample}", "{sample}.bcerror.tsv.gz"
            ),
            sample=samples.keys(),
        ),
        # TODO: odds_ratios temporarily disabled
        # odds_ratios=expand(
        #     os.path.join(
        #         outdir,
        #         "summary",
        #         "tables",
        #         "{sample}",
        #         "{sample}.odds_ratios.tsv.gz",
        #     ),
        #     sample=samples.keys(),
        # ),
    output:
        html=os.path.join(outdir, "reports", "qc_report.html"),
    log:
        os.path.join(outdir, "logs", "report", "qc_report.log"),
    params:
        template=os.path.join(SNAKEFILE_DIR, "report", "qc-report.qmd"),
        config_file=workflow.configfiles[0],
        ml_threshold=config.get("ml-threshold", 200),
        custom_include=config.get("report", {}).get("custom_include", ""),
    shell:
        """
        custom_include="{params.custom_include}"
        custom_target="$(dirname {params.template})/_custom.qmd"
        if [ -n "$custom_include" ] && [ -f "$custom_include" ]; then
            cp "$custom_include" "$custom_target"
        else
            touch "$custom_target"
        fi

        CONFIG_ABS=$(realpath {params.config_file})
        OUTPUT_DIR_ABS=$(realpath $(dirname {output.html}))
        LOG_ABS=$(realpath {log})
        cd $(dirname {params.template}) && \
        quarto render $(basename {params.template}) \
            -P config_file:$CONFIG_ABS \
            -P ml_threshold:{params.ml_threshold} \
            --output-dir $OUTPUT_DIR_ABS \
            --output $(basename {output.html}) \
            2>&1 | tee $LOG_ABS

        rm -f "$custom_target"
        """
