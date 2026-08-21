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
        edx_concordance=(
            [os.path.join(outdir, "summary", "edx", "edx_concordance.tsv.gz")]
            if config.get("edx", {}).get("enabled", False)
            else []
        ),
        wdx_summary=(
            expand(
                os.path.join(
                    outdir,
                    "demux",
                    "read_ids",
                    "{run_id}",
                    "demux_summary.tsv.gz",
                ),
                run_id=(
                    list(
                        {
                            info["run_id"]
                            for info in samples.values()
                            if info.get("barcode") and info.get("run_id")
                        }
                    )
                    if is_demux_enabled()
                    else []
                ),
            )
        ),
    output:
        html=os.path.join(outdir, "reports", "qc_report.html"),
    log:
        os.path.join(outdir, "logs", "report", "qc_report.log"),
    params:
        template=os.path.join(SNAKEFILE_DIR, "report", "qc-report.qmd"),
        config_file=workflow.configfiles[0],
        # `charging.ml_threshold`, NOT a top-level `ml-threshold` -- that key has
        # never existed, so this silently fell back to the literal 200 and the QC
        # report disagreed with the CPM tables for anyone who moved the
        # threshold. The two must read the same value: get_cca_trna_cpm takes it
        # from config["charging"]["ml_threshold"] (aatrnaseq-charging.smk).
        ml_threshold=config.get("charging", {}).get("ml_threshold", 200),
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
        cd $(dirname {params.template}) \
            && quarto render $(basename {params.template}) \
                -P config_file:$CONFIG_ABS \
                -P ml_threshold:{params.ml_threshold} \
                --output-dir $OUTPUT_DIR_ABS \
                --output $(basename {output.html}) \
                2>&1 | tee $LOG_ABS

        rm -f "$custom_target"
        """
