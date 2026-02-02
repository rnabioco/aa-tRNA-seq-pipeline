# Rules for pairwise modification odds ratio analysis
# Computes odds ratios between all pairs of positions (+ charging) across tRNA references


rule compute_odds_ratios:
    """
    Compute pairwise modification odds ratios across tRNA references.

    For each pair of positions (plus charging status), builds a 2x2
    contingency table across all tRNA genes in the sample and computes
    odds ratios with Fisher's exact test.
    """
    input:
        bcerror=rules.base_calling_error.output.tsv,
        charging=rules.get_cca_trna.output.charging_tab,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.odds_ratios.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "odds_ratios", "{sample}"),
    params:
        src=SCRIPT_DIR,
        mod_thresh=config.get("odds_ratios", {}).get("mod_threshold", 0.3),
        ml_thresh=config.get("odds_ratios", {}).get("ml_threshold", 200),
        min_cov=config.get("odds_ratios", {}).get("min_coverage", 10),
    shell:
        """
        python {params.src}/compute_odds_ratios.py \
            --bcerror {input.bcerror} \
            --charging {input.charging} \
            --output {output.tsv} \
            --mod-threshold {params.mod_thresh} \
            --ml-threshold {params.ml_thresh} \
            --min-coverage {params.min_cov} \
            2>&1 | tee {log}
        """
