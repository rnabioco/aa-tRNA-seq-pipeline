# Rules for per-tRNA pairwise modification odds ratio analysis
# Computes odds ratios between all pairs of positions (+ charging) within each tRNA


rule compute_odds_ratios:
    """
    Compute per-tRNA pairwise modification odds ratios.

    For each tRNA, uses individual reads as the unit of observation to test
    whether modification at position X is correlated with modification at
    position Y (and with charging status) via 2x2 contingency tables,
    odds ratios, and Fisher's exact test.
    """
    input:
        modkit=rules.modkit_extract_calls.output.tsv,
        charging=rules.get_cca_trna.output.charging_tab,
    output:
        tsv=os.path.join(
            outdir, "summary", "tables", "{sample}", "{sample}.odds_ratios.tsv.gz"
        ),
    log:
        os.path.join(outdir, "logs", "odds_ratios", "{sample}"),
    params:
        src=SCRIPT_DIR,
        ml_thresh=config.get("odds_ratios", {}).get("ml_threshold", 200),
        min_cov=config.get("odds_ratios", {}).get("min_coverage", 10),
    shell:
        """
        python {params.src}/compute_odds_ratios.py \
            --modkit {input.modkit} \
            --charging {input.charging} \
            --output {output.tsv} \
            --ml-threshold {params.ml_thresh} \
            --min-coverage {params.min_cov} \
            2>&1 | tee {log}
        """
