"""
Rules for amino acid identity classification using leech multiclass models.

Runs ``leech predict --bundle --all`` on the final BAMs to produce per-sample
TSV files with amino acid identity predictions. Only loaded when
``classify_aa.enabled`` is ``true`` in config.
"""


rule classify_aa:
    """
    Run leech multiclass bundle on final BAMs to predict amino acid identity.

    Uses reference-anchored signal extraction over the CCA 3' end junction.
    Outputs per-read predictions as gzipped TSV with probability columns for
    each amino acid class, plus the CL (charging) tag from the BAM.
    """
    input:
        bam=rules.finalize_bam.output.bam,
        bai=rules.finalize_bam.output.bai,
        pod5=get_classification_pod5,
        ref=get_validated_reference(),
    output:
        predictions=os.path.join(
            outdir,
            "summary",
            "tables",
            "{sample}",
            "{sample}.aa_classify.tsv.gz",
        ),
    log:
        os.path.join(outdir, "logs", "classify_aa", "{sample}"),
    params:
        bundle=config["classify_aa"]["bundle"],
        device=config.get("classify_aa", {}).get("device", "cuda"),
        batch_size=config.get("classify_aa", {}).get("batch_size", 1024),
        base_justify=config.get("classify_aa", {}).get("base_justify", "end"),
    shell:
        """
        leech predict \
            --bundle {params.bundle} \
            --all \
            --reference-anchored \
            --base-justify {params.base_justify} \
            --pod5 {input.pod5} \
            --bam {input.bam} \
            --output {output.predictions} \
            --device {params.device} \
            --batch-size {params.batch_size} \
            --workers 0 \
            2>&1 | stdbuf -oL -eL tee {log}
        """
