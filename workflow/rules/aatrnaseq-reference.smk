"""
Rules for tRNA reference validation and building.

Ensures reference FASTA has correct adapter structure before alignment.
The CCAGGC junction (CCA from tRNA + GGC from 3' adapter) is required
for the Remora charging classification model.

Modes:
  validate: Check existing adapted reference (default)
  build: Create adapted reference from raw tRNA sequences
"""


def get_adapter_5p():
    """Get 5' adapter sequence from config with default fallback."""
    return config.get("adapters", {}).get("five_prime", "CCTAAGAGCAAGAAGAAGCCTGG")


def get_adapter_3p():
    """Get 3' adapter sequence from config with default fallback."""
    return config.get("adapters", {}).get(
        "three_prime", "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"
    )


def get_reference_mode():
    """Get reference processing mode (validate or build)."""
    return config.get("reference", {}).get("mode", "validate")


def get_validated_reference():
    """
    Return path to validated/built reference based on mode.
    This is used by downstream rules (bwa_idx, bwa_align, etc.).
    """
    mode = get_reference_mode()
    if mode == "build":
        return os.path.join(outdir, "reference", "adapted.fa")
    return os.path.join(outdir, "reference", "validated.fa")


rule validate_reference:
    """
    Validate that an existing reference FASTA has correct adapter structure.

    Checks:
    - All sequences have correct 5' adapter prefix
    - All tRNA portions end with CCA
    - All sequences have correct 3' adapter suffix (starting with GGC)
    - CCAGGC junction exists for charging classification
    - No duplicate sequence names

    Pipeline fails if validation fails.
    """
    input:
        fasta=config["fasta"],
    output:
        validated=os.path.join(outdir, "reference", "validated.fa"),
        report=os.path.join(outdir, "reference", "validation_report.txt"),
    log:
        os.path.join(outdir, "logs", "reference", "validate.log"),
    params:
        script=os.path.join(SCRIPT_DIR, "build_trna_reference.py"),
        adapter_5p=get_adapter_5p(),
        adapter_3p=get_adapter_3p(),
    shell:
        """
        python {params.script} \
            --mode validate \
            --input {input.fasta} \
            --output {output.validated} \
            --report {output.report} \
            --adapter-5p "{params.adapter_5p}" \
            --adapter-3p "{params.adapter_3p}" \
            2>&1 | tee {log}
        """


rule build_reference:
    """
    Build an adapted tRNA reference from raw tRNA sequences.

    Input: Raw tRNA FASTA (sequences without adapters)
    Output: Adapted reference FASTA with 5' and 3' adapters

    Steps:
    1. Check for CCA endings - add CCA if missing (with warning)
    2. Prepend 5' adapter to each tRNA
    3. Append 3' adapter after CCA
    4. Verify CCAGGC junction is created
    5. Write adapted FASTA
    """
    input:
        raw_fasta=lambda wildcards: config["reference"]["raw_fasta"],
    output:
        adapted=os.path.join(outdir, "reference", "adapted.fa"),
        report=os.path.join(outdir, "reference", "build_report.txt"),
    log:
        os.path.join(outdir, "logs", "reference", "build.log"),
    params:
        script=os.path.join(SCRIPT_DIR, "build_trna_reference.py"),
        adapter_5p=get_adapter_5p(),
        adapter_3p=get_adapter_3p(),
    shell:
        """
        python {params.script} \
            --mode build \
            --input {input.raw_fasta} \
            --output {output.adapted} \
            --report {output.report} \
            --adapter-5p "{params.adapter_5p}" \
            --adapter-3p "{params.adapter_3p}" \
            2>&1 | tee {log}
        """
