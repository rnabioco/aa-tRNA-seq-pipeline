"""
Rules for tRNA reference validation and building.

Ensures reference FASTA has correct adapter structure before alignment.
The CCAGGC junction (CCA from tRNA + GGC from 3' adapter) is required
for the Remora charging classification model.

Modes:
  validate: Check existing adapted reference (default)
  build: Create adapted reference from raw tRNA sequences
  skip: Copy reference as-is without validation
"""


def get_adapter_5p():
    """Get 5' adapter sequence from config with default fallback."""
    return config.get("adapters", {}).get("five_prime", "CCTAAGAGCAAGAAGAAGCCTGG")


def get_adapter_3p():
    """Get 3' adapter sequence from config with default fallback.

    For backward compatibility, returns the first adapter sequence if multiple
    are configured. Use get_adapter_3p_list() to get all adapters.
    """
    three_prime = config.get("adapters", {}).get(
        "three_prime", "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"
    )
    if isinstance(three_prime, str):
        return three_prime
    # List format - return first adapter's sequence
    return three_prime[0]["seq"]


def get_adapter_3p_list():
    """Get list of (name, seq) tuples for 3' adapters.

    Supports both string format (backward compatible) and list format
    for multiple adapter versions.

    Returns:
        List of (name, sequence) tuples
    """
    default_seq = "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"
    three_prime = config.get("adapters", {}).get("three_prime", default_seq)

    if isinstance(three_prime, str):
        return [("default", three_prime)]

    # List format: [{name: "edx01", seq: "..."}, ...]
    return [(a["name"], a["seq"]) for a in three_prime]


def get_reference_mode():
    """Get reference processing mode (validate, build, or skip)."""
    return config.get("reference", {}).get("mode", "validate")


def get_5p_offset():
    """Total 5' prefix length: adapter (23bp) + N variable position (1bp)."""
    return len(get_adapter_5p()) + 1


def get_3p_offset():
    """3' adapter length."""
    return len(get_adapter_3p())


def get_validated_reference():
    """
    Return path to validated/built reference based on mode.
    This is used by downstream rules (bwa_idx, bwa_align, etc.).
    """
    mode = get_reference_mode()
    if mode == "build":
        return os.path.join(outdir, "reference", "adapted.fa")
    elif mode == "skip":
        return os.path.join(outdir, "reference", "reference.fa")
    return os.path.join(outdir, "reference", "validated.fa")


def get_raw_reference():
    """
    Return path to raw/input reference fasta based on mode.
    For build mode, returns the raw_fasta from reference config.
    For other modes, returns the main fasta config value.
    """
    mode = get_reference_mode()
    if mode == "build":
        return config["reference"]["raw_fasta"]
    return config["fasta"]


def get_trna_fasta():
    """Return path to tRNA-only FASTA (adapters stripped)."""
    return os.path.join(outdir, "reference", "trna_only.fa")


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
        adapter_3p_args=lambda wc: " ".join(
            f'--adapter-3p "{seq}"' for _, seq in get_adapter_3p_list()
        ),
    shell:
        """
        python {params.script} \
            --mode validate \
            --input {input.fasta} \
            --output {output.validated} \
            --report {output.report} \
            --adapter-5p "{params.adapter_5p}" \
            {params.adapter_3p_args} \
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


rule skip_reference_validation:
    """
    Skip validation and copy reference as-is.

    Use this mode when the reference has non-standard adapters but you want
    to proceed without validation. The reference is copied to the output
    directory without any checks.

    WARNING: Charging classification may not work correctly if the reference
    does not have the expected CCAGGC junction structure.
    """
    input:
        fasta=config["fasta"],
    output:
        reference=os.path.join(outdir, "reference", "reference.fa"),
        report=os.path.join(outdir, "reference", "skip_report.txt"),
    log:
        os.path.join(outdir, "logs", "reference", "skip.log"),
    shell:
        """
        cp {input.fasta} {output.reference}
        echo "Reference validation skipped." > {output.report}
        echo "Input: {input.fasta}" >> {output.report}
        echo "Output: {output.reference}" >> {output.report}
        echo "WARNING: No adapter structure validation performed." >> {output.report}
        echo "Reference copied without validation." | tee {log}
        """


rule trim_reference:
    """
    Produce tRNA-only FASTA by stripping adapter sequences.

    Removes the 5' adapter + N position from the start and 3' adapter
    from the end of each reference sequence. The output is used by
    clover for MODOMICS annotation and tRNA structure visualization.

    Note: only the first 3' adapter is used here because trimming is
    length-based (fixed offset), not sequence-based. This is safe because
    validate_reference enforces that all configured 3' adapters have equal
    length. If adapters of different lengths are ever needed, this rule
    must be updated to handle per-adapter offsets.
    """
    input:
        fasta=get_validated_reference(),
    output:
        trna_fasta=os.path.join(outdir, "reference", "trna_only.fa"),
    log:
        os.path.join(outdir, "logs", "reference", "trim.log"),
    params:
        script=os.path.join(SCRIPT_DIR, "build_trna_reference.py"),
        adapter_5p=get_adapter_5p(),
        adapter_3p=get_adapter_3p(),
    shell:
        """
        python {params.script} \
            --mode trim \
            --input {input.fasta} \
            --output {output.trna_fasta} \
            --adapter-5p "{params.adapter_5p}" \
            --adapter-3p "{params.adapter_3p}" \
            2>&1 | tee {log}
        """
