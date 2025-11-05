# Issue: Add Primer/Adapter Trimming Tags to BAM Files

## Summary

Add custom BAM tags to store tRNA boundaries and adapter coordinates, enabling downstream analysis to distinguish between adapter and tRNA sequences while preserving the full alignment needed for Remora classification.

## Motivation

The pipeline aligns reads to references containing both adapters and tRNA sequences:
- 5' adapter: 24 bp
- tRNA: variable length (~70-90 bp)
- 3' adapter: 40 bp

While full alignments are necessary for Remora (which analyzes signal over the CCA-adapter junction), downstream analysis would benefit from explicit adapter boundary information. Currently, adapter coordinates must be inferred from reference structure and alignment positions.

## Proposed Solution

### Add Four Custom BAM Tags

```
ts:i:X    # tRNA start - reference coordinate where tRNA begins (0-based)
te:i:X    # tRNA end - reference coordinate where tRNA ends (0-based)
a5:i:X    # actual 5' adapter length in this alignment
a3:i:X    # actual 3' adapter length in this alignment
```

### Example

**Full-length read:**
```
Alignment: reference positions 0-150
Reference: [24bp 5'adapter][86bp tRNA][40bp 3'adapter]

Tags:
ts:i:24    # tRNA starts at ref position 24
te:i:110   # tRNA ends at ref position 110
a5:i:24    # Full 5' adapter
a3:i:40    # Full 3' adapter
```

**Truncated 5' adapter:**
```
Alignment: reference positions 5-150 (missing 5bp of 5' adapter)
Reference: [19bp 5'adapter][86bp tRNA][40bp 3'adapter]

Tags:
ts:i:24    # tRNA still starts at ref position 24
te:i:110   # tRNA ends at ref position 110
a5:i:19    # Only 19bp of 5' adapter present
a3:i:40    # Full 3' adapter
```

## Implementation Plan

### 1. Create New Script: `workflow/scripts/add_adapter_tags.py`

```python
#!/usr/bin/env python
"""Add adapter boundary tags to BAM file."""

import pysam
import argparse

def add_adapter_tags(input_bam, output_bam, adapter_5p_len, adapter_3p_len):
    """
    Add tags indicating tRNA boundaries and adapter lengths.

    Tags added:
    - ts:i - tRNA start position (reference coordinate, 0-based)
    - te:i - tRNA end position (reference coordinate, 0-based)
    - a5:i - actual 5' adapter length in alignment
    - a3:i - actual 3' adapter length in alignment
    """
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
            for read in infile:
                # Skip unmapped, secondary, supplementary
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    outfile.write(read)
                    continue

                # Get reference length for this alignment target
                ref_length = infile.get_reference_length(read.reference_name)

                # Calculate tRNA boundaries (standard positions)
                trna_start = adapter_5p_len
                trna_end = ref_length - adapter_3p_len

                # Calculate actual adapter lengths in this alignment
                actual_5p = min(read.reference_start, adapter_5p_len)
                actual_3p = min(ref_length - read.reference_end, adapter_3p_len)

                # Add tags
                read.set_tag("ts", trna_start, "i")
                read.set_tag("te", trna_end, "i")
                read.set_tag("a5", adapter_5p_len - actual_5p, "i")
                read.set_tag("a3", adapter_3p_len - actual_3p, "i")

                outfile.write(read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add adapter tags to BAM")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--adapter-5p", type=int, required=True)
    parser.add_argument("--adapter-3p", type=int, required=True)
    args = parser.parse_args()

    add_adapter_tags(args.input, args.output, args.adapter_5p, args.adapter_3p)
```

### 2. Add New Snakemake Rule

In `workflow/rules/aatrnaseq-process.smk`, add after `transfer_bam_tags`:

```python
rule add_adapter_tags:
    """
    Add tags indicating tRNA boundaries and adapter lengths
    """
    input:
        bam=rules.transfer_bam_tags.output.classified_bam,
        bai=rules.transfer_bam_tags.output.classified_bam_bai,
    output:
        tagged_bam=os.path.join(outdir, "bam", "final", "{sample}.tagged.bam"),
        tagged_bai=os.path.join(outdir, "bam", "final", "{sample}.tagged.bam.bai"),
    params:
        src=SCRIPT_DIR,
        adapter_5p_len=config["adapter_5p_length"],
        adapter_3p_len=config["adapter_3p_length"],
    log:
        os.path.join(outdir, "logs", "add_adapter_tags", "{sample}"),
    shell:
        """
        python {params.src}/add_adapter_tags.py \
          --input {input.bam} \
          --output {output.tagged_bam} \
          --adapter-5p {params.adapter_5p_len} \
          --adapter-3p {params.adapter_3p_len}

        samtools index {output.tagged_bam}
        """
```

### 3. Update Configuration

Add to `config/config-base.yml`:

```yaml
# Adapter lengths for calculating tRNA boundaries
adapter_5p_length: 24  # 5' adapter length in reference
adapter_3p_length: 40  # 3' adapter length in reference
```

### 4. Update Downstream Rules

Modify summary rules to use the new tagged BAM:
- `get_cca_trna`: Change input from `rules.transfer_bam_tags` to `rules.add_adapter_tags`
- `base_calling_error`: Change input similarly
- `align_stats`: Update classified BAM input
- `bam_to_coverage`: Update input
- `remora_signal_stats`: Update input
- `modkit_*` rules: Update inputs

### 5. Update `pipeline_outputs()` in `common.smk`

If exposing tagged BAMs as final outputs, add to the outputs list.

## Benefits

1. **Non-destructive**: Preserves alignments for Remora
2. **Explicit**: No need to infer adapter positions
3. **Handles truncation**: Stores actual adapter lengths found
4. **Flexible**: Downstream tools can use tags to extract trimmed sequences
5. **Informative**: Useful for QC and filtering

## Edge Cases Handled

- Unmapped reads: Skip tag addition
- Secondary/supplementary: Skip tag addition
- Truncated adapters: Store reduced adapter length
- Short reads: Still calculate appropriate tags

## Testing Plan

1. Run on test dataset
2. Verify all reads have expected tags
3. Check truncated reads have reduced `a5`/`a3` values
4. Ensure downstream rules still work
5. Validate tag values match expected boundaries

## Documentation Updates

- [ ] Update `CLAUDE.md` with new tags
- [ ] Update `README.md` to mention adapter tags
- [ ] Add docstring to new script
- [ ] Update workflow DAG if needed

## Questions

1. Should we also integrate `filter_reads.py` into the pipeline?
2. Should adapter lengths be reference-specific or global config?
3. Do we need additional QC metrics based on these tags?

## Alternative Considered

**Hard/soft clipping**: Rejected because it would interfere with Remora analysis and complicate the workflow.

**Separate trimmed BAM**: Rejected due to storage overhead and potential confusion.

## References

- Current filtering logic: `workflow/scripts/filter_reads.py`
- Transfer tags implementation: `workflow/scripts/transfer_tags.py`
- SAM/BAM spec: https://samtools.github.io/hts-specs/SAMv1.pdf
