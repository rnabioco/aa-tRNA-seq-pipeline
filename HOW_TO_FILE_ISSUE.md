# How to File the Primer Trimming Tags Issue on GitHub

Since the GitHub CLI is not available in this environment, please file the issue manually using the content in `ISSUE_PRIMER_TRIMMING_TAGS.md`.

## Steps to File Issue

1. **Navigate to the repository**: https://github.com/rnabioco/aa-tRNA-seq-pipeline/issues

2. **Click "New Issue"**

3. **Title**: `Add Primer/Adapter Trimming Tags to BAM Files`

4. **Labels to add** (if available):
   - `enhancement`
   - `pipeline`
   - `documentation`

5. **Body**: Copy the entire content of `ISSUE_PRIMER_TRIMMING_TAGS.md`

## Quick Summary for Issue Description

If you prefer a shorter issue, here's a condensed version:

---

### Title
Add Primer/Adapter Trimming Tags to BAM Files

### Description

**Problem**: The pipeline aligns reads to references containing adapters, but adapter boundaries are not explicitly stored in BAM tags, requiring downstream tools to infer them from reference structure.

**Proposed Solution**: Add four custom BAM tags to store tRNA boundaries and adapter information:
- `ts:i` - tRNA start position (reference coordinate)
- `te:i` - tRNA end position (reference coordinate)
- `a5:i` - actual 5' adapter length in alignment
- `a3:i` - actual 3' adapter length in alignment

**Benefits**:
- Non-destructive (preserves alignments needed for Remora)
- Explicitly stores adapter boundaries
- Handles truncated adapters correctly
- Enables easier downstream analysis

**Implementation**: See `ISSUE_PRIMER_TRIMMING_TAGS.md` for detailed implementation plan including:
- New script: `workflow/scripts/add_adapter_tags.py`
- New rule: `add_adapter_tags`
- Config updates: Add `adapter_5p_length` and `adapter_3p_length` parameters
- Update downstream rules to use tagged BAM

See `ISSUE_PRIMER_TRIMMING_TAGS.md` and `.github-issue-primer-trimming-tags.md` for complete analysis and implementation details.

---

## Files in This Branch

- **ISSUE_PRIMER_TRIMMING_TAGS.md**: Complete issue description with implementation plan (use this as the issue body)
- **.github-issue-primer-trimming-tags.md**: Detailed technical analysis of different approaches
- **HOW_TO_FILE_ISSUE.md**: This file

## Next Steps After Filing Issue

1. Get feedback from maintainers on the proposed approach
2. Discuss any concerns about tag naming or format
3. Confirm adapter length configuration approach
4. Implement the solution once approved
5. Add tests and documentation
