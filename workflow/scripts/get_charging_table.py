#! /usr/bin/env python

"""
Generate table of read id, ref, value of charging tag, and the read's tie set.

`tRNA` is the read's primary reference -- the one the classifier scored it
against -- so the table stays one row per scored read. `tie_refs` lists the
OTHER references `escpod align` found tied with it at the best score (its `XA`
tag), `;`-separated, empty for a uniquely placed read. escpod's primary is just
the lowest-indexed reference of a tie, so counting reads by `tRNA` alone moves
13.0% total variation distance in per-reference counts against bwa on
adat2ko-pool1, where splitting each tied read evenly across its tie set brings
it to 1.8% (issue #200). get_trna_charging_cpm.py does that split, using
reference_weights() below.
"""

import argparse
import csv
import gzip
import sys

import pysam


def charging_tag_value(read, tag):
    """
    Return the scalar charging tag for one read, or None if it has none.

    Raises ValueError if the tag holds a multi-element array, which is not a
    single charging score (e.g. the dorado `ML` mod-base tag holds one
    probability per modified base) and so cannot be collapsed to one value.

    Note that a tag value of 0 is a valid, maximally-confident *uncharged*
    call, so callers must test for None rather than truthiness.
    """
    tags_dict = dict(read.tags)
    tag_raw = tags_dict.get(tag)

    # Fallback to uppercase tag for backward compat with older BAMs
    # TODO: remove fallback once all BAMs have been reprocessed
    if tag_raw is None and tag.islower():
        tag_raw = tags_dict.get(tag.upper())

    if tag_raw is None:
        return None

    if hasattr(tag_raw, "__len__") and not isinstance(tag_raw, str):
        if len(tag_raw) > 1:
            raise ValueError(f"tag '{tag}' holds a multi-element array")
        return tag_raw[0]

    return tag_raw


def tie_references(read):
    """The references tied with `read`'s primary, from its `XA` tag.

    `XA` is bwa's format, `ref,+pos,CIGAR,NM;` per entry, which is what
    `escpod align` writes -- one entry per tied reference besides the
    primary. Returned in `XA` order, de-duplicated, and without the primary
    reference itself even if an entry names it (bwa's own `XA` can list a
    second position on the same reference, which is not a second reference).
    """
    if not read.has_tag("XA"):
        return []
    primary = read.reference_name
    refs = []
    for entry in str(read.get_tag("XA")).split(";"):
        ref = entry.split(",", 1)[0].strip()
        if ref and ref != primary and ref not in refs:
            refs.append(ref)
    return refs


def reference_weights(reference, tie_refs=()):
    """{reference: weight} for one read: 1/(n+1) to each of the primary and its
    `n` tied references, so every read contributes exactly 1 in total."""
    refs = [reference] + [r for r in tie_refs if r != reference]
    refs = list(dict.fromkeys(refs))
    weight = 1.0 / len(refs)
    return dict.fromkeys(refs, weight)


def extract_tag(bam_file, output_tsv, tag):
    open_func = gzip.open if output_tsv.endswith(".gz") else open
    mode = "wt" if output_tsv.endswith(".gz") else "w"

    n_written = 0
    n_multi_skipped = 0

    with (
        pysam.AlignmentFile(bam_file, "rb") as bam,
        open_func(output_tsv, mode) as tsvfile,
    ):
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["read_id", "tRNA", "charging_likelihood", "tie_refs"])

        for read in bam.fetch():
            read_id = read.query_name
            reference = read.reference_name or "*"

            # Handle both scalar (cl:i:200) and array (CL:B:C:200) tag values.
            # Count and warn on a multi-element array rather than dropping it
            # silently: a silent drop here biases the charging fraction and CPM
            # denominator exactly like the ML==0 bug.
            try:
                tag_value = charging_tag_value(read, tag)
            except ValueError:
                n_multi_skipped += 1
                continue

            if tag_value is None:
                continue

            # Write on tag PRESENCE, not truthiness: a charging tag of 0 is a
            # valid, maximally-confident *uncharged* call (ML score range is
            # 0-255, >=200 = charged). `if tag_value` would silently drop
            # ML==0 reads, biasing charging fraction upward and shrinking the
            # CPM denominator downstream.
            if tag_value is not None and reference != "*":
                writer.writerow(
                    [read_id, reference, tag_value, ";".join(tie_references(read))]
                )
                n_written += 1

    if n_multi_skipped:
        print(
            f"WARNING: skipped {n_multi_skipped} read(s) whose '{tag}' tag was a "
            f"multi-element array (not a single charging score). If you meant to "
            f"extract the charging tag, pass --tag cl.",
            file=sys.stderr,
        )
    if n_written == 0:
        print(
            f"WARNING: wrote 0 reads to {output_tsv}; no read carried a usable "
            f"scalar '{tag}' tag.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract a specified tag from a BAM file and write to TSV."
    )
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument(
        "output_tsv", help="Output TSV file (can be .gz for compression)"
    )
    parser.add_argument("--tag", default="ML", help="BAM tag to extract (default: ML)")

    args = parser.parse_args()
    extract_tag(args.bam_file, args.output_tsv, args.tag)
