import os
import sys
import argparse
import logging

import pod5
import pysam
import numpy as np
import polars as pl

import remora
from remora import io, refine_signal_map, util

#######
# Requires forked version of remora, using the "metrics_missing_ok" branch
# https://github.com/rnabioco/remora/tree/metrics_missing_ok
# to install
# pip install git+https://github.com/rnabioco/remora.git@metrics_missing_ok
#######

# silence Remora DEBUG messages
# logging.getLogger("Remora").setLevel(logging.DEBUG)


def get_metric_data(
    bam_fh,
    pod5_obj,
    level_table,
    metric,
    sample_name,
    chrom,
    strand,
    start,
    end,
    skip_refine_signal=False,
    signal_norm_method="norm",
    scale_iters_opt=0,
):
    s_name = sample_name
    ref_chr = chrom
    ref_strand = strand
    ref_start = start  # start is 0-based
    ref_end = end  # end is half-open
    ref_reg = io.RefRegion(ctg=ref_chr, strand=ref_strand, start=ref_start, end=ref_end)

    if skip_refine_signal:
        sig_map_refiner = None
    else:
        sig_map_refiner = refine_signal_map.SigMapRefiner(
            kmer_model_filename=level_table,
            do_rough_rescale=True,
            scale_iters=scale_iters_opt,
            do_fix_guage=True,
        )

    samples_metrics, all_bam_reads = io.get_ref_reg_samples_metrics(
        ref_reg,
        [(pod5_obj, bam_fh)],
        metric=metric,
        sample_names=[s_name],
        sig_map_refiner=sig_map_refiner,
        reverse_signal=True,
        missing_ok=True,
        signal_type=signal_norm_method,
    )

    return samples_metrics, all_bam_reads, ref_reg


def iter_metrics(sample_name, samples_metrics, all_bam_reads, ref_reg):
    for metrics, reads in zip(samples_metrics, all_bam_reads):
        n_reads, n_pos = next(iter(metrics.values())).shape

        for i in range(n_reads):
            read_name = reads[i].query_name
            for reg_pos in range(n_pos):
                row = (sample_name, ref_reg.ctg, reg_pos + ref_reg.start + 1, read_name)
                for metric in metrics.values():
                    val = metric[i, reg_pos]
                    if np.isnan(val):
                        row += ("NA",)
                        continue
                    row += (val,)
                yield row


def window_regions(regions, window_length):
    win_regions = []
    for chrom, start, end, strand in regions:
        for i in range(start, end, window_length):
            win_regions.append((chrom, i, min(i + window_length, end), strand))
    return win_regions


def get_regions_to_query(bam_fh, bed_file=None, region=None, window_length=100):
    regions = []
    if region:
        strand = "+"
        if region.endswith("+") or region.endswith("-"):
            strand = region[-1]
            region = region[:-2]
        _, tid, start, end = bam_fh.parse_region(region=region)
        chrom = bam_fh.get_reference_name(tid)
        regions.append((chrom, start, end, strand))

    elif bed_file:
        with open(bed_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                chrom, start, end = fields[:3]

                if len(fields) >= 6:
                    strand = fields[5]
                else:
                    strand = "+"

                regions.append((chrom, int(start), int(end), strand))
    else:
        for chrom in bam_fh.references:
            regions.append((chrom, 0, bam_fh.get_reference_length(chrom), "+"))
    if window_length > 0:
        regions = window_regions(regions, window_length)

    return regions


def main(args):
    pod5_dir = os.path.abspath(args.pod5_dir)
    bam = os.path.abspath(args.bam)
    level_tab = os.path.abspath(args.kmer)

    metric = args.metric
    sample_name = args.sample_name

    strand = "+"

    b = pysam.AlignmentFile(bam, "rb")

    regions = get_regions_to_query(b, args.bed, args.region, args.window)

    bam_obj = io.ReadIndexedBam(bam)
    pod5_obj = pod5.DatasetReader(pod5_dir)

    print("Sample\tContig\tReference_Position\tRead_id\t", "\t".join(metric.split("_")))

    for region in regions:
        chrom, start, end, strand = region

        n_mapped_reads = b.count(contig=chrom, start=start, stop=end)

        if n_mapped_reads == 0:
            continue

        try:
            samples_metrics, all_bam_reads, ref_reg = get_metric_data(
                bam_obj,
                pod5_obj,
                level_tab,
                metric,
                sample_name,
                chrom,
                strand,
                start,
                end,
                args.skip_refine_signal,
                args.signal_norm,
                args.scale_iters,
            )

        except Exception as e:
            print(f"Unable to process region {region}: {e}", file=sys.stderr)
            continue

        for metric_row in iter_metrics(
            sample_name, samples_metrics, all_bam_reads, ref_reg
        ):
            print("\t".join(map(str, metric_row)))

    b.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Extract signal metrics using the Remora API. Output will be generated
        across all regions with read coverage, restricting to reads mapped on 
        the positive strand, as it is expected that the reads are aligned 
        against transcripts rather than a genome reference. Output can be restricted
        to a specific region using the --region option, or to a set of regions
        using the --bed option. 
        The output is TSV text with the following columns:
        Sample\tContig\tReference_Position\tRead_id\tMetric1\tMetric2\tMetric3\t...
        The Reference_Position is 1-based.
        """
    )
    parser.add_argument(
        "--pod5_dir",
        type=str,
        help="Directory, where merged pod5 is stored (only one file per dir)",
        required=True,
    )
    parser.add_argument(
        "--bam", type=str, help="Path to mapped bam file", required=True
    )
    parser.add_argument(
        "--kmer", type=str, help="Path to k-mers level table", required=True
    )
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=0,
        help="""
        Window size used for signal extraction. Regions will be chunked into windows of this length prior
        to processing. Use this option if you want to extract data from large regions (e.g. regions >> than the read length).
        Without this option the entire region will be processed at once, which for e.g. chromosomes or long RNAs 
        would use excessive memory. Set this to the median of the read lengths in the dataset. Setting to 0 disables this
        option, which is the default. Default: 0
        """,
    )
    parser.add_argument(
        "--metric",
        type=str,
        help="Metric that will be calculated. Possible choices:\ndwell\ndwell_mean\ndwell_mean_sd\ndwell_trimmean\ndwell_trimmean_trimsd\nDefault: dwell_trimmean_trimsd",
        required=False,
        default="dwell_trimmean_trimsd",
    )
    parser.add_argument(
        "--signal_norm",
        type=str,
        help="Signal type, one of norm, pa, and dac\nDefault: norm",
        required=False,
        default="norm",
    )
    parser.add_argument(
        "--skip_refine_signal",
        help="If set, skip signal refinement",
        action="store_true",
    )
    parser.add_argument(
        "--scale_iters",
        help="scale iters arg for signal refinement",
        type=int,
        default=0,
        required=False,
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        help="Sample name, default: sample",
        required=False,
        default="sample",
    )
    parser.add_argument(
        "--region",
        type=str,
        help="""
        If supplied, operate only on the specified region. Uses samtools style region string (e.g. tRNA:1-20)
        A strand can be provided using the format tRNA:1-20:+ or tRNA:1-20:-, and will be "+" if not supplied.
        """,
        required=False,
    )

    parser.add_argument(
        "--bed",
        help="""
        If supplied, produce output from the regions specified in supplied bed file
        A strand can specified by including the strand in the 6th column.
        """,
        required=False,
    )

    args = parser.parse_args()
    main(args)
