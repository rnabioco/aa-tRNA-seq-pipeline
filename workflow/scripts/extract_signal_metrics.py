import os
import sys
import argparse
import logging
from multiprocessing import Pool
from functools import partial

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


def get_metric_data(bam_fh, pod5_obj, level_table, metric, sample_name, chrom,
        strand, start, end, skip_refine_signal = False, signal_norm_method = "norm",
        scale_iters_opt = 0):
    
    s_name = sample_name
    ref_chr = chrom
    ref_strand = strand
    ref_start = start
    ref_end = end+1
    ref_reg = io.RefRegion(
        ctg=ref_chr, strand=ref_strand, start=ref_start, end=ref_end
    )

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
        signal_type=signal_norm_method
    )

    return samples_metrics,  all_bam_reads, ref_reg

def iter_metrics(samples_metrics, all_bam_reads, ref_reg):

    for samp_metrics, samp_bam_reads in zip(
        samples_metrics, all_bam_reads
    ):
        for metric, reads_metrics in samp_metrics.items():
            for bam_read, read_metrics in zip(samp_bam_reads, reads_metrics):
                for reg_pos, metric_value in enumerate(read_metrics):
                    if np.isnan(metric_value):
                        print(ref_reg.ctg,"\t",reg_pos + ref_reg.start,"\t", bam_read.query_name,"\t",metric,"\t",'NA')
                        continue
                    yield ref_reg.ctg, reg_pos + ref_reg.start, bam_read.query_name, metric, metric_value


def main(args):
    pod5_dir = os.path.abspath(args.pod5_dir)
    bam = os.path.abspath(args.bam)
    level_tab = os.path.abspath(args.level_tab)
    
    metric = args.metric
    sample_name = args.sample_name
    
    strand = args.strand

    b = pysam.AlignmentFile(bam, "rb")

    if args.chrom:
        chromosomes = [args.chrom]
    else:
        chromosomes = b.references

    start = 0
    if args.start:
        start = args.start
    
    end = None
    if args.end:
        end = args.end

    bam_obj = io.ReadIndexedBam(bam)
    pod5_obj = pod5.DatasetReader(pod5_dir)  

    print("Sample\tRead_id\tReference_Position\tMetric\tValue")

    for chrom in chromosomes:

        if end is None:
            end = b.get_reference_length(chrom)
        
        n_mapped_reads = b.count(contig = chrom, start=start, stop=end)

        if n_mapped_reads == 0:
            continue
    
        samples_metrics, all_bam_reads, ref_reg = get_metric_data(bam_obj,
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
                args.scale_iters)
        
        for metrics,reads in zip(samples_metrics, all_bam_reads):

            for i in range(len(reads)):
                read_name = reads[i].query_name,
                for metric in metrics:
                    vals = metric.values()[i, ]
                    for j in vals:
                        print(read_name,)

        sys.exit()
        for metric_row in iter_metrics(
            (sample_name,), samples_metrics, all_bam_reads, ref_reg
        ):
            print("\t".join(map(str, metric_row)))

    b.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get dwell time and signal mean for all reads (printed to stdout)"
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
        "--level_tab", type=str, help="Path to k-mers level table", required=True
    )
    parser.add_argument(
        "--metric",
        type=str,
        help="Metric that will be calculated. Possible choices:\ndwell\ndwell_mean\ndwell_mean_sd\ndwell_trimmean\ndwell_trimmean_trimsd\nDefault: dwell_mean",
        required=False,
        default="dwell_mean"
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
        "--sample_name", type=str, help="Sample name, default: sample", required=False,
        default="sample"
    )
    parser.add_argument(
        "--chrom", type=str, help="Chromosome, default: run on all chromosomes", required=False
    )
    parser.add_argument(
        "--strand", type=str, help="Strand, default: +", required=False, default="+"
    )

    parser.add_argument(
        "--start", type=int, help="Start position, default: 0", required=False, default=0
    )
    parser.add_argument(
        "--end", type=int, help="End position, default: end of chromosome", required=False
    )
    args = parser.parse_args()
    main(args)
