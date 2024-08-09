
rule merge_pods:
  """
  merge all fast5/pod5s into a single pod5
  """
  input:
    get_raw_inputs
  output:
    os.path.join(rbc_outdir, "{sample}", "{sample}.pod5")
  log:
    os.path.join(outdir, "logs", "merge_pods", "{sample}")
  params:
    is_fast5 = config["input_format"],
  shell:
    """
    if [ "{params.is_fast5}" == "FAST5" ]; then
      pod5 convert fast5 -f --output {output} {input}
    else
      pod5 merge -f -o {output} {input} 
    fi
    """

rule rebasecall:
  """
  rebasecall using different accuracy model
  requires a GPU
  """
  input:
    rules.merge_pods.output
  output:
    protected(os.path.join(rbc_outdir, "{sample}", "{sample}.unmapped.bam"))
  log:
    os.path.join(outdir, "logs", "rebasecall", "{sample}")
  params:
    model = config["base_calling_model"],
    is_fast5 = config["input_format"],
    raw_data_dir = get_basecalling_dir,
    temp_pod5 = os.path.join(rbc_outdir, "{sample}", "{sample}.pod5"),
    dorado_opts = config["opts"]["dorado"], 
  shell:
    """
    if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
      echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
      export CUDA_VISIBLE_DEVICES 
    fi

    dorado basecaller {params.dorado_opts} -v {params.model} {input} > {output}
    """

def get_optional_bam_inputs(wildcards):
  sample = wildcards.sample

  if config["input_format"] == "BAM":
    return samples[sample]["raw_files"]
  else:
    return os.path.join(rbc_outdir, sample, sample + ".unmapped.bam") 

rule ubam_to_fq:
  """
  extract reads from bam into FASTQ format for alignment
  """
  input:
    get_optional_bam_inputs
  output:
    os.path.join(outdir, "fastqs", "{sample}.fastq.gz")
  log:
    os.path.join(outdir, "logs", "ubam_to_fq", "{sample}") 
  shell:
    """
    samtools fastq -T "*" {input} | gzip > {output}
    """

rule bwa_idx:
  input:
    config["fasta"]
  output:
    multiext(config["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa")
  log:
    os.path.join(outdir, "logs", "bwa_idx", "log") 
  shell:
    """
    bwa index {input}
    """ 
    
rule bwa:
  """
  align reads to tRNA references with bwa mem
  """
  input:
    reads = rules.ubam_to_fq.output,
    idx = rules.bwa_idx.output
  output:
    bam = os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.unfiltered.bam"),
    bai = os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.unfiltered.bam.bai"),
  params:
    index = config["fasta"],
    bwa_opts = config["opts"]["bwa"],
  log:
    os.path.join(outdir, "logs", "bwa", "{sample}") 
  threads: 12 
  shell:
    """
    bwa mem -C -t {threads} {params.bwa_opts} {params.index} {input.reads} \
        | samtools sort -o {output.bam}

    samtools index {output.bam}
    """

def trna_filter_opts():
  trna_table = config["trna_table"]
  if trna_table is not None and trna_table != "":
    opts = f" -r -t {config['trna_table']}"
  else:
    opts = ""
  return opts

rule filter_bwa:
  """
  filter bwa mem reads
  """
  input:
    reads = rules.bwa.output.bam,
  output:
    bam = os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.bam"),
    bai = os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.bam.bai"),
    failed_bam = os.path.join(outdir, "bams", "{sample}", "{sample}.bwa.failed.bam")
  params:
    src = SCRIPT_DIR,
    trna_table = trna_filter_opts(),
    bf_opts = config["opts"]["bam_filter"],
  log:
    os.path.join(outdir, "logs", "bwa", "{sample}_filter") 
  shell:
    """
    python {params.src}/filter_reads.py \
      {params.bf_opts} \
      -i {input.reads} \
      -o {output.bam} \
      -f {output.failed_bam} \
      {params.trna_table} 

      
    samtools index {output.bam}
    samtools index {output.failed_bam}
    """

rule bcerror:
  """
  extract base calling error metrics to tsv file
  """
  input:
    bam = rules.filter_bwa.output.bam, 
    bai = rules.filter_bwa.output.bai  
  output:
    tsv = os.path.join(outdir, "tables", "{sample}.bwa.bcerror.tsv"), 
  log:
    os.path.join(outdir, "logs", "bcerror", "{sample}.bwa") 
  params:
    src = SCRIPT_DIR,
    fa  = config["fasta"]
  shell:
    """
    python {params.src}/get_bcerror_freqs.py \
      {input.bam} \
      {params.fa} \
      {output}
    """

rule align_stats:
  """
  extract alignment stats
  """
  input:
    unmapped = get_optional_bam_inputs,
    unfiltered = rules.bwa.output.bam,
    mapped = rules.filter_bwa.output.bam 
  output:
    tsv = os.path.join(outdir, "tables", "{sample}.bwa.align_stats.tsv"),
  log:
    os.path.join(outdir, "logs", "stats", "{sample}.bwa") 
  params:
    src = SCRIPT_DIR
  shell:
    """
    python {params.src}/get_align_stats.py \
      -o {output.tsv} \
      -a unmapped unfiltered filtered \
      -i {wildcards.sample} \
      -b {input.unmapped} \
         {input.unfiltered} \
         {input.mapped} 
    """

rule combine_sample_stats:
  """
  extract alignment stats
  """
  input:
    expand(os.path.join(outdir, "tables", "{sample}.bwa.align_stats.tsv"),
        sample = samples.keys())   
  output:
    os.path.join(outdir, "tables", "align_stats.tsv"),
  log:
    os.path.join(outdir, "logs", "combine_stats", "combine_stats") 
  shell:
    """
    # don't combine the headers from every file
    head -n 1 {input[0]} > {output}
    tail -n +2 -q {input} >> {output}
    """

rule bam_to_coverage:
  input:
    bam = rules.filter_bwa.output.bam,
    bai = rules.filter_bwa.output.bai,
  output:
    counts = os.path.join(outdir, "tables", "{sample}.bwa.counts.bg"),
    cpm = os.path.join(outdir, "tables", "{sample}.bwa.cpm.bg")
  params:
    bg_opts = config["opts"]["coverage"]
  log:
   os.path.join(outdir, "logs", "bg", "{sample}.bwa.txt")
  threads: 4 
  shell:
    """
    bamCoverage \
      -b {input.bam} \
      -o {output.cpm} \
      --normalizeUsing CPM \
      --outFileFormat bedgraph \
      -bs 1 \
      -p {threads} \
      {params.bg_opts}

    bamCoverage \
      -b {input.bam} \
      -o {output.counts} \
      --outFileFormat bedgraph \
      -bs 1 \
      -p {threads} \
      {params.bg_opts}

    """

rule remora:
  """
  run remora to get signal stats
  """
  input:
    bam = rules.filter_bwa.output.bam,
    bai = rules.filter_bwa.output.bai,
    pod5 = rules.merge_pods.output
  output:
    tsv = os.path.join(outdir, "tables", "{sample}.bwa.remora.tsv.gz"),
  log:
    os.path.join(outdir, "logs", "remora", "{sample}") 
  params:
    src = SCRIPT_DIR,
    kmer = config["remora_kmer_table"],
    opts = config["opts"]["remora"],
  shell:
    """
    python {params.src}/extract_signal_metrics.py \
      --pod5_dir {input.pod5} \
      --bam {input.bam} \
      --kmer {params.kmer} \
      --sample_name {wildcards.sample} \
      {params.opts} \
      | gzip > {output}
    """