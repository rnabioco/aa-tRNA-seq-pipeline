rule parasail:
  """
  align reads to tRNA references with parasail routine 
  """
  input:
    reads = rules.ubam_to_fq.output,
  output:
    bam = os.path.join(outdir, "bams", "{sample}", "{sample}.bam"),
    bai = os.path.join(outdir, "bams", "{sample}", "{sample}.bam.bai"),
  params:
    index = config["fasta"],
    src = config["src"],
    outpre = os.path.join(outdir, "bams", "{sample}", "{sample}")
  log:
    os.path.join(outdir, "logs", "parasail", "{sample}") 
  threads: 4
  shell:
    """
    python {params.src}/parasail.py \
        -p 0.95 \
        -t {threads} \
        -o {params.outpre} \
        {params.index} \
        {input.reads}
    """

