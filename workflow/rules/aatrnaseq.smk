
rule rebasecall:
  """
  rebasecall using different accuracy model
  requires a GPU
  """
  input:
    get_basecalling_inputs
  output:
    os.path.join(rbc_outdir, "{sample}", "{sample}.unmapped.bam")
  log:
    os.path.join(outdir, "logs", "rebasecall", "{sample}")
  params:
    model = config["base_calling_model"],
    is_fast5 = config["input_format"],
    raw_data_dir = get_basecalling_dir,
    temp_pod5 = os.path.join(rbc_outdir, "{sample}", "{sample}.pod5") 
  shell:
    """
    if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
      echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
      export CUDA_VISIBLE_DEVICES 
    fi

    if [ "{params.is_fast5}" == "FAST5" ]; then
      pod5 convert fast5 {params.raw_data_dir} --output {params.temp_pod5}
      od=$(dirname {params.temp_pod5})
      dorado basecaller --emit-moves -v {params.model} $od > {output} 
      rm {params.temp_pod5}
    else 
      dorado basecaller --emit-moves -v {params.model} {params.raw_data_dir} > {output}
    fi
    """


rule ubam_to_fq:
  """
  extract reads from bam into FASTQ format for alignment
  """
  input:
    rules.rebasecall.output
  output:
    os.path.join(outdir, "fastqs", "{sample}.fastq.gz")
  log:
    os.path.join(outdir, "logs", "ubam_to_fq", "{sample}") 
  shell:
    """
    samtools fastq {input} | gzip > {output}
    """

# note that the mapped input for merge_bams is conditional on config["aligner"] settings
# either parasail.smk or bwa.smk are used to generate the mapped bam file.

rule merge_bams:
  """
  merge move table data from unmapped bam into mapped bam 
  note that output bam will have .sorted.bam appended to 
  the  output bam file listed in the "output" params
  """
  input:
    unmapped = rules.rebasecall.output,
    mapped = os.path.join(outdir, "bams", "{sample}", "{sample}.{aligner}.bam")
  output:
    bam = os.path.join(outdir, "bams", "{sample}", "{sample}.{aligner}.merged.sorted.bam"),
  log:
    os.path.join(outdir, "logs", "merge_bams", "{sample}.{aligner}") 
  params:
    src = config["src"]
  shell:
    """
    python {params.src}/merge_bams.py \
        --ubam {input.unmapped} \
        --mbam {input.mapped} \
        --merged {output.bam}
    """

rule calc_samples_per_base:
  """
  calculate samples per base metric for all mapped positions
  """
  input:
    rules.merge_bams.output,
  output:
    os.path.join(outdir, "bams", "{sample}", "{sample}.{aligner}.merged.sb.bam"),
  log:
    os.path.join(outdir, "logs", "calc_samples_per_base", "{sample}.{aligner}")  
  params:
    root =  os.path.join(outdir, "bams", "{sample}"), 
    input_bam = "{sample}.{aligner}.merged.sorted.bam",
    output_bam = "{sample}.{aligner}.merged.sb.bam", 
    src = config["src"]
  shell:
    """
    python {params.src}/add_sb_tags.py \
        --input {input} \
        --output {output}
    """

rule extract_sb_tag:
  """
  extract metrics to tsv file
  """
  input:
    expand(os.path.join(outdir, "bams", "{sample}", "{sample}.{aligner}.merged.sb.bam"),
      sample = samples,
      aligner = config["aligner"])  
  output:
    os.path.join(outdir, "tables", "sb_values.tsv"),
  params:
    src = config["src"],
  log:
    os.path.join(outdir, "logs", "extract_sb_tag", "sb") 
  threads:
    4
  shell:
    """
    python {params.src}/process_sb_with_insertions.py \
        --output {output} \
        --proc {threads} \
        {input} 
    """


rule bcerror:
  """
  extract base calling error  metrics to tsv file
  """
  input:
    rules.calc_samples_per_base.output  
  output:
    tsv = os.path.join(outdir, "tables", "{sample}.{aligner}.bcerror.tsv"), 
  log:
    os.path.join(outdir, "logs", "bcerror", "{sample}.{aligner}") 
  params:
    src = config["src"],
    fa  = config["fasta"]
  shell:
    """
    python {params.src}/get_bcerror_freqs.py \
      {input} \
      {params.fa} \
      {output}
    """

