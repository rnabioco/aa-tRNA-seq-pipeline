# Configuring the pipeline with config.yml

Edit config.yml to specify the following parameters. 

- `samples`: Provide a path to a TSV file which indicates the samples to process (e.g. `config/samples.tsv`). The two-column
  TSV file should have (1) a unique id for the sample and (2 a path to the sequencing run folder which has `pod5_pass`, `pod5`, `pod5_fail`, `fast5_pass`, or `fast5_fail` subdirectories containing raw data.
  The pipeline will recursively search for POD5 files to process within the specified directory.

  - a unique id for the sample
  - a path to the sequencing run folder with `pod5_pass` and `pod5_fail` subdirectories containing raw data.

- `base_calling_model`: Path to the dorado basecalling model to use for rebasecalling. We use `rna004_130bps_sup@v5.0.0` for now, will evaluate newer model soon.

- `input_format`: A string, either "FAST5" or "POD5", if FAST5 then these files will be converted to pod5 before rebasecalling

- `fasta`: A path to the reference fasta file to use for bwa alignment. A BWA index will be built if it does not exist for this fasta file

- `remora_kmer_table`: Path to a table of expected normalized signal intensites for each kmer, provided by ONT at [nanoporetech/kmer_models](https://github.com/nanoporetech/kmer_models)  

- `trna_table`: A path to table with trna isodecoder + sequencing adapter annotation from the fasta reference file.

  The format is four whitespace deliminated columns, no header: 
    - sequence name of uncharged tRNA: Name of the uncharged tRNA sequence, must match the fasta entry (e.g. tRNA-Ala-AGC-1-1-uncharged)
    - sequence name of charged tRNA: Name of the charged tRNA sequence, must match the fasta entry (e.g. tRNA-Ala-AGC-1-1-charged)
    - isodecoder: The isodecoder family of the tRNA sequence (e.g. Ala-AGC)
    - tRNA gene name: tRNA gene name that is used to represent the tRNA, can be any string, and doesn't have to match a sequence name in the fasta file.

  As this pipeline current performs charging classification at the level of nanopore current signal (using Remora) rather than adapter sequences, this table is currently optional, but may be useful for filtering based on additional alignment-level information.

- `opts`: Customized options for some commands. Note that the bam_filter option controls full-length read filtering options. 
