# Configuring the pipeline with config.yml

Edit config.yml to specify the following parameters. 

## samples
Provide a path to a TSV file which indicates the samples to process (e.g. `config/samples.tsv`). The 
TSV file should have two columns containing the following information.

  - a unique id for the sample
  - a path to the sequencing run folder which has `pod5_pass` or `fast5_pass` subdirectories containing raw data 
  OR a path to an unmapped BAM file. **If an unmapped BAM file is supplied, rebasecalling will not be performed.**

The pipeline will recursively search for pod5 or fast5 files to process within the specified directory and subdirectories 

## output_directory
A path to an output directory for files produced by pipeline

## rebasecalled_bam_directory
A path to an output directory where the rebasecalled BAM files and merged pod5 files will be written

## base_calling_model
A path to the dorado basecalling model to use for rebasecalling

## input_format
A string, either "FAST5" or "POD5", if FAST5 then these files will be converted to pod5 before rebasecalling

## fasta
A path to the reference fasta file to use for bwa alignment. A BWA index will be built if it does not exist for this fasta file

## Additional options
These should not be edited at this time. 