# Log Summary

This file summarizes the status of the main Snakemake rules based on their log files.

## fastqc

- **Log file**: `logs/fastqc_NA12878.log`
- **Status**: OK (FastQC completed and produced HTML + ZIP outputs)
- **Notes**: Raw FASTQ file `data/raw/NA12878.fastq.gz` was successfully processed.

## index_reference

- **Log file**: `logs/index_reference.log`
- **Status**: OK
- **Notes**: BWA index files (`chr20.fa.*`) were created for `data/ref/chr20.fa`.

## bwa_map

- **Log file**: `logs/bwa_NA12878.log`
- **Status**: OK
- **Notes**: Reads from `NA12878.fastq.gz` were aligned to `chr20.fa` and written to `results/mapping/NA12878.bam`.

## sort_bam

- **Log file**: `logs/sort_NA12878.log`
- **Status**: OK
- **Notes**: `NA12878.bam` was successfully sorted into `NA12878.sorted.bam`.

## index_bam

- **Log file**: `logs/index_NA12878.log`
- **Status**: OK
- **Notes**: `NA12878.sorted.bam.bai` index was created.

## call_variants

- **Log file**: `logs/calling_NA12878.log`
- **Status**: OK
- **Notes**: `bcftools mpileup` and `bcftools call` completed. Outputs:
  - `results/variants/NA12878.bcf`
  - `results/variants/NA12878.vcf.gz`
  - `results/variants/NA12878.vcf.gz.csi`

