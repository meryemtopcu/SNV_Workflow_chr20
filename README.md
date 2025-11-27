# SNV Calling Workflow (NA12878 – chr20)

This project implements a minimal Snakemake workflow for single nucleotide variant (SNV) calling using short-read data.

## Data

- Sample: NA12878 (subset, single FASTQ)
- Reference: `data/ref/chr20.fa` (toy chr20 sequence)
- Raw reads: `data/raw/NA12878.fastq.gz`

## Workflow

Main steps:

1. **QC (FastQC)**
   - Input: `data/raw/NA12878.fastq.gz`
   - Output: `results/qc/NA12878_fastqc.html`, `results/qc/NA12878_fastqc.zip`
   - Rule file: `workflow/rules/qc.smk`

2. **Mapping (BWA + samtools)**
   - Input: `data/raw/NA12878.fastq.gz`, `data/ref/chr20.fa`
   - Output:
     - `results/mapping/NA12878.sorted.bam`
     - `results/mapping/NA12878.sorted.bam.bai`
   - Rule file: `workflow/rules/mapping.smk`

3. **Variant Calling (bcftools)**
   - Input:
     - `results/mapping/NA12878.sorted.bam`
     - `results/mapping/NA12878.sorted.bam.bai`
     - `data/ref/chr20.fa`
   - Output:
     - `results/variants/NA12878.bcf`
     - `results/variants/NA12878.vcf.gz`
     - `results/variants/NA12878.vcf.gz.csi`
   - Rule file: `workflow/rules/calling.smk`

## Snakemake Concepts

- **Wildcards**: `{sample}` for samples (e.g. `NA12878`).
- **Config file**: `config/config.yaml` stores sample list and paths.
- **Conda environment**: `envs/snv_env.yaml` defines all tools (snakemake, bwa, samtools, bcftools, fastqc).
- **DAG**: generated via  
  `snakemake -n --dag | dot -Tpng > workflow_dag.png`
- **Logs**: per-step logs in `logs/` (e.g. `logs/bwa_NA12878.log`, `logs/calling_NA12878.log`).
