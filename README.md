# Project: WES Tumor-Only Analysis Pipeline (Breast Cancer Case Study)

## Language Selection / Выбор языка

[English](#readme-english) | [Русский](#readme-русский)

---

## **README (English)**

### Project Overview
Tumor-Only WES Analysis Pipeline for Breast Cancer

#Description
This repository contains a complete **whole exome sequencing (WES)** pipeline for tumor-only analysis, specifically designed for breast cancer. The pipeline includes multiple steps, from raw data download to clinical report generation.

This analysis uses tools like **GATK**, **BWA-MEM2**, **bcftools**, **VEP**, and others, to process raw NGS data, call mutations, and provide clinically relevant results.

How to Run the Pipeline

Prerequisites
Before running the pipeline, ensure you have the following installed:
- `GATK`
- `BWA-MEM2`
- `bcftools`
- `VEP`
- `samtools`
- `fastp`
- `fastqc`
- `python3` (optional, if you need to use custom scripts)

Step 1: Download data from SRA


prefetch SRP291993

Step 2: Convert SRA to FASTQ Format

fastq-dump --split-files --gzip ~/1/SRR13018652.sra

Explanation:
Converts the SRA file (from Step 1) into FASTQ format. The --split-files flag ensures that paired-end data is separated into two files. The --gzip flag compresses the output files to save space.

Step 3: Decompress the FASTQ Files

gunzip ~/1/SRR13018652_1.fastq.gz

gunzip ~/1/SRR13018652_2.fastq.gz

Explanation:
Decompress the gzipped FASTQ files so they can be used in the next steps. The two FASTQ files correspond to the paired-end reads.

Step 4: Quality Control (QC) Before Trimming

fastqc ~/1/SRR13018652_1.fastq ~/1/SRR13018652_2.fastq -o ~/1/fastqc_reports

Explanation:
Runs FastQC on the raw FASTQ files to assess their quality. This generates an HTML report which gives a visual summary of the sequencing quality and identifies any issues, such as low-quality reads or adapter contamination.

Step 5: Trim Low-Quality Bases Using fastp

fastp \
  -i ~/1/SRR13018652_1.fastq \
  -I ~/1/SRR13018652_2.fastq \
  -o ~/1/SRR13018652_1.trim.fastq \
  -O ~/1/SRR13018652_2.trim.fastq \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --cut_front \
  --cut_tail \
  --cut_mean_quality 20 \
  --length_required 30 \
  --html ~/1/fastp_report.html \
  --json ~/1/fastp_report.json

Explanation:
This command runs fastp, a tool for quality control and trimming. It:
Removes adapter sequences from the paired-end reads.
Trims low-quality bases from the beginning and end of the reads.
Ensures that reads have a minimum length of 30 bases and a mean quality score above 20.
The resulting trimmed FASTQ files are saved in the specified output files.

Step 6: QC After Trimming

fastqc ~/1/SRR13018652_1.trim.fastq ~/1/SRR13018652_2.trim.fastq -o ~/1/fastqc_reports

Explanation:
Runs FastQC again on the trimmed FASTQ files to ensure that the trimming process did not introduce any problems. This step provides insight into how effective the trimming process was.

Step 7: Align Reads to the Reference Genome Using BWA-MEM2

bwa-mem2 mem -t 16 \
  -R '@RG\tID:SRR13018652\tSM:SRR13018652\tPL:ILLUMINA' \
  ~/hg38.fa \
  ~/1/SRR13018652_1.trim.fastq \
  ~/1/SRR13018652_2.trim.fastq \

Explanation:
Aligns the trimmed reads to the hg38 reference genome using BWA-MEM2. The -t 16 flag uses 16 threads to speed up the alignment process, and the -R flag specifies read group information (important for downstream analysis).

Step 8: Sort BAM File by Coordinate
  
~/1/SRR13018652.sam

Explanation:
Aligns the trimmed reads to the hg38 reference genome using BWA-MEM2. The -t 16 flag uses 16 threads to speed up the alignment process, and the -R flag specifies read group information (important for downstream analysis).

Step 9: Mark Duplicates
gatk MarkDuplicates \
  -I ~/1/SRR13018652.sorted.bam \
  -O ~/1/SRR13018652.dedup.bam \
  -M ~/1/SRR13018652.metrics.txt
Explanation:
Marks duplicate reads in the aligned BAM file using GATK MarkDuplicates. Duplicate reads are typically artifacts of the sequencing process and should be removed to avoid false variant calls.

Step 10: Index BAM File
samtools index ~/1/SRR13018652.dedup.bam
Explanation:
Indexes the deduplicated BAM file. BAM files need to be indexed for efficient access during variant calling.

Step 11: BaseRecalibrator (Base Quality Score Recalibration)
gatk BaseRecalibrator \
  -I ~/1/SRR13018652.dedup.bam \
  -R hg38.fa \
  --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O SRR13018652.recal.table
Explanation:
This step uses GATK BaseRecalibrator to perform base quality score recalibration (BQSR). The purpose of BQSR is to adjust the base quality scores in the sequencing data to account for systematic errors caused by the sequencing technology. This improves the accuracy of variant calling later in the pipeline.
--known-sites: This option provides known variant sites (such as SNPs and indels) from databases like dbSNP and Mills, which are used to recalibrate the base quality scores.
-O SRR13018652.recal.table: This produces an output recalibration table (SRR13018652.recal.table) that will be used in the next step to adjust the base quality scores in the BAM file.

Step 12: Apply Base Quality Score Recalibration (BQSR)
gatk ApplyBQSR \
  -R hg38.fa \
  -I ~/1/SRR13018652.dedup.bam \
  --bqsr-recal-file SRR13018652.recal.table \
  -O SRR13018652.recal.bam
Explanation:
In this step, GATK ApplyBQSR applies the recalibrated base quality scores to the deduplicated BAM file, which was produced in Step 9. This helps improve variant calling accuracy by correcting systematic biases in the base quality scores.
--bqsr-recal-file SRR13018652.recal.table: This option applies the recalibration table (SRR13018652.recal.table) generated in the previous step.
-O SRR13018652.recal.bam: This produces a new BAM file (SRR13018652.recal.bam) with the adjusted base quality scores, ready for variant calling.

