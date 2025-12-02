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

#How to Run the Pipeline

#Prerequisites
Before running the pipeline, ensure you have the following installed:
- `GATK`
- `BWA-MEM2`
- `bcftools`
- `VEP`
- `samtools`
- `fastp`
- `fastqc`
- `python3` (optional, if you need to use custom scripts)

#Step 1: Download data from SRA
prefetch SRP291993

#Step 2: Convert SRA to FASTQ Format
fastq-dump --split-files --gzip ~/1/SRR13018652.sra
Explanation:
Converts the SRA file (from Step 1) into FASTQ format. The --split-files flag ensures that paired-end data is separated into two files. The --gzip flag compresses the output files to save space.

Step 3: Decompress the FASTQ Files
