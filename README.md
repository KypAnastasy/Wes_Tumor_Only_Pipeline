# Project: WES Tumor-Only Analysis Pipeline (Breast Cancer Case Study)

## Language Selection / Выбор языка

[English](#readme-english) | [Русский](#readme-русский)

---

## **README (English)**

### Project Overview
Tumor-Only WES Analysis Pipeline for Breast Cancer

Description
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

>prefetch SRP291993

Step 2: Convert SRA to FASTQ Format

>fastq-dump --split-files --gzip ~/1/SRR13018652.sra

Explanation:
Converts the SRA file (from Step 1) into FASTQ format. The --split-files flag ensures that paired-end data is separated into two files. The --gzip flag compresses the output files to save space.

Step 3: Decompress the FASTQ Files

>gunzip ~/1/SRR13018652_1.fastq.gz

>gunzip ~/1/SRR13018652_2.fastq.gz

Explanation:
Decompress the gzipped FASTQ files so they can be used in the next steps. The two FASTQ files correspond to the paired-end reads.

Step 4: Quality Control (QC) Before Trimming

>fastqc ~/1/SRR13018652_1.fastq ~/1/SRR13018652_2.fastq -o ~/1/fastqc_reports

Explanation:
Runs FastQC on the raw FASTQ files to assess their quality. This generates an HTML report which gives a visual summary of the sequencing quality and identifies any issues, such as low-quality reads or adapter contamination.

Step 5: Trim Low-Quality Bases Using fastp

>fastp \
 >-i ~/1/SRR13018652_1.fastq \
  >-I ~/1/SRR13018652_2.fastq \
  >-o ~/1/SRR13018652_1.trim.fastq \
  >-O ~/1/SRR13018652_2.trim.fastq \
  >--detect_adapter_for_pe \
  >--trim_poly_g \
  >--cut_front \
  >--cut_tail \
  >--cut_mean_quality 20 \
  >--length_required 30 \
  >--html ~/1/fastp_report.html \
  >--json ~/1/fastp_report.json

Explanation:
This command runs fastp, a tool for quality control and trimming. It:
Removes adapter sequences from the paired-end reads.
Trims low-quality bases from the beginning and end of the reads.
Ensures that reads have a minimum length of 30 bases and a mean quality score above 20.
The resulting trimmed FASTQ files are saved in the specified output files.

Step 6: QC After Trimming

>fastqc ~/1/SRR13018652_1.trim.fastq ~/1/SRR13018652_2.trim.fastq -o ~/1/fastqc_reports

Explanation:
Runs FastQC again on the trimmed FASTQ files to ensure that the trimming process did not introduce any problems. This step provides insight into how effective the trimming process was.

Step 7: Align Reads to the Reference Genome Using BWA-MEM2

>bwa-mem2 mem -t 16 \
  >-R '@RG\tID:SRR13018652\tSM:SRR13018652\tPL:ILLUMINA' \
  >~/hg38.fa \
  >~/1/SRR13018652_1.trim.fastq \
  >~/1/SRR13018652_2.trim.fastq \

Explanation:
Aligns the trimmed reads to the hg38 reference genome using BWA-MEM2. The -t 16 flag uses 16 threads to speed up the alignment process, and the -R flag specifies read group information (important for downstream analysis).

Step 8: Sort BAM File by Coordinate
  
>~/1/SRR13018652.sam

Explanation:
Aligns the trimmed reads to the hg38 reference genome using BWA-MEM2. The -t 16 flag uses 16 threads to speed up the alignment process, and the -R flag specifies read group information (important for downstream analysis).

Step 9: Mark Duplicates
>gatk MarkDuplicates \
  >-I ~/1/SRR13018652.sorted.bam \
  >-O ~/1/SRR13018652.dedup.bam \
  >-M ~/1/SRR13018652.metrics.txt
Explanation:
Marks duplicate reads in the aligned BAM file using GATK MarkDuplicates. Duplicate reads are typically artifacts of the sequencing process and should be removed to avoid false variant calls.

Step 10: Index BAM File

>samtools index ~/1/SRR13018652.dedup.bam

Explanation:
Indexes the deduplicated BAM file. BAM files need to be indexed for efficient access during variant calling.

Step 11: BaseRecalibrator (Base Quality Score Recalibration)

>gatk BaseRecalibrator \
  >-I ~/1/SRR13018652.dedup.bam \
  >-R hg38.fa \
  >--known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  >--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  >-O SRR13018652.recal.table

Explanation:
This step uses GATK BaseRecalibrator to perform base quality score recalibration (BQSR). The purpose of BQSR is to adjust the base quality scores in the sequencing data to account for systematic errors caused by the sequencing technology. This improves the accuracy of variant calling later in the pipeline.

--known-sites: This option provides known variant sites (such as SNPs and indels) from databases like dbSNP and Mills, which are used to recalibrate the base quality scores.

-O SRR13018652.recal.table: This produces an output recalibration table (SRR13018652.recal.table) that will be used in the next step to adjust the base quality scores in the BAM file.

Step 12: Apply Base Quality Score Recalibration (BQSR)

>gatk ApplyBQSR \
  >-R hg38.fa \
  >-I ~/1/SRR13018652.dedup.bam \
  >--bqsr-recal-file SRR13018652.recal.table \
  >-O SRR13018652.recal.bam
>
Explanation:
In this step, GATK ApplyBQSR applies the recalibrated base quality scores to the deduplicated BAM file, which was produced in Step 9. This helps improve variant calling accuracy by correcting systematic biases in the base quality scores.

--bqsr-recal-file SRR13018652.recal.table: This option applies the recalibration table (SRR13018652.recal.table) generated in the previous step.

-O SRR13018652.recal.bam: This produces a new BAM file (SRR13018652.recal.bam) with the adjusted base quality scores, ready for variant calling.

Step 13: Mutect2 (Somatic Variant Calling)

>gatk Mutect2 \
  >-R hg38.fa \
  >-I ~/1/SRR13018652.recal.bam \
  -tumor SRR13018652 \
  >--germline-resource af-only-gnomad.hg38.vcf.gz \
  >--f1r2-tar-gz ~/1/SRR13018652.f1r2.tar.gz \
  >-O ~/1/SRR13018652.unfiltered.vcf.gz

Explanation:

Here, GATK Mutect2 is used to call somatic variants in the tumor sample in tumor-only mode. The tool identifies mutations specific to the tumor, such as single nucleotide variants (SNVs) and indels.

--germline-resource af-only-gnomad.hg38.vcf.gz: This provides the germline variant resource to differentiate somatic mutations from normal variants.

--f1r2-tar-gz ~/1/SRR13018652.f1r2.tar.gz: This option provides the F1R2 files, which are used for orientation bias correction.

The output file (SRR13018652.unfiltered.vcf.gz) contains all the variants identified by Mutect2, but these have not yet been filtered.

Step 14: Learn Read Orientation Model

>gatk LearnReadOrientationModel \
  >-I ~/1/SRR13018652.f1r2.tar.gz \
  >-O ~/1/SRR13018652.read-orientation-model.tar.gz

Explanation:

This step trains a read orientation model using the F1R2 files. These files are essential for detecting and correcting read orientation biases that can occur during sequencing, which could lead to incorrect variant calls.

-O ~/1/SRR13018652.read-orientation-model.tar.gz: The model is saved to the specified output file (read-orientation-model.tar.gz), which will be used later for filtering variants in the next steps.

Step 15: GetPileupSummaries (Contamination Estimation)
>gatk GetPileupSummaries \
  >-I ~/1/SRR13018652.recal.bam \
  >-V af-only-gnomad.hg38.vcf.gz \
  >-L af-only-gnomad.hg38.vcf.gz \
  >-O ~/1/SRR13018652.pileups.table

Explanation:

This step calculates pileup summaries using GATK GetPileupSummaries, which estimates the contamination in the tumor sample. Tumor samples can be contaminated with normal cells, which could affect variant calling. This step helps assess the level of normal cell contamination in the sample.

-V af-only-gnomad.hg38.vcf.gz: The known variant sites from the gnomAD database are used to help estimate contamination levels.

-O ~/1/SRR13018652.pileups.table: The output file (pileups.table) contains the contamination estimates that will be used for variant filtering in the next steps.

Step 16: Calculate Contamination

>gatk CalculateContamination \
  >-I ~/1/SRR13018652.pileups.table \
  >-O ~/1/SRR13018652.contamination.table \
  >--tumor-segmentation ~/1/SRR13018652.segments.table

Explanation:

In this step, GATK CalculateContamination computes the contamination level of the tumor sample using the pileup summary table generated in the previous step. The contamination table is essential for filtering out variants that may be due to contamination with normal cells.

-O ~/1/SRR13018652.contamination.table: This generates the contamination table (contamination.table), which will be used in the next step for filtering variants.

--tumor-segmentation ~/1/SRR13018652.segments.table: This provides the segmentation information of the tumor, which is necessary for accurate contamination estimation.

Step 17: FilterMutectCalls (Variant Filtering)

>gatk FilterMutectCalls \
  >-V ~/1/SRR13018652.unfiltered.vcf.gz \
  >-R hg38.fa \
  >--contamination-table ~/1/SRR13018652.contamination.table \
  >--tumor-segmentation ~/1/SRR13018652.segments.table \
  >--ob-priors ~/1/SRR13018652.read-orientation-model.tar.gz \
 > -O ~/1/SRR13018652.filtered.vcf.gz

Explanation:

This step uses GATK FilterMutectCalls to filter the variants that were called in Step 13 (SRR13018652.unfiltered.vcf.gz). The filtering is based on several criteria, including contamination levels, read orientation biases, and tumor segmentation.

--contamination-table ~/1/SRR13018652.contamination.table: The contamination information is applied here to correct for possible contamination from normal cells.

--ob-priors ~/1/SRR13018652.read-orientation-model.tar.gz: The orientation bias model is applied for further filtering of variants.

The output file (SRR13018652.filtered.vcf.gz) contains the filtered variants.

Step 18: Hard Filtering (Apply Custom Filters)

>bcftools filter -i \
  >'FORMAT/DP>=20 && FORMAT/AD[0:1]>=5 && (FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]))>=0.05 && FILTER="PASS"' \
  >~/1/SRR13018652.filtered.vcf.gz \
  >-Oz -o ~/1/SRR13018652.final.vcf.gz

Explanation:

This step uses bcftools to apply hard filters on the variants. The filters are based on criteria such as:

Depth of coverage (DP),

Allelic depth (AD),

Allele frequency (AF).

The result is a highly filtered VCF file (SRR13018652.final.vcf.gz) with variants that meet all the specified criteria.

Step 19: Decompress Final VCF File

>gunzip -c /home/ser/1/SRR13018652.final.vcf.gz > /home/ser/1/SRR13018652.final.vcf

Explanation:

This step decompresses the final VCF file, making it available as a standard text file (SRR13018652.final.vcf) for further analysis or visualization.

Step 20: VEP Annotation (Variant Effect Predictor)

>vep -i ~/1/SRR13018652.final.vcf \
    >-o ~/1/SRR13018652.annotated.vcf \
    >--cache \
    >--dir_cache /home/ser/.vep \
    >--species homo_sapiens \
    >--assembly GRCh38 \
    >--symbol \
    >--canonical \
    >--protein \
    >--af \
    >--sift b \
    >--polyphen b \
    >--mane \
    >--vcf \
    >--force_overwrite


Explanation:

In this step, VEP (Variant Effect Predictor) is used to annotate the final VCF file with detailed information about the functional consequences of the variants.

--cache: Uses a local cache to speed up the annotation process.

--species homo_sapiens: Specifies that the annotations are for human variants.

--assembly GRCh38: Specifies the reference genome version (GRCh38).

--symbol: Includes gene symbols in the annotation output.

--canonical: Selects the canonical transcript for each gene.

--protein: Includes protein information like protein domains and functional effects.

--af: Adds allele frequency data to the annotations.

--sift b and --polyphen b: These options use SIFT and PolyPhen to predict the impact of the variants on protein function (b = both prediction scores).

--mane: Includes MANE (Mutant Allele-specific Nucleotide Environment) annotations for transcript accuracy.

--vcf: Outputs the results in VCF format.

--force_overwrite: Overwrites the output file if it already exists.

The output (SRR13018652.annotated.vcf) contains the annotated variants, providing detailed information on their functional impact.

Step 21: Create Breast Cancer Genes List

>cat > /home/1/breast_cancer_genes.txt << 'EOF'
>BRCA1
>BRCA2
>TP53
>PTEN
>PIK3CA
>AKT1
>ESR1
>ERBB2
>CDH1
>STK11
>PALB2
>CHEK2
>ATM
>BRIP1
>EOF

Explanation:

This step creates a text file (breast_cancer_genes.txt) containing a list of known breast cancer-related genes. These genes are selected because mutations in them are known to be associated with breast cancer.

The list includes key tumor suppressor genes (e.g., BRCA1, BRCA2) and oncogenes (e.g., PIK3CA, AKT1).

This list will later be used to filter variants to focus specifically on mutations in breast cancer-related genes.

Step 22: Filter Variants in Breast Cancer Genes

>bcftools view -h /home/1/SRR13018652.annotated.vcf > /home/ser/1/header.vcf

>bcftools view -H /home/ser/1/SRR13018652.annotated.vcf | \
  >grep -E "$(paste -s -d '|' /home/ser/1/breast_cancer_genes.txt)" | \
  >cat /home/ser/1/header.vcf - > /home/ser/1/SRR13018652.breast_cancer_genes.vcf

Explanation:

This step filters the annotated variants to retain only those in the genes listed in breast_cancer_genes.txt.

The first command (bcftools view -h) extracts the header from the VCF file and saves it to header.vcf.

The second command filters the variants using grep and the gene list (breast_cancer_genes.txt), and then merges the header with the filtered variants.

The final result is stored in SRR13018652.breast_cancer_genes.vcf, which contains only variants in breast cancer-related genes.

Step 23: Analyze Results (Variant Count and Distribution by Gene)

>echo "Total variants in breast cancer genes:"
>grep -v "^#" /home/ser/1/SRR13018652.breast_cancer_genes.vcf | wc -l

>echo "Distribution by genes:"
>grep -v "^#" /home/ser/1/SRR13018652.breast_cancer_genes.vcf | \
  >grep -o "CSQ=[^;]*" | \
  >tr ',' '\n' | \
  >cut -d'|' -f4 | \
  >sort | uniq -c | sort -nr

Explanation:

This step provides an analysis of the variants in the breast cancer-related genes:

Total Variants: It counts the total number of variants in the filtered VCF (SRR13018652.breast_cancer_genes.vcf), excluding comment lines (lines starting with #).

Distribution by Gene: It analyzes the distribution of the variants across the genes:

grep -o "CSQ=[^;]*" extracts the Consequence Annotation (CSQ) field from the VCF.

cut -d'|' -f4 extracts the gene symbol from the CSQ field.

sort | uniq -c | sort -nr counts the occurrences of each gene and sorts them in descending order.

Step 24: Extract Key Mutation Details

>echo "Detailed variants information:"
>grep -v "^#" /home/ser/1/SRR13018652.breast_cancer_genes.vcf | \
  >awk -F'\t' '{
    >printf "%-10s %-10s %-15s", $1, $2, $4 ">" $5;
    >if ($8 ~ /CSQ=/) {
      >split($8, info, ";");
      >for (i in info) {
        >if (info[i] ~ /^CSQ=/) {
          >split(info[i], csq, "|");
          >printf " %s (%s)", csq[4], csq[2];
        >}
      >}
    >}
    >print "";
  >}'


Explanation:

This step extracts detailed information about each key mutation:

It processes each variant line from SRR13018652.breast_cancer_genes.vcf, extracting details such as:

Chromosome ($1),

Position ($2),

Reference and alternate alleles ($4 and $5).

It also includes the Consequence Annotation (CSQ), which provides information on the gene symbol, mutation type, and the functional effect (e.g., missense, frameshift).

The output displays these details in a readable format.

Step 25: Determine Patient Sex

>echo "Determining patient sex:"
>y_reads=$(samtools idxstats /home/ser/1/SRR13018652.recal.bam | awk '$1 == "chrY" {print $3}')
>echo "Y chromosome reads: $y_reads"

>if [ "$y_reads" -gt 1000 ]; then
    >echo "PATIENT SEX: MALE"
    >sex="MALE"
>else
    >echo "PATIENT SEX: FEMALE"
    >sex="FEMALE"
>fi

Explanation:

This step determines the biological sex of the patient based on the presence of Y chromosome reads in the BAM file:

samtools idxstats provides statistics about the reference chromosomes in the BAM file.

The script checks the number of reads aligned to the Y chromosome (chrY):

If the number of Y chromosome reads is greater than 1000, the patient is classified as male.

Otherwise, the patient is classified as female.

The determined sex (MALE or FEMALE) is stored in the sex variable.


