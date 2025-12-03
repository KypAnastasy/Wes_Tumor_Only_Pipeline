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

>gunzip -c /home/1/SRR13018652.final.vcf.gz > /home/1/SRR13018652.final.vcf

Explanation:

This step decompresses the final VCF file, making it available as a standard text file (SRR13018652.final.vcf) for further analysis or visualization.

Step 20: VEP Annotation (Variant Effect Predictor)

vep -i ~/1/SRR13018652.final.vcf \
    -o ~/1/SRR13018652.annotated.vcf \
    --cache \
    --dir_cache /home/.vep \
    --species homo_sapiens \
    --assembly GRCh38 \
    --symbol \
    --canonical \
    --protein \
    --af \
    --sift b \
    --polyphen b \
    --mane \
    --vcf \
    --force_overwrite


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

>bcftools view -h /home/1/SRR13018652.annotated.vcf > /home/1/header.vcf

>bcftools view -H /home/1/SRR13018652.annotated.vcf | \
  >grep -E "$(paste -s -d '|' /home/1/breast_cancer_genes.txt)" | \
  >cat /home/1/header.vcf - > /home/1/SRR13018652.breast_cancer_genes.vcf

Explanation:

This step filters the annotated variants to retain only those in the genes listed in breast_cancer_genes.txt.

The first command (bcftools view -h) extracts the header from the VCF file and saves it to header.vcf.

The second command filters the variants using grep and the gene list (breast_cancer_genes.txt), and then merges the header with the filtered variants.

The final result is stored in SRR13018652.breast_cancer_genes.vcf, which contains only variants in breast cancer-related genes.

Step 23: Analyze Results (Variant Count and Distribution by Gene)

>echo "Total variants in breast cancer genes:"
>grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | wc -l

>echo "Distribution by genes:"
>grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | \
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
>grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | \
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
>y_reads=$(samtools idxstats /home/1/SRR13018652.recal.bam | awk '$1 == "chrY" {print $3}')
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

Step 26: Detailed Analysis of Key Mutations

>cat > /home/1/analyze_key_mutations.sh << 'EOF'
>#!/bin/bash

>echo "Detailed analysis of key mutations"

>mutations=("chr17 7676154" "chr17 43063913" "chr3 179218303")

>for mutation in "${mutations[@]}"; do
    IFS=" " read -r chrom pos <<< "$mutation"
    
    echo ""
    echo "Analyzing mutation: $chrom $pos"
    echo "------------------------"
    
    variant_line=$(awk -v chrom="$chrom" -v pos="$pos" '$1 == chrom && $2 == pos' /home/ser/1/SRR13018652.breast_cancer_genes.vcf)
    
    if [ -n "$variant_line" ]; then
        chrom=$(echo "$variant_line" | awk '{print $1}')
        pos=$(echo "$variant_line" | awk '{print $2}')
        ref=$(echo "$variant_line" | awk '{print $4}')
        alt=$(echo "$variant_line" | awk '{print $5}')
        qual=$(echo "$variant_line" | awk '{print $6}')
        filter=$(echo "$variant_line" | awk '{print $7}')
        info=$(echo "$variant_line" | awk '{print $8}')
        sample=$(echo "$variant_line" | awk '{print $10}')
        
        csq_info=$(echo "$variant_line" | grep -o "CSQ=[^;]*")
        first_csq=$(echo "$csq_info" | cut -d',' -f1)
        
        symbol=$(echo "$first_csq" | cut -d'|' -f4)
        consequence=$(echo "$first_csq" | cut -d'|' -f2)
        impact=$(echo "$first_csq" | cut -d'|' -f3)
        protein_change=$(echo "$first_csq" | cut -d'|' -f11)
        amino_acids=$(echo "$first_csq" | cut -d'|' -f12)
        
        echo "Gene: $symbol"
        echo "Change: $ref -> $alt" 
        echo "Type: $consequence"
        echo "Impact: $impact"
        if [ -n "$protein_change" ] && [ "$protein_change" != "" ]; then
            echo "Protein change: $protein_change"
        fi
        if [ -n "$amino_acids" ] && [ "$amino_acids" != "" ]; then
            echo "Amino acids: $amino_acids"
        fi
        echo "Filter: $filter"
        echo "Quality: $qual"
        
        dp=$(echo "$info" | grep -o "DP=[0-9]*" | cut -d'=' -f2)
        af=$(echo "$info" | grep -o "AF=[0-9.]*" | cut -d'=' -f2)
        tlod=$(echo "$info" | grep -o "TLOD=[0-9.]*" | cut -d'=' -f2)
        
        echo "Depth (INFO): $dp"
        echo "Allele frequency (INFO): $af"
        if [ -n "$tlod" ]; then
            echo "Tumor LOD: $tlod"
        fi
        
        gt=$(echo "$sample" | cut -d':' -f1)
        ad=$(echo "$sample" | cut -d':' -f2)
        sample_af=$(echo "$sample" | cut -d':' -f3)
        sample_dp=$(echo "$sample" | cut -d':' -f4)
        
        echo "Genotype: $gt"
        echo "Allelic Depths (AD): $ad"
        echo "Allele frequency (FORMAT): $sample_af"
        echo "Depth (FORMAT): $sample_dp"
        
        ref_count=$(echo "$ad" | cut -d',' -f1)
        alt_count=$(echo "$ad" | cut -d',' -f2)
        if [ -n "$ref_count" ] && [ -n "$alt_count" ] && [ "$ref_count" -gt 0 ]; then
            vaf=$(echo "scale=4; $alt_count / ($ref_count + $alt_count)" | bc)
            echo "VAF (calculated): $vaf"
        fi
        
        echo ""
        echo "Additional information:"
        echo "$info" | tr ';' '\n' | grep -E "^(DP|AF|TLOD|MBQ|MMQ|MPOS)="
        
    else
        echo "Mutation $chrom $pos not found"
    fi
    echo ""
>done

>EOF

>chmod +x /home/1/analyze_key_mutations.sh
>/home/1/analyze_key_mutations.sh > /home/1/key_mutations_detailed_analysis.txt


Explanation:

This script performs a detailed analysis of specific key mutations in the breast cancer-related genes list. It examines specific mutations (given by chromosome and position) and extracts detailed variant information such as:

Gene symbol, type, and impact of mutation.

Genotype and depth of sequencing (AD, DP, AF).

Variant allele frequency (VAF).

Additional details like tumor log odds (TLOD), depth (DP), and allele frequency (AF) from the INFO field.

After running the script, the output is saved in key_mutations_detailed_analysis.txt.

Step 27: Create Final Clinical Report

>cat > /home/ser/1/final_clinical_report.md << EOF
># GENOMIC ANALYSIS REPORT
>## Whole Exome Sequencing - Tumor Only Analysis

>### PATIENT INFORMATION
>- Sample ID: SRR13018652
>- Biological Sex: $sex
>- Analysis Date: $(date +"%Y-%m-%d")

>### EXECUTIVE SUMMARY
>Whole exome sequencing identified three clinically significant mutations in key cancer genes.

>### KEY MUTATIONS IDENTIFIED

>#### 1. TP53 Mutation
>- Location: chr17:7676154
>- Variant: G>C (Missense)
>- VAF: 8.9%

>#### 2. BRCA1 Mutation
>- Location: chr17:43063913
>- Variant: G>C (Missense) 
>- VAF: 15.1%

>#### 3. PIK3CA Mutation
>- Location: chr3:179218303
>- Variant: G>A (Missense)
>- VAF: 36.9%

>### CLINICAL RECOMMENDATIONS

>#### Genetic Counseling & Testing
>1. Confirm germline status of BRCA1 mutation
>2. Family member testing recommended
>3. Discuss reproductive implications

>#### Cancer Surveillance
>1. Enhanced breast cancer screening
>2. $([ "$sex" = "MALE" ] && echo "Prostate cancer screening starting at age 40")
>3. Consider pancreatic cancer screening

>#### Therapeutic Considerations
>1. PARP inhibitors for BRCA1-related cancers
>2. PI3K inhibitors for PIK3CA-mutated cancers
>3. Clinical trial eligibility assessment

>### METHODS
>- Sequencing Technology: Whole Exome Sequencing
>- Reference Genome: hg38
>- Variant Caller: GATK Mutect2 (tumor-only mode)
>- Annotation: Ensembl VEP

>### LIMITATIONS
>- Tumor-only analysis without matched normal
- Germline vs somatic status undetermined
>>- Functional validation required for missense variants

>Analysis completed: $(date)
>EOF

Explanation:

This step creates a clinical report in markdown format. The report includes:

Patient Information: Sample ID and biological sex.

Executive Summary: A brief overview of the key mutations found in the analysis (e.g., TP53, BRCA1, PIK3CA).

Clinical Recommendations: Suggestions for genetic counseling, cancer surveillance, and possible therapeutic considerations based on the detected mutations.

Methods: Details on the sequencing technology, reference genome, variant caller, and annotation methods used in the analysis.

Limitations: Acknowledgement of the limitations of the analysis, such as the tumor-only approach without a matched normal sample.

The clinical report is saved as final_clinical_report.md.

Step 28: Prostate Cancer Analysis for Male Patients

>if [ "$sex" = "MALE" ]; then
    >cat > /home/1/prostate_cancer_genes.txt << 'EOF'
>BRCA1
>BRCA2
>ATM
>CHEK2
>TP53
>HOXB13
>MLH1
>MSH2
>MSH6
>PMS2
>EOF

    bcftools view -h /home/1/SRR13018652.annotated.vcf > /home/1/header_prostate.vcf
    
    bcftools view -H /home/1/SRR13018652.annotated.vcf | \
      grep -E "$(paste -s -d '|' /home/1/prostate_cancer_genes.txt)" | \
      cat /home1/header_prostate.vcf - > /home/1/SRR13018652.prostate_cancer_genes.vcf

    prostate_variants=$(grep -v "^#" /home/1/SRR13018652.prostate_cancer_genes.vcf | wc -l)
    echo "Variants in prostate cancer genes: $prostate_variants"
    
    echo "" >> /home/1/final_clinical_report.md
    echo "### PROSTATE CANCER RISK ASSESSMENT" >> /home/1/final_clinical_report.md
    echo "Based on BRCA1 and TP53 mutations:" >> /home/1/final_clinical_report.md
    echo "- Lifetime prostate cancer risk: 20-25%" >> /home/1/final_clinical_report.md
    echo "- Consider PSA screening starting at age 40" >> /home/1/final_clinical_report.md
>fi

Explanation:
This step checks if the patient is male and, if so, analyzes prostate cancer-related mutations (in genes like BRCA1, TP53, etc.). It adds information about prostate cancer risk and screening recommendations to the final clinical report.

Step 29: Create Results Summary Table

>cat > /home/1/results_summary.csv << EOF
>Category,Value
>Sample_ID,SRR13018652
>Patient_Sex,$sex
>Total_Variants,$(grep -v "^#" /home/1/SRR13018652.final.vcf | wc -l)
>Breast_Cancer_Gene_Variants,$(grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | wc -l)
Key_Driver_Mutations,3
>TP53_VAF,0.089
>BRCA1_VAF,0.151
>PIK3CA_VAF,0.369
>Analysis_Date,$(date +"%Y-%m-%d")

Explanation:

Create a CSV file summarizing the analysis results:

This step generates a summary table in CSV format, which includes key metrics such as:

Sample ID (SRR13018652).

Patient sex ($sex).

The total number of variants found in the final VCF (Total_Variants).

The number of variants found in breast cancer genes (Breast_Cancer_Gene_Variants).

The count of key driver mutations (hardcoded to 3 in this case).

Variant Allele Frequencies (VAF) for TP53, BRCA1, and PIK3CA mutations.

The date of analysis.

Dynamic Variables:

The table pulls dynamic data such as variant counts and analysis date by running commands within the $(...) syntax (e.g., counting the lines in the VCF files for total variants).

File Created:

A file results_summary.csv is created to provide a quick overview of the analysis results.

Step 30: Organize Final Results

>mkdir -p /home/1/final_results/01_vcf_files
>mkdir -p /home/1/final_results/02_reports
>mkdir -p /home/1/final_results/03_quality_metrics

>cp /home/1/SRR13018652.final.vcf /home/1/final_results/01_vcf_files/
>cp /home/1/SRR13018652.annotated.vcf /home/1/final_results/01_vcf_files/
>cp /home/1/SRR13018652.breast_cancer_genes.vcf /home/1/final_results/01_vcf_files/
>[ "$sex" = "MALE" ] && cp /home/1/SRR13018652.prostate_cancer_genes.vcf /home/1/final_results/01_vcf_files/

>cp /home/1/final_clinical_report.md /home/1/final_results/02_reports/
>cp /home/1/results_summary.csv /home/1/final_results/02_reports/
>cp /home/1/key_mutations_detailed_analysis.txt /home/1/final_results/02_reports/

>cp -r /home/1/fastqc_reports /home/1/final_results/03_quality_metrics/
>cp /home/1/fastp_report.html /home/1/final_results/03_quality_metrics/

>echo "ANALYSIS COMPLETED"
>echo "Results available in: /home/1/final_results/"
>echo "Clinical report: /home/1/final_results/02_reports/final_clinical_report.md"

Explanation:

Create Directories for Final Results:

mkdir -p /home/1/final_results/01_vcf_files creates a directory for VCF files.

mkdir -p /home/1/final_results/02_reports creates a directory for reports.

mkdir -p /home/1/final_results/03_quality_metrics creates a directory for quality metrics (QC).

Organize Files into Proper Directories:

VCF files: The final VCF files (final.vcf, annotated.vcf, breast_cancer_genes.vcf, and optionally, prostate_cancer_genes.vcf for male patients) are copied into the 01_vcf_files directory.

Reports: The clinical report (final_clinical_report.md), results summary (results_summary.csv), and detailed mutations analysis (key_mutations_detailed_analysis.txt) are copied into the 02_reports directory.

Quality Metrics: QC reports (fastqc_reports and fastp_report.html) are copied into the 03_quality_metrics directory.

Completion Message:

After organizing all the final results into the proper directories, the script prints a completion message with the paths to the results.

Final Results Path:

The results can be found in the /home/1/final_results/ directory, and the clinical report can be accessed at /home/1/final_results/02_reports/final_clinical_report.md.

## **README (Русский)**

### Обзор проекта
Тумор-Онли WES Анализ для Рака Молочной Железы

Описание
Этот репозиторий содержит полный пиплайн для анализа всего экзома (WES) в режиме "тумор-онли" (анализ только опухоли), специально разработанный для рака молочной железы. Пиплайн включает несколько этапов, начиная с загрузки необработанных данных и заканчивая генерацией клинического отчета.

В анализе используются такие инструменты, как GATK, BWA-MEM2, bcftools, VEP и другие, для обработки необработанных NGS-данных, вызова мутаций и предоставления клинически значимых результатов.

Как запустить пиплайн

Предварительные требования

Перед запуском пиплайна убедитесь, что у вас установлены следующие инструменты:
- `GATK`
- `BWA-MEM2`
- `bcftools`
- `VEP`
- `samtools`
- `fastp`
- `fastqc`
- `python3` (опционально, если нужны кастомные скрипты)

Шаг 1: Загрузка данных из SRA

>prefetch SRP291993

Шаг 2: Конвертация SRA в формат FASTQ

>fastq-dump --split-files --gzip ~/1/SRR13018652.sra

Объяснение:
Конвертирует файл SRA (из Шага 1) в формат FASTQ. Флаг --split-files гарантирует, что данные с парными ридами будут разделены на два файла. Флаг --gzip сжимает выходные файлы, чтобы сэкономить место.

Шаг 3: Объяснение:

Конвертирует файл SRA (из Шага 1) в формат FASTQ. Флаг --split-files гарантирует, что данные с парными ридами будут разделены на два файла. Флаг --gzip сжимает выходные файлы, чтобы сэкономить место.

Шаг 3: Распаковка файлов FASTQ

>gunzip ~/1/SRR13018652_1.fastq.gz

>gunzip ~/1/SRR13018652_2.fastq.gz

Объяснение:
Распаковывает сжатые файлы FASTQ, чтобы они могли быть использованы на следующих шагах. Два файла FASTQ соответствуют парным ридам.

Шаг 4: Контроль качества (QC) до тримминга (обрезки)

>fastqc ~/1/SRR13018652_1.fastq ~/1/SRR13018652_2.fastq -o ~/1/fastqc_reports

Объяснение:
Запускает FastQC на исходных файлах FASTQ для оценки их качества. Это генерирует HTML-отчет, который дает визуальное представление о качестве секвенирования и выявляет возможные проблемы, такие как низкокачественные риды или загрязнение адаптерами.

Шаг 5: Обрезка низкокачественных оснований с помощью fastp

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

Объяснение:
Эта команда запускает fastp — инструмент для контроля качества и обрезки. Он:

Удаляет последовательности адаптеров из парных ридов.

Обрезает низкокачественные основания с начала и конца ридов.

Обеспечивает, чтобы риды имели минимальную длину 30 оснований и средний балл качества выше 20.

Результирующие обрезанные файлы FASTQ сохраняются в указанные выходные файлы.

Шаг 6: Контроль качества (QC) после обрезки

>fastqc ~/1/SRR13018652_1.trim.fastq ~/1/SRR13018652_2.trim.fastq -o ~/1/fastqc_reports

Объяснение:
Снова запускает FastQC на обрезанных файлах FASTQ, чтобы убедиться, что процесс обрезки не привел к появлению проблем. Этот шаг дает представление о том, насколько эффективно была выполнена обрезка.

Шаг 7: Выравнивание ридов с референсным геномом с помощью BWA-MEM2

>bwa-mem2 mem -t 16 \
  >-R '@RG\tID:SRR13018652\tSM:SRR13018652\tPL:ILLUMINA' \
  >~/hg38.fa \
  >~/1/SRR13018652_1.trim.fastq \
  >~/1/SRR13018652_2.trim.fastq \

Объяснение:
Выравнивает обрезанные риды с референсным геномом hg38 с использованием BWA-MEM2. Флаг -t 16 использует 16 потоков для ускорения процесса выравнивания, а флаг -R указывает информацию о группе ридов (что важно для дальнейшего анализа)

Шаг 8: Сортировка BAM-файла по координатам
  
>~/1/SRR13018652.sam

Объяснение:
Сортирует BAM-файл по координатам, чтобы риды были упорядочены в соответствии с их положением на референсном геноме. Это важно для последующих шагов анализа, таких как выявление дубликатов и выравнивание.

Шаг 9: Отметка дубликатов

>gatk MarkDuplicates \
  >-I ~/1/SRR13018652.sorted.bam \
  >-O ~/1/SRR13018652.dedup.bam \
  >-M ~/1/SRR13018652.metrics.txt

Объяснение:
Отмечает дубликаты ридов в выровненном BAM-файле с помощью GATK MarkDuplicates. Дубликаты ридов обычно являются артефактами процесса секвенирования, и их следует удалить, чтобы избежать ложных вызовов вариантов.

Шаг 10: Индексация BAM-файла

>samtools index ~/1/SRR13018652.dedup.bam

Объяснение:
Индексирует BAM-файл после удаления дубликатов. BAM-файлы должны быть проиндексированы для эффективного доступа при вызове вариантов.

Шаг 11: BaseRecalibrator (Перекалибровка баллов качества оснований)

>gatk BaseRecalibrator \
  >-I ~/1/SRR13018652.dedup.bam \
  >-R hg38.fa \
  >--known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  >--known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  >-O SRR13018652.recal.table

Объяснение:
Этот шаг использует GATK BaseRecalibrator для выполнения перекалибровки баллов качества оснований (BQSR). Цель BQSR — откорректировать баллы качества оснований в данных секвенирования, чтобы учесть систематические ошибки, вызванные технологией секвенирования. Это улучшает точность вызова вариантов на следующих этапах пайплайна.

--known-sites: Эта опция предоставляет известные участки вариантов (такие как SNP и инделы) из баз данных, таких как dbSNP и Mills, которые используются для перекалибровки баллов качества оснований.

-O SRR13018652.recal.table: Этот шаг генерирует выходной файл таблицы перекалибровки (SRR13018652.recal.table), который будет использован на следующем этапе для корректировки баллов качества оснований в BAM-файле.

Шаг 12: Применение повторной калибровки оценок качества оснований (BQSR)

>gatk ApplyBQSR \
  >-R hg38.fa \
  >-I ~/1/SRR13018652.dedup.bam \
  >--bqsr-recal-file SRR13018652.recal.table \
  >-O SRR13018652.recal.bam
>
Объяснение:
На этом шаге инструмент GATK ApplyBQSR применяет повторно откалиброванные оценки качества оснований к файлу BAM с удалёнными дубликатами, который был получен на Шаге 9. Это помогает повысить точность выявления вариантов (вариант-коллинга) за счёт корректировки систематических смещений в оценках качества оснований.

--bqsr-recal-file SRR13018652.recal.table: This option applies the recalibration table (SRR13018652.recal.table) generated in the previous step.

-O SRR13018652.recal.bam: This produces a new BAM file (SRR13018652.recal.bam) with the adjusted base quality scores, ready for variant calling.

Шаг 13: Mutect2 (вызов соматического варианта)

>gatk Mutect2 \
  >-R hg38.fa \
  >-I ~/1/SRR13018652.recal.bam \
  -tumor SRR13018652 \
  >--germline-resource af-only-gnomad.hg38.vcf.gz \
  >--f1r2-tar-gz ~/1/SRR13018652.f1r2.tar.gz \
  >-O ~/1/SRR13018652.unfiltered.vcf.gz

Объяснение:
В данном случае GATK Mutect2 используется для выявления соматических вариантов в образце опухоли в режиме «только опухоль». Инструмент выявляет мутации, специфичные для опухоли, такие как однонуклеотидные варианты (SNV) и инделы.

--germline-resource af-only-gnomad.hg38.vcf.gz: Этот шаг предоставляет ресурс герминальных вариантов для дифференциации соматических мутаций от нормальных вариантов.

--f1r2-tar-gz ~/1/SRR13018652.f1r2.tar.gz: Эта опция предоставляет файлы F1R2, которые используются для коррекции смещения по ориентации.

Выходной файл (SRR13018652.unfiltered.vcf.gz) содержит все варианты, выявленные Mutect2, но они ещё не были отфильтрованы.

Шаг 14: Построение модели ориентации ридов (чтений)

>gatk LearnReadOrientationModel \
  >-I ~/1/SRR13018652.f1r2.tar.gz \
  >-O ~/1/SRR13018652.read-orientation-model.tar.gz

Объяснение:

На этом шаге обучается модель ориентации ридов с использованием файлов F1R2. Эти файлы необходимы для выявления и коррекции смещений по ориентации прочтений, которые могут возникнуть при секвенировании и привести к ошибочному определению вариантов.

-O ~/1/SRR13018652.read-orientation-model.tar.gz: Модель сохраняется в указанный выходной файл (read-orientation-model.tar.gz), который будет использоваться позже для фильтрации вариантов на следующих шагах.

Шаг 15: Получение сводных данных о нуклеотидных частотах (Оценка контаминации)

>gatk GetPileupSummaries \
  >-I ~/1/SRR13018652.recal.bam \
  >-V af-only-gnomad.hg38.vcf.gz \
  >-L af-only-gnomad.hg38.vcf.gz \
  >-O ~/1/SRR13018652.pileups.table

Объяснение:

Этот шаг вычисляет сводки о нуклеотидных частотах с помощью инструмента GATK GetPileupSummaries для оценки уровня контаминации в опухолевом образце. Опухолевые образцы могут быть загрязнены нормальными клетками, что может повлиять на выявление вариантов. Этот шаг помогает оценить уровень такого загрязнения.

-V af-only-gnomad.hg38.vcf.gz: Используются известные сайты вариантов из базы данных gnomAD для помощи в оценке уровня контаминации.

-O ~/1/SRR13018652.pileups.table: Выходной файл (pileups.table) содержит оценки контаминации, которые будут использоваться для фильтрации вариантов на следующих шагах.

Шаг 16: Расчет уровня контаминации

>gatk CalculateContamination \
  >-I ~/1/SRR13018652.pileups.table \
  >-O ~/1/SRR13018652.contamination.table \
  >--tumor-segmentation ~/1/SRR13018652.segments.table

Объяснение:

На этом шаге инструмент GATK CalculateContamination вычисляет уровень контаминации опухолевого образца, используя таблицу сводок о нуклеотидных частотах, созданную на предыдущем шаге. Таблица контаминации необходима для фильтрации вариантов, которые могут быть вызваны загрязнением нормальными клетками.

-O ~/1/SRR13018652.contamination.table: Создаёт таблицу контаминации (contamination.table), которая будет использоваться на следующем шаге для фильтрации вариантов.

--tumor-segmentation ~/1/SRR13018652.segments.table: Предоставляет информацию о сегментации опухоли, необходимую для точной оценки контаминации.

Шаг 17: Фильтрация вызовов Mutect (Фильтрация вариантов)

>gatk FilterMutectCalls \
  >-V ~/1/SRR13018652.unfiltered.vcf.gz \
  >-R hg38.fa \
  >--contamination-table ~/1/SRR13018652.contamination.table \
  >--tumor-segmentation ~/1/SRR13018652.segments.table \
  >--ob-priors ~/1/SRR13018652.read-orientation-model.tar.gz \
 > -O ~/1/SRR13018652.filtered.vcf.gz

Объяснение:

Этот шаг использует инструмент GATK FilterMutectCalls для фильтрации вариантов, выявленных на Шаге 13 (файл SRR13018652.unfiltered.vcf.gz). Фильтрация основана на нескольких критериях, включая уровни контаминации, смещения по ориентации ридов и сегментацию опухоли.

--contamination-table ~/1/SRR13018652.contamination.table: Информация о контаминации применяется здесь для коррекции возможного загрязнения нормальными клетками.

--ob-priors ~/1/SRR13018652.read-orientation-model.tar.gz: Модель смещения по ориентации применяется для дальнейшей фильтрации вариантов.

Выходной файл (SRR13018652.filtered.vcf.gz) содержит отфильтрованные варианты.

Шаг 18: Строгая фильтрация (Применение пользовательских фильтров)

>bcftools filter -i \
  >'FORMAT/DP>=20 && FORMAT/AD[0:1]>=5 && (FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]))>=0.05 && FILTER="PASS"' \
  >~/1/SRR13018652.filtered.vcf.gz \
  >-Oz -o ~/1/SRR13018652.final.vcf.gz

Объяснение:

На этом шаге используется инструмент bcftools для применения строгих фильтров к вариантам. Фильтры основаны на таких критериях, как:

Глубина покрытия (DP),

Аллельная глубина (AD),

Частота аллели (AF).

Результатом является строго отфильтрованный файл VCF (SRR13018652.final.vcf.gz), содержащий варианты, соответствующие всем заданным критериям.

Шаг 19: Распаковка финального VCF файла

>gunzip -c /home/1/SRR13018652.final.vcf.gz > /home/1/SRR13018652.final.vcf

Объяснение:

Этот шаг распаковывает финальный VCF файл, делая его доступным в виде стандартного текстового файла (SRR13018652.final.vcf) для дальнейшего анализа или визуализации.

Шаг 20: Аннотация с помощью VEP (Variant Effect Predictor / Предиктор эффекта вариантов)

vep -i ~/1/SRR13018652.final.vcf \
    -o ~/1/SRR13018652.annotated.vcf \
    --cache \
    --dir_cache /home/.vep \
    --species homo_sapiens \
    --assembly GRCh38 \
    --symbol \
    --canonical \
    --protein \
    --af \
    --sift b \
    --polyphen b \
    --mane \
    --vcf \
    --force_overwrite


Объяснение:

На этом шаге используется VEP (Variant Effect Predictor / Предиктор эффекта вариантов) для аннотирования финального VCF файла подробной информацией о функциональных последствиях вариантов.

--cache: Использует локальный кэш для ускорения процесса аннотации.

--species homo_sapiens: Указывает, что аннотации предназначены для вариантов человека.

--assembly GRCh38: Указывает версию референсного генома (GRCh38).

--symbol: Включает символы генов в выходные данные аннотации.

--canonical: Выбирает канонический транскрипт для каждого гена.

--protein: Включает информацию о белке, такую как белковые домены и функциональные эффекты.

--af: Добавляет данные о частоте аллелей к аннотациям.

--sift b и --polyphen b: Эти параметры используют SIFT и PolyPhen для прогнозирования влияния вариантов на функцию белка (b = оба показателя предсказания).

--mane: Включает аннотации MANE (Mutant Allele-specific Nucleotide Environment) для точности транскриптов.

--vcf: Выводит результаты в формате VCF.

--force_overwrite: Перезаписывает выходной файл, если он уже существует.

Выходной файл (SRR13018652.annotated.vcf) содержит аннотированные варианты с подробной информацией об их функциональном воздействии.

Шаг 21: Создание списка генов рака молочной железы

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

Объяснение:

На этом шаге создаётся текстовый файл (breast_cancer_genes.txt), содержащий список известных генов, связанных с раком молочной железы. Эти гены выбраны потому, что мутации в них известны своей ассоциацией с данным заболеванием.

Список включает ключевые гены-супрессоры опухолей (например, BRCA1, BRCA2) и онкогены (например, PIK3CA, AKT1).

Этот список позже будет использован для фильтрации вариантов, чтобы сфокусироваться именно на мутациях в генах, связанных с раком молочной железы.

Шаг 22: Фильтрация вариантов в генах рака молочной железы

>bcftools view -h /home/1/SRR13018652.annotated.vcf > /home/1/header.vcf

>bcftools view -H /home/1/SRR13018652.annotated.vcf | \
  >grep -E "$(paste -s -d '|' /home/1/breast_cancer_genes.txt)" | \
  >cat /home/1/header.vcf - > /home/1/SRR13018652.breast_cancer_genes.vcf

Объяснение:

На этом шаге фильтруются аннотированные варианты, чтобы оставить только те, что находятся в генах из списка breast_cancer_genes.txt.

Первая команда (bcftools view -h) извлекает заголовок из VCF файла и сохраняет его в header.vcf.

Вторая команда фильтрует варианты с помощью grep и списка генов (breast_cancer_genes.txt), а затем объединяет заголовок с отфильтрованными вариантами.

Итоговый результат сохраняется в файл SRR13018652.breast_cancer_genes.vcf, который содержит только варианты в генах, связанных с раком молочной железы.

Шаг 23: Анализ результатов (подсчёт вариантов и их распределение по генам)

>echo "Total variants in breast cancer genes:"
>grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | wc -l

>echo "Distribution by genes:"
>grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | \
  >grep -o "CSQ=[^;]*" | \
  >tr ',' '\n' | \
  >cut -d'|' -f4 | \
  >sort | uniq -c | sort -nr

Объяснение:

На этом шаге проводится анализ вариантов в генах, связанных с раком молочной железы:

Общее количество вариантов: Подсчитывается общее число вариантов в отфильтрованном VCF-файле (SRR13018652.breast_cancer_genes.vcf), исключая строки-комментарии (строки, начинающиеся с #).

Распределение по генам: Анализируется распределение вариантов по генам.

grep -o "CSQ=[^;]*" извлекает из VCF поля аннотации эффекта вариантов (CSQ).

cut -d'|' -f4 извлекает символ гена из поля CSQ.

sort | uniq -c | sort -nr подсчитывает вхождения каждого гена и сортирует их по убыванию.

Шаг 24: Извлечение ключевых деталей мутаций

>echo "Detailed variants information:"
>grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | \
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

Объяснение:

На этом шаге извлекается подробная информация о каждой ключевой мутации:

Обрабатывается каждая строка варианта из файла SRR13018652.breast_cancer_genes.vcf, извлекая такие детали, как:

Хромосома ($1),

Позиция ($2),

Референсный и альтернативный аллели ($4 и $5).

Также включается поле аннотации эффекта вариантов (CSQ), которое содержит информацию о символе гена, типе мутации и её функциональном эффекте (например, миссенс-мутация, фреймшифт).

Результат выводится в удобочитаемом формате.

Шаг 25: Определение пола пациента

>echo "Determining patient sex:"
>y_reads=$(samtools idxstats /home/1/SRR13018652.recal.bam | awk '$1 == "chrY" {print $3}')
>echo "Y chromosome reads: $y_reads"

>if [ "$y_reads" -gt 1000 ]; then
    >echo "PATIENT SEX: MALE"
    >sex="MALE"
>else
    >echo "PATIENT SEX: FEMALE"
    >sex="FEMALE"
>fi

Объяснение:

На этом шаге определяется биологический пол пациента на основе наличия ридов Y-хромосомы в BAM-файле:

samtools idxstats предоставляет статистику по референсным хромосомам в BAM-файле.

Скрипт проверяет количество ридов, картированных на Y-хромосому (chrY):

Если количество ридов на Y-хромосоме больше 1000, пациент классифицируется как мужского пола (MALE).

В противном случае пациент классифицируется как женского пола (FEMALE).

Определённый пол (MALE или FEMALE) сохраняется в переменной sex.

Шаг 26: Подробный анализ ключевых мутаций

cat > /home/1/analyze_key_mutations.sh << 'EOF'
#!/bin/bash

echo "Detailed analysis of key mutations"

mutations=("chr17 7676154" "chr17 43063913" "chr3 179218303")

for mutation in "${mutations[@]}"; do
    IFS=" " read -r chrom pos <<< "$mutation"
    
    echo ""
    echo "Analyzing mutation: $chrom $pos"
    echo "------------------------"
    
    variant_line=$(awk -v chrom="$chrom" -v pos="$pos" '$1 == chrom && $2 == pos' /home/ser/1/SRR13018652.breast_cancer_genes.vcf)
    
    if [ -n "$variant_line" ]; then
        chrom=$(echo "$variant_line" | awk '{print $1}')
        pos=$(echo "$variant_line" | awk '{print $2}')
        ref=$(echo "$variant_line" | awk '{print $4}')
        alt=$(echo "$variant_line" | awk '{print $5}')
        qual=$(echo "$variant_line" | awk '{print $6}')
        filter=$(echo "$variant_line" | awk '{print $7}')
        info=$(echo "$variant_line" | awk '{print $8}')
        sample=$(echo "$variant_line" | awk '{print $10}')
        
        csq_info=$(echo "$variant_line" | grep -o "CSQ=[^;]*")
        first_csq=$(echo "$csq_info" | cut -d',' -f1)
        
        symbol=$(echo "$first_csq" | cut -d'|' -f4)
        consequence=$(echo "$first_csq" | cut -d'|' -f2)
        impact=$(echo "$first_csq" | cut -d'|' -f3)
        protein_change=$(echo "$first_csq" | cut -d'|' -f11)
        amino_acids=$(echo "$first_csq" | cut -d'|' -f12)
        
        echo "Gene: $symbol"
        echo "Change: $ref -> $alt" 
        echo "Type: $consequence"
        echo "Impact: $impact"
        if [ -n "$protein_change" ] && [ "$protein_change" != "" ]; then
            echo "Protein change: $protein_change"
        fi
        if [ -n "$amino_acids" ] && [ "$amino_acids" != "" ]; then
            echo "Amino acids: $amino_acids"
        fi
        echo "Filter: $filter"
        echo "Quality: $qual"
        
        dp=$(echo "$info" | grep -o "DP=[0-9]*" | cut -d'=' -f2)
        af=$(echo "$info" | grep -o "AF=[0-9.]*" | cut -d'=' -f2)
        tlod=$(echo "$info" | grep -o "TLOD=[0-9.]*" | cut -d'=' -f2)
        
        echo "Depth (INFO): $dp"
        echo "Allele frequency (INFO): $af"
        if [ -n "$tlod" ]; then
            echo "Tumor LOD: $tlod"
        fi
        
        gt=$(echo "$sample" | cut -d':' -f1)
        ad=$(echo "$sample" | cut -d':' -f2)
        sample_af=$(echo "$sample" | cut -d':' -f3)
        sample_dp=$(echo "$sample" | cut -d':' -f4)
        
        echo "Genotype: $gt"
        echo "Allelic Depths (AD): $ad"
        echo "Allele frequency (FORMAT): $sample_af"
        echo "Depth (FORMAT): $sample_dp"
        
        ref_count=$(echo "$ad" | cut -d',' -f1)
        alt_count=$(echo "$ad" | cut -d',' -f2)
        if [ -n "$ref_count" ] && [ -n "$alt_count" ] && [ "$ref_count" -gt 0 ]; then
            vaf=$(echo "scale=4; $alt_count / ($ref_count + $alt_count)" | bc)
            echo "VAF (calculated): $vaf"
        fi
        
        echo ""
        echo "Additional information:"
        echo "$info" | tr ';' '\n' | grep -E "^(DP|AF|TLOD|MBQ|MMQ|MPOS)="
        
    else
        echo "Mutation $chrom $pos not found"
    fi
    echo ""
done

EOF

chmod +x /home/1/analyze_key_mutations.sh
/home/1/analyze_key_mutations.sh > /home/1/key_mutations_detailed_analysis.txt


Объяснение:

Этот скрипт выполняет подробный анализ конкретных ключевых мутаций из списка генов, связанных с раком молочной железы. Он исследует определённые мутации (заданные по хромосоме и позиции) и извлекает детальную информацию о вариантах, такую как:

Символ гена, тип и значимость мутации.

Генотип и глубина секвенирования (AD, DP, AF).

Частота вариантного аллеля (VAF).

Дополнительные детали из поля INFO, такие как логарифм отношения правдоподобия для опухоли (TLOD), глубина покрытия (DP) и частота аллели (AF).

После запуска скрипта результат сохраняется в файл key_mutations_detailed_analysis.txt.

Шаг 27: Создание итогового клинического отчёта

>cat > /home/ser/1/final_clinical_report.md << EOF
># GENOMIC ANALYSIS REPORT
>## Whole Exome Sequencing - Tumor Only Analysis

>### PATIENT INFORMATION
>- Sample ID: SRR13018652
>- Biological Sex: $sex
>- Analysis Date: $(date +"%Y-%m-%d")

>### EXECUTIVE SUMMARY
>Whole exome sequencing identified three clinically significant mutations in key cancer genes.

>### KEY MUTATIONS IDENTIFIED

>#### 1. TP53 Mutation
>- Location: chr17:7676154
>- Variant: G>C (Missense)
>- VAF: 8.9%

>#### 2. BRCA1 Mutation
>- Location: chr17:43063913
>- Variant: G>C (Missense) 
>- VAF: 15.1%

>#### 3. PIK3CA Mutation
>- Location: chr3:179218303
>- Variant: G>A (Missense)
>- VAF: 36.9%

>### CLINICAL RECOMMENDATIONS

>#### Genetic Counseling & Testing
>1. Confirm germline status of BRCA1 mutation
>2. Family member testing recommended
>3. Discuss reproductive implications

>#### Cancer Surveillance
>1. Enhanced breast cancer screening
>2. $([ "$sex" = "MALE" ] && echo "Prostate cancer screening starting at age 40")
>3. Consider pancreatic cancer screening

>#### Therapeutic Considerations
>1. PARP inhibitors for BRCA1-related cancers
>2. PI3K inhibitors for PIK3CA-mutated cancers
>3. Clinical trial eligibility assessment

>### METHODS
>- Sequencing Technology: Whole Exome Sequencing
>- Reference Genome: hg38
>- Variant Caller: GATK Mutect2 (tumor-only mode)
>- Annotation: Ensembl VEP

>### LIMITATIONS
>- Tumor-only analysis without matched normal
- Germline vs somatic status undetermined
>>- Functional validation required for missense variants

>Analysis completed: $(date)
>EOF

Объяснение:

На этом шаге создаётся клинический отчёт в формате Markdown. Отчёт включает:

Информация о пациенте: Идентификатор образца и биологический пол.

Резюме: Краткий обзор ключевых мутаций, обнаруженных в ходе анализа (например, TP53, BRCA1, PIK3CA).

Клинические рекомендации: Предложения по генетическому консультированию, наблюдению за онкологическим статусом и возможным терапевтическим соображениям на основе выявленных мутаций.

Методы: Подробности об использованной технологии секвенирования, референсном геноме, программе для вызова вариантов и методах аннотации.

Ограничения: Указание на ограничения проведённого анализа, такие как подход без парного нормального образца (tumor-only).

Клинический отчёт сохраняется как final_clinical_report.md.

Шаг 28: Анализ рака простаты для пациентов мужского пола

>if [ "$sex" = "MALE" ]; then
    >cat > /home/1/prostate_cancer_genes.txt << 'EOF'
>BRCA1
>BRCA2
>ATM
>CHEK2
>TP53
>HOXB13
>MLH1
>MSH2
>MSH6
>PMS2
>EOF

    bcftools view -h /home/1/SRR13018652.annotated.vcf > /home/1/header_prostate.vcf
    
    bcftools view -H /home/1/SRR13018652.annotated.vcf | \
      grep -E "$(paste -s -d '|' /home/1/prostate_cancer_genes.txt)" | \
      cat /home1/header_prostate.vcf - > /home/1/SRR13018652.prostate_cancer_genes.vcf

    prostate_variants=$(grep -v "^#" /home/1/SRR13018652.prostate_cancer_genes.vcf | wc -l)
    echo "Variants in prostate cancer genes: $prostate_variants"
    
    echo "" >> /home/1/final_clinical_report.md
    echo "### PROSTATE CANCER RISK ASSESSMENT" >> /home/1/final_clinical_report.md
    echo "Based on BRCA1 and TP53 mutations:" >> /home/1/final_clinical_report.md
    echo "- Lifetime prostate cancer risk: 20-25%" >> /home/1/final_clinical_report.md
    echo "- Consider PSA screening starting at age 40" >> /home/1/final_clinical_report.md
>fi

Объяснение:
Этот шаг проверяет, является ли пациент мужчиной, и если да, анализирует мутации, связанные с раком простаты (в генах, таких как BRCA1, TP53 и др.). Он добавляет информацию о риске рака простаты и рекомендации по скринингу в итоговый клинический отчёт.

Шаг 29: Создание сводной таблицы результатов

>cat > /home/1/results_summary.csv << EOF
>Category,Value
>Sample_ID,SRR13018652
>Patient_Sex,$sex
>Total_Variants,$(grep -v "^#" /home/1/SRR13018652.final.vcf | wc -l)
>Breast_Cancer_Gene_Variants,$(grep -v "^#" /home/1/SRR13018652.breast_cancer_genes.vcf | wc -l)
Key_Driver_Mutations,3
>TP53_VAF,0.089
>BRCA1_VAF,0.151
>PIK3CA_VAF,0.369
>Analysis_Date,$(date +"%Y-%m-%d")

Объяснение:

Создание CSV-файла, суммирующего результаты анализа:

На этом шаге генерируется сводная таблица в формате CSV, которая включает ключевые метрики, такие как:

Идентификатор образца (SRR13018652).

Пол пациента ($sex).

Общее количество вариантов, найденных в финальном VCF (Total_Variants).

Количество вариантов, найденных в генах рака молочной железы (Breast_Cancer_Gene_Variants).

Количество ключевых драйверных мутаций (в данном случае жёстко задано как 3).

Частоты вариантных аллелей (VAF) для мутаций в TP53, BRCA1 и PIK3CA.

Дата анализа.

Динамические переменные:
Таблица извлекает динамические данные, такие как количество вариантов и дата анализа, выполняя команды внутри синтаксиса $(...) (например, подсчитывая строки в VCF-файлах для общего числа вариантов).

Созданный файл:
Файл results_summary.csv создаётся для быстрого обзора результатов анализа.

Шаг 30: Организация финальных результатов

>mkdir -p /home/1/final_results/01_vcf_files
>mkdir -p /home/1/final_results/02_reports
>mkdir -p /home/1/final_results/03_quality_metrics

>cp /home/1/SRR13018652.final.vcf /home/1/final_results/01_vcf_files/
>cp /home/1/SRR13018652.annotated.vcf /home/1/final_results/01_vcf_files/
>cp /home/1/SRR13018652.breast_cancer_genes.vcf /home/1/final_results/01_vcf_files/
>[ "$sex" = "MALE" ] && cp /home/1/SRR13018652.prostate_cancer_genes.vcf /home/1/final_results/01_vcf_files/

>cp /home/1/final_clinical_report.md /home/1/final_results/02_reports/
>cp /home/1/results_summary.csv /home/1/final_results/02_reports/
>cp /home/1/key_mutations_detailed_analysis.txt /home/1/final_results/02_reports/

>cp -r /home/1/fastqc_reports /home/1/final_results/03_quality_metrics/
>cp /home/1/fastp_report.html /home/1/final_results/03_quality_metrics/

>echo "ANALYSIS COMPLETED"
>echo "Results available in: /home/1/final_results/"
>echo "Clinical report: /home/1/final_results/02_reports/final_clinical_report.md"

Объяснение:

Создание директорий для финальных результатов:

mkdir -p /home/1/final_results/01_vcf_files создаёт директорию для VCF-файлов.
mkdir -p /home/1/final_results/02_reports создаёт директорию для отчётов.
mkdir -p /home/1/final_results/03_quality_metrics создаёт директорию для метрик качества (QC).

Организация файлов по соответствующим директориям:

VCF-файлы: Финальные VCF-файлы (final.vcf, annotated.vcf, breast_cancer_genes.vcf и, при наличии, prostate_cancer_genes.vcf для пациентов мужского пола) копируются в директорию 01_vcf_files.

Отчёты: Клинический отчёт (final_clinical_report.md), сводка результатов (results_summary.csv) и детальный анализ мутаций (key_mutations_detailed_analysis.txt) копируются в директорию 02_reports.

Метрики качества: Отчёты по контролю качества (fastqc_reports и fastp_report.html) копируются в директорию 03_quality_metrics.

Сообщение о завершении:

После организации всех финальных результатов по соответствующим директориям скрипт выводит сообщение о завершении с указанием путей к результатам.

Путь к финальным результатам:

Результаты можно найти в директории /home/1/final_results/, а клинический отчёт доступен по пути /home/1/final_results/02_reports/final_clinical_report.md.
