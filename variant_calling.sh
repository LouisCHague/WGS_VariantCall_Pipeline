#!/bin/bash

# Yeast Pipeline

# Download Data
echo "Downloading Data"

# Download paired reads from SRA
echo "Downloading SRA data..."
# Uncomment the lines below to download and convert SRA to FASTQ
prefetch SRR30661808
fastq-dump --split-files --gzip --outdir /app/pipeline_folder/reads ./SRR30661808/SRR30661808.sra

# Prep Files
echo "Downloading reference genome"
wget -P /app/pipeline_folder/supporting_files/r64 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
gunzip /app/pipeline_folder/supporting_files/r64/GCF_000146045.2_R64_genomic.fna.gz

echo "Indexing reference genome"
samtools faidx /app/pipeline_folder/supporting_files/r64/GCF_000146045.2_R64_genomic.fna

echo "Creating Dictionary File"
gatk CreateSequenceDictionary R=/app/pipeline_folder/supporting_files/r64/GCF_000146045.2_R64_genomic.fna O=/app/pipeline_folder/supporting_files/r64/GCF_000146045.2_R64_genomic.dict

# Directories for variant calling
ref="/app/pipeline_folder/supporting_files/r64/GCF_000146045.2_R64_genomic.fna"
aligned_reads="/app/pipeline_folder/aligned_reads"
reads="/app/pipeline_folder/reads"
results="/app/pipeline_folder/results"

# -----------------------
# STEP 1: QC - Run FastQC 
# -----------------------

echo "STEP 1: QC - Run fastqc"
fastqc ${reads}/SRR30661808_1.fastq.gz -o ${reads}/
fastqc ${reads}/SRR30661808_2.fastq.gz -o ${reads}/

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"
bwa index ${ref}
bwa mem -t 4 -R "@RG\tID:SRR30661808\tPL:ILLUMINA\tLB:T9\tSM:T9\tPU:SRR30661808" ${ref} ${reads}/SRR30661808_1.fastq.gz ${reads}/SRR30661808_2.fastq.gz > ${aligned_reads}/SRR30661808.paired.sam

# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR30661808.paired.sam -O ${aligned_reads}/SRR30661808_sorted_dedup_reads.bam

# -----------------------------------------------
# STEP 4: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

echo "STEP 4: Collect Alignment & Insert Size Metrics"
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR30661808_sorted_dedup_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR30661808_sorted_dedup_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

# -------------------------------------------
# STEP 5: Variant Calling using GATK HaplotypeCaller
# -------------------------------------------

echo "STEP 5: Variant Calling"
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR30661808_sorted_dedup_reads.bam -O ${results}/SRR30661808_variants.vcf
