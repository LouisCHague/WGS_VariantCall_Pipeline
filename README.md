# WGS Variant Calling Pipeline
This pipeline performs germline variant calling for Whole Genome Sequencing (WGS) data, specifically targeting short variant discovery of SNPs and Indels. The workflow is designed to follow the [GATK Best Practices for variant calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels), focusing on a single-sample analysis pipeline that aligns, pre-processes, and generates a final VCF file of variant calls using GATK HaplotypeCaller. This pipeline has been tested on [Baker's yeast paired-end sequencing data](https://www.ncbi.nlm.nih.gov/sra/SRX26080880[accn]).

## Workflow
The pipeline includes the following main steps (I was unable to perform base quality recalibration as I could not find a variant file for Baker's Yeast):

* Quality Control: Evaluate raw sequencing data quality using FastQC.
* Alignment: Map trimmed reads to a reference genome using BWA-MEM.
* Duplicate Marking: Mark duplicate reads using GATK MarkDuplicatesSpark.
* Variant Calling: Detect SNPs and Indels using GATK HaplotypeCaller.

## Installation

### Clone the repository:
git clone https://github.com/username/WGS-variant-calling-pipeline.git
cd WGS-variant-calling-pipeline

#### Requirements
* GATK: Version 4.x or higher
* FastQC: Version 0.11.x or higher
* BWA: Version 0.7.17 or higher
* Samtools: Version 1.9 or higher

### Optional: Use the provided Docker container to avoid installing software dependencies individually:
docker build -t gatk_pipeline .
docker run -v C:\Users\louis\gatk_pipeline\pipeline_folder:/app/pipeline_folder -it gatk_pipeline

## Additional Information
Pipeline_folder should contain the directories: aligned_reads, data, reads, results, scripts, supporting_files.



