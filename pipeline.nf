#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Data
params.reads = "$projectDir/reads/SRR*_{1,2}.fastq.gz"
params.genome = "$projectDir/genome/GCF_000146045.2_R64_genomic.fna"

// FastQC Reads
process FastQC {
    input:
    tuple val(sra_id), path(pair_path)

    output:
    path "fastqc_${sra_id}_logs"

    script:
    """
    fastq.sh "$sra_id" "$pair_path"
    """
}

// Index genome
process indexGenome {
    input:
    path genome_path

    output:
    path "${genome_path}.fai"

    script:
    """
    samtools faidx $genome_path
    """
}

// Generate genome dictionary
process createDict {
    input:
    path genome_path

    output:
    path "${genome_path}.dict"

    script:
    """
    gatk CreateSequenceDictionary R=${genome_path} O=${genome_path}.dict
    """
}

workflow {
    Channel
    .fromFilePairs(params.reads, checkIfExists:true)
    .set { read_pairs_ch }

    //read_pairs_ch.view()
    FastQC(read_pairs_ch)

    // Index genome and create dictionary
    indexGenome(params.genome)
    createDict(params.genome)
}
