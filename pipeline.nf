#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sra_id = "SRR30661808"
params.genome = "/mnt/c/Users/louis/nextflow_test/genome/GCF_000146045.2_R64_genomic.fna"

// Download paired-end reads from SRA
process PrefetchAndConvert {
    input:
    val sra_id

    output:
    path "./reads/${sra_id}_1.fastq.gz"
    path "./reads/${sra_id}_2.fastq.gz"

    script:
    """
    fastq-dump --split-files --gzip --outdir ./reads/ ${sra_id}
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



workflow {
    // Install paired-end DNA reads
    p_reads = PrefetchAndConvert(params.sra_id)
    // Index genome 
    p_index = indexGenome(params.genome)

}