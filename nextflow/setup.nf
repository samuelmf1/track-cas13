#!/usr/bin/env nextflow
// setup.nf

process loadTranscriptomeFASTA {
    tag "GENCODE ${params.gencode_version}"
    debug false
    publishDir "datasets/annotations", mode: 'copy', overwrite: true
    
    output:
    path "gencode.${params.gencode_version}.transcripts.fa", emit: transcripts_fa 

    script:
    def version_num = params.gencode_version.replaceAll('v', '')
    """
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version_num}/gencode.${params.gencode_version}.transcripts.fa.gz
    gunzip gencode.${params.gencode_version}.transcripts.fa.gz
    echo "Downloaded gencode.${params.gencode_version}.transcripts.fa"
    seqkit seq -u gencode.${params.gencode_version}.transcripts.fa > gencode.${params.gencode_version}.transcripts.uc.fa
    echo " - Uppercased all nucleotides for use as Sassy reference"
    mv gencode.${params.gencode_version}.transcripts.uc.fa gencode.${params.gencode_version}.transcripts.fa
    """
}


process loadGTF {
    tag "GENCODE ${params.gencode_version}"
    debug false
    publishDir "datasets/annotations", mode: 'copy', overwrite: true

    output:
    path "gencode.${params.gencode_version}.annotation.gtf", emit: annotation_gtf
    
    script:
    def version_num = params.gencode_version.replaceAll('v', '')
    """
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${version_num}/gencode.${params.gencode_version}.annotation.gtf.gz
    gzip -dc gencode.${params.gencode_version}.annotation.gtf.gz > gencode.${params.gencode_version}.annotation.gtf
    echo "Downloaded gencode.${params.gencode_version}.annotation.gtf"
    """
}


workflow SETUP {
    main:
    transcripts_ch = loadTranscriptomeFASTA()
    annotation_ch = loadGTF()
    
    emit:
    gtf = annotation_ch
    transcripts = transcripts_ch
}
