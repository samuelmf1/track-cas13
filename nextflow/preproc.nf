#!/usr/bin/env nextflow
// preproc.nf

process findCDSBounds {
    publishDir "datasets/annotations", mode: 'copy', overwrite: true
    cache true
    debug true

    input:
    path preproc_cds_boundaries_script
    path gtf_file

    output:
    path 'transcript_region_boundaries.tsv', emit: cds_boundaries
    path 'transcript_region_boundaries_validation_report.txt', emit: validation_report

    script:
    """
    python ${preproc_cds_boundaries_script} \
    -i ${gtf_file} \
    -o transcript_region_boundaries.tsv \
    --validate --validation-report transcript_region_boundaries_validation_report.txt
    """
}

process findExonJuncts {
    publishDir "datasets/annotations", mode: 'copy', overwrite: true

    input:
    path preproc_exon_script
    path gtf_file
    path transcripts_fa

    output:
    path 'exon_sequences.tsv', emit: exons

    script:
    """
    python ${preproc_exon_script} -f ${transcripts_fa} -g ${gtf_file} -o exon_sequences.tsv
    """
}

process filtFASTA {
    publishDir "datasets/filtered", mode: 'copy', overwrite: true
    cache true
    debug true

    input:
    path transcripts_fa
    path cds_boundaries
    path gtf_file
    path preproc_reduce_fasta_script
    path preproc_reduce_fasta_validation_script

    output:
    path "gencode.${params.gencode_version}.transcripts.subset.fasta", emit: subset_fasta
    path "gencode.${params.gencode_version}.transcripts.subset.summary.tsv", emit: summary
    path "validation.log", emit: validation_log

    script:
    """
    python ${preproc_reduce_fasta_script} \
      ${transcripts_fa} \
      gencode.${params.gencode_version}.transcripts.subset.fasta \
      -b ${cds_boundaries} \
      -g ${gtf_file} \
      -s gencode.${params.gencode_version}.transcripts.subset.summary.tsv

    python ${preproc_reduce_fasta_validation_script} \
      gencode.${params.gencode_version}.transcripts.subset.fasta \
      gencode.${params.gencode_version}.transcripts.subset.summary.tsv
    """
}

process expression {
    publishDir "datasets/expression", mode: 'copy', overwrite: true
    debug false
    
    input:
    path expr_file
    path model_file
    path transcript_summary_file
    path preproc_expr_script
    
    output:
    path "gencode.${params.gencode_version}.transcripts.subset.expr.tsv", emit: subset_expr
    path "gencode.${params.gencode_version}.transcripts.subset.expr.log", emit: expression_log
    
    script:
    """
    python ${preproc_expr_script} \
      --expr ${expr_file} \
      --model ${model_file} \
      --transcript-summary ${transcript_summary_file} \
      --skip-models
    """
}

workflow PREPROC {
    take:
    gtf_file
    transcripts_fa

    main:
    preproc_cds_boundaries_script = file(params.preproc_cds_boundaries_script)
    cds_boundaries_ch = findCDSBounds(preproc_cds_boundaries_script, gtf_file)

    preproc_exon_script = file(params.preproc_exon_script)
    exons = findExonJuncts(preproc_exon_script, gtf_file, transcripts_fa)

    preproc_reduce_fasta_script = file(params.preproc_reduce_fasta_script)
    preproc_reduce_fasta_validation_script = file(params.preproc_reduce_fasta_validation_script)
    subset_fasta_ch = filtFASTA(transcripts_fa, cds_boundaries_ch.cds_boundaries, gtf_file, preproc_reduce_fasta_script, preproc_reduce_fasta_validation_script)

    preproc_expr_script = file(params.preproc_expr_script)
    model_file = file(params.model_file)
    expr_file = file(params.expr_file)
    subset_expr_ch = expression(expr_file, model_file, subset_fasta_ch.summary, preproc_expr_script)
    
    emit:
    subset_fasta = subset_fasta_ch.subset_fasta
    gencode_summary = subset_fasta_ch.summary
    cds_boundaries = cds_boundaries_ch.cds_boundaries
    subset_expr = subset_expr_ch.subset_expr
    exons = exons
}