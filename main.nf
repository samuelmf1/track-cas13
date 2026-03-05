#!/usr/bin/env nextflow

include { SETUP } from './nextflow/setup.nf'
include { PREPROC } from './nextflow/preproc.nf'
include { ONTARGET } from './nextflow/ontarget.nf'
include { OFFTARGET } from './nextflow/offtarget.nf'
include { REDUCE } from './nextflow/reduce.nf'

workflow {

    setup_ch = SETUP()

    preproc_ch = PREPROC(
        setup_ch.gtf,
        setup_ch.transcripts
    )

    ontarget_ch = ONTARGET(
        preproc_ch.subset_fasta
    )

    offtarget_ch = OFFTARGET(
        setup_ch.transcripts,
        preproc_ch.gencode_summary,
        ontarget_ch.ontarget_chunks
    )

    reduce_ch = REDUCE(
        setup_ch.transcripts,
        setup_ch.gtf,
        preproc_ch.exons,
        preproc_ch.cds_boundaries,
        preproc_ch.gencode_summary,
        preproc_ch.subset_expr,
        offtarget_ch.cleared_chunks
    )

}