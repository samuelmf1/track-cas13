#!/usr/bin/env nextflow
// reduce.nf

process annotate {
    tag "${chunk_id}"
    
    input:
    tuple val(chunk_id), path(cleared_csv)
    path exons
    path cds_boundaries
    path gencode_summary
    path subset_expr
    path annot_script

    output:
    // tuple val(chunk_id), path("${cleared_csv.baseName}.a1.csv"),     emit: chunk_a1
    // tuple val(chunk_id), path("${cleared_csv.baseName}.a2.csv"),     emit: chunk_a2
    tuple val(chunk_id), path("${cleared_csv.baseName}.annots.csv"), emit: chunk_annots

    script:
    """
    python ${annot_script} \
        --csv ${cleared_csv} \
        --exons ${exons} \
        --cds ${cds_boundaries} \
        --summary ${gencode_summary} \
        --expr ${subset_expr} \
        --out_prefix ${cleared_csv.baseName}
    """
}

process collapse {
    tag "${chunk_id}"
    debug false

    input:
    tuple val(chunk_id), path(cleared_csv)

    output:
    tuple val(chunk_id), path("*.fin.csv"), emit: chunk
    tuple val(chunk_id), path("*.fin.collapsed.csv"), emit: collapsed_chunk
    tuple val(chunk_id), path("*.fin.collapsed.utr.csv"), emit: collapsed_utr
    tuple val(chunk_id), path("*.fin.collapsed.cds.csv"), emit: collapsed_cds

    script:
    basename = cleared_csv.baseName.replaceAll(/\.cleared\.bt\.sassy$/, '')
    """
    echo "Filtering genes/transcripts with fewer than 5 guides >23 bp apart"
    
    python ${params.reduce_collapse_script} \
        -g ${cleared_csv} \
        -o ${basename}.fin.csv
    """
}

process constitutive {
    publishDir "datasets/chunks_${params.gencode_version}", mode: 'copy', overwrite: true
    
    input:
    tuple val(chunk_id), path(guides_annot)
    path gencode_summary
    path constitutive_script
    
    output:
    tuple val(chunk_id), path("*.selected.constitutive.tsv"), emit: guides_annot_selected
    tuple val(chunk_id), path("*.selected.constitutive.summary.tsv"), emit: guides_annot_selected_summary
    tuple val(chunk_id), path("*.excluded.constitutive.tsv"), emit: guides_annot_excluded
    tuple val(chunk_id), path("*.excluded.constitutive.summary.tsv"), emit: guides_annot_excluded_summary

    script:
    def summary_file = "${guides_annot.baseName}.selected.constitutive.summary.tsv"
    def selected_guides = "${guides_annot.baseName}.selected.constitutive.tsv"
    def excluded_guides = "${guides_annot.baseName}.excluded.constitutive.tsv"
    def excluded_summary = "${guides_annot.baseName}.excluded.constitutive.summary.tsv"

    """
    python ${constitutive_script} \
        --tsv ${guides_annot} \
        --metadata ${gencode_summary} \
        --output ${selected_guides} \
        --lookahead 0

    echo "Selected ${selected_guides} / generated ${excluded_guides}"
    """
    // --sassy_ref ${params.sassy_ref} \
    // --sassy_cache ${params.sassy_cache}
}

process utr {
    publishDir "datasets/chunks_${params.gencode_version}", mode: 'copy', overwrite: true
    
    input:
    tuple val(chunk_id), path(guides_annot)
    path gencode_summary
    path utr_script
    
    output:
    tuple val(chunk_id), path("*.selected.utr.tsv"), emit: guides_annot_selected
    tuple val(chunk_id), path("*.selected.utr.summary.tsv"), emit: guides_annot_selected_summary
    tuple val(chunk_id), path("*.excluded.utr.tsv"), emit: guides_annot_excluded
    tuple val(chunk_id), path("*.excluded.utr.summary.tsv"), emit: guides_annot_excluded_summary

    script:
    def summary_file = "${guides_annot.baseName}.selected.utr.summary.tsv"
    def selected_guides = "${guides_annot.baseName}.selected.utr.tsv"
    def excluded_guides = "${guides_annot.baseName}.excluded.utr.tsv"
    def excluded_summary = "${guides_annot.baseName}.excluded.utr.summary.tsv"

    """
    python ${utr_script} \\
        --tsv ${guides_annot} \\
        --metadata ${gencode_summary} \\
        --output ${selected_guides} \\
        --lookahead 0

    echo "Selected ${selected_guides} / generated ${excluded_guides}"
    """
}



workflow REDUCE {
    take:
    transcripts_fa
    gtf_file
    exons
    cds_boundaries
    gencode_summary
    subset_expr
    cleared_chunks

    main:
    annot_script = file(params.reduce_annot_script)
    annotated_results = annotate(
        cleared_chunks,
        exons,
        cds_boundaries,
        gencode_summary,
        subset_expr,
        annot_script
    )

    final_res = collapse(
        annotated_results.chunk_annots
    )

    full_dataset = final_res.collapsed_chunk
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.ot.pt.bt.sa.collapsed.tsv",
            storeDir: "datasets",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\d+/) { it.padLeft(10, '0') } }
        )

    constitutive_script = file(params.reduce_constitutive_script)
    constitutive_chunks = constitutive(final_res.collapsed_cds, gencode_summary, constitutive_script)

    constitutive_library = constitutive_chunks.guides_annot_selected
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.selected.constitutive.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\d+/) { it.padLeft(10, '0') } }
        )

    constitutive_summary = constitutive_chunks.guides_annot_selected_summary
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.selected.constitutive.summary.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    constitutive_excluded = constitutive_chunks.guides_annot_excluded
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.excluded.constitutive.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    constitutive_excluded_summary = constitutive_chunks.guides_annot_excluded_summary
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.excluded.constitutive.summary.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    utr_script = file(params.reduce_utr_script ?: "\${projectDir}/src/reduce_utr_library.py")
    utr_chunks = utr(final_res.collapsed_utr, gencode_summary, utr_script)

    utr_library = utr_chunks.guides_annot_selected
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.selected.utr.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    utr_summary = utr_chunks.guides_annot_selected_summary
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.selected.utr.summary.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    utr_excluded = utr_chunks.guides_annot_excluded
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.excluded.utr.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    utr_excluded_summary = utr_chunks.guides_annot_excluded_summary
        .map { chunk_id, tsv -> tsv }
        .collectFile(
            name: "gencode.${params.gencode_version}.excluded.utr.summary.tsv",
            storeDir: "datasets/libraries",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\\d+/) { it.padLeft(10, '0') } }
        )

    emit:
    constitutive_library = constitutive_library
    constitutive_summary = constitutive_summary
    constitutive_excluded = constitutive_excluded
    constitutive_excluded_summary = constitutive_excluded_summary
    utr_library = utr_library
    utr_summary = utr_summary
    utr_excluded = utr_excluded
    utr_excluded_summary = utr_excluded_summary
}