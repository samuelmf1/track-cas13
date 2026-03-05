// offtarget.nf

process BT_INDEX {
    tag "GENCODE ${params.gencode_version}"
    publishDir "datasets/index", mode: 'copy'
    cache true

    input:
    path transcripts
    path summary

    output:
    path "*.ebwt", emit: fm_index
    path "gencode.${params.gencode_version}.filtered.fa", emit: filtered_fa

    script:
    """
    awk -F'\\t' '\$7 == "scaffold_patch_haplotype" {print \$2}' $summary > remove.txt
    seqkit grep -v -f remove.txt --id-regexp "^([^|]+)" $transcripts -o gencode.${params.gencode_version}.filtered.fa.tmp
    seqkit grep -v -s -r -p '"N{23,}"' gencode.${params.gencode_version}.filtered.fa.tmp -o gencode.${params.gencode_version}.filtered.fa
    bowtie-build --threads ${task.cpus} gencode.${params.gencode_version}.filtered.fa gencode.${params.gencode_version}.transcripts
    """
}

process BT_ALIGN {
    // publishDir "datasets/chunks_${params.gencode_version}", mode: 'copy'

    input:
    tuple val(id), path(csv)
    path index

    output:
    tuple val(id), path("*.bt.bam"), path(csv), emit: bam_and_csv

    script:
    def base = csv.baseName.replaceAll(/\.ont$/, '')
    def idx = index.find { it.name.endsWith('.1.ebwt') }.name.replaceAll(/\.1\.ebwt$/, '')
    def script_fa = params.offtarget_create_fasta_script
    """
    awk -f $script_fa < $csv | \
    bowtie -a -v 2 -p ${task.cpus} --norc -x $idx -f - -S | \
    samtools view -b -F 20 - > ${base}.bt.bam
    """
}

process BT_FILTER {
    // publishDir "datasets/chunks_${params.gencode_version}", mode: 'copy'
    
    input:
    tuple val(id), path(bam), path(ont_csv)
    val level

    output:
    tuple val(id), path("*.bt.csv"), emit: chunks

    script:
    def base = bam.baseName.replaceAll(/\.bt$/, '')
    """
    set -e
    set -o pipefail
    
    mkdir -p ./tmp
    export TMPDIR=\$PWD/tmp
    export LC_ALL=C

    samtools view -F 4 $bam | \
    awk -f ${params.offtarget_condense_bt_script} -v mode=$level | \
    sort -t',' -k1,1n -k2,2n -T \$TMPDIR \
        --buffer-size=${params.sort_buffer_size} \
        --parallel=${task.cpus} \
        --compress-program=gzip > intermediate_unq.csv

    bash ${params.offtarget_postfilter_bt_script} $ont_csv intermediate_unq.csv $base
    # outputs ${base}.bt.csv (cleared) and ${base}.bt.ot.csv (contains off-target)
    """
}

process SASSY_SPLIT_REF {
    input:
    path ref

    output:
    path "chunks/*", emit: chunks

    script:
    """
    seqkit split2 --by-part 75 $ref -O chunks -f
    """
}

process SASSY_PREP {
    publishDir "datasets/sassy", mode: 'copy'
    debug true
    
    input:
    path "chunks/*"  // Accepts a directory/list of files

    output:
    path "patterns.txt", emit: patterns

    script:
    """
    export LC_ALL=C
    # Extract 3rd column from all chunks, skip header on each, and sort unique
    # Using 'tail -q -n +2' skips the first line of every file in the list
    tail -q -n +2 chunks/* | awk -F, '\$3 != "" {print \$3}' | sort -u > patterns.txt
    
    COUNT=\$(wc -l < patterns.txt)
    echo "\$COUNT unique guides to check with Sassy"
    """
}

process SASSY_SEARCH_CHUNK {
    tag "${chunk.baseName.split('_').last()}"
    
    input:
    path patterns
    path bin
    path chunk

    output:
    path "chunk_output.txt"

    script:
    """
    set -e
    set -o pipefail
    # Copy patterns to ensure fully isolated file for each task
    cp $patterns patterns_task.txt
    $bin search -j ${task.cpus} -k 2 --no-rc --max-n-frac 0 --pattern-file patterns_task.txt $chunk | cat > chunk_output.txt
    """
}

process SASSY_MERGE_CHUNKS {
    input:
    path "outputs/*"

    output:
    path "sassy_output.txt"

    script:
    """
    # 1. Grab the header from the first file
    head -n 1 \$(ls outputs/* | head -n 1) > sassy_output.txt

    # 2. Append all files, skipping their first lines to avoid duplicate headers
    tail -q -n +2 outputs/* >> sassy_output.txt
    """
}

process SASSY_BLACKLIST {
    publishDir "datasets/sassy", mode: 'copy'
    debug true
    
    input:
    path sassy_output
    path patterns_file

    output:
    path "blacklist_sassy.txt", emit: sassy_blacklist
    path "blacklist_sassy_seqs.txt", emit: sassy_blacklist_seqs

    script:
    def script_cigar_filt = params.offtarget_filter_cigar_script
    """
    set -e
    set -o pipefail

    # 1. Filter out some non-worrisome indel sequences from the Sassy output
    python3 $script_cigar_filt $sassy_output > blacklist_sassy.txt

    # 2. Add the sequences to the blacklist
    awk 'BEGIN {FS=OFS="\\t"} 
        NR==FNR { line[\$1]=\$0; next } 
        (FNR in line) { print line[FNR], \$1 }' blacklist_sassy.txt $patterns_file > blacklist_sassy_seqs.txt

    # 3. Clean format: Sequence [Space] ENSG_ID
    # We split the second column by "|" and grab the second part (ENSG)
    awk -F'\\t' '{
        split(\$2, ids, "|"); 
        print \$NF, ids[2]
    }' blacklist_sassy_seqs.txt > blacklist_sassy_seqs.txt.tmp
    mv blacklist_sassy_seqs.txt.tmp blacklist_sassy_seqs.txt
    """
}

process SASSY_FILTER {
    tag "$id"

    input:
    tuple val(id), path(bt_csv)
    path blacklist

    output:
    tuple val(id), path("${id}.bt.sa.csv"), emit: chunks

    script:
    def filter_script = params.offtarget_sassy_filter_script
    """
    awk -f $filter_script $blacklist $bt_csv > ${id}.bt.sa.csv
    """
}

process MERGE_BT_RESULTS {
    publishDir "datasets", mode: 'copy'

    input:
    path "inputs/*"

    output:
    path "gencode.${params.gencode_version}.ot.pt.bt.csv", emit: merged_csv

    script:
    """
    # Grab header from the first file alphabetically
    FIRST_FILE=\$(ls inputs/* | head -n 1)
    head -n 1 "\$FIRST_FILE" > gencode.${params.gencode_version}.ot.pt.bt.csv

    # Append all data, skipping headers, and sort if necessary 
    # (Using simple sort here to mimic collectFile's sort logic)
    tail -q -n +2 inputs/* >> gencode.${params.gencode_version}.ot.pt.bt.csv
    """
}

process MERGE_SA_RESULTS {
    publishDir "datasets", mode: 'copy'

    input:
    path "inputs/*"

    output:
    path "gencode.${params.gencode_version}.ot.pt.bt.sa.csv", emit: merged_csv

    script:
    """
    FIRST_FILE=\$(ls inputs/* | head -n 1)
    head -n 1 "\$FIRST_FILE" > gencode.${params.gencode_version}.ot.pt.bt.sa.csv
    tail -q -n +2 inputs/* >> gencode.${params.gencode_version}.ot.pt.bt.sa.csv
    """
}

workflow OFFTARGET {
    take:
    transcripts
    gencode_summary
    ontarget_chunks

    main:
    index_ch = BT_INDEX(transcripts, gencode_summary)
    bt_align = BT_ALIGN(ontarget_chunks, index_ch.fm_index.collect())
    bt_postfilt = BT_FILTER(bt_align.bam_and_csv, params.level)

    full_dataset_bt_ch = MERGE_BT_RESULTS(
        bt_postfilt.chunks.map { it[1] }.collect()
    )

    // sassy_ref_chunks = SASSY_SPLIT_REF(index_ch.filtered_fa)
    // sassy_prep = SASSY_PREP(bt_postfilt.chunks.map { it[1] }.collect())
    
    // sassy_search_results = SASSY_SEARCH_CHUNK(
    //     sassy_prep.patterns.collect(),
    //     params.sassy_bin,
    //     sassy_ref_chunks.chunks.flatten()
    // )

    // sassy_merged = SASSY_MERGE_CHUNKS(sassy_search_results.collect())
    // sassy_blacklist = SASSY_BLACKLIST(sassy_merged, sassy_prep.patterns)
    // sassy_postfilt = SASSY_FILTER(
    //     bt_postfilt.chunks,
    //     sassy_blacklist.sassy_blacklist_seqs
    // )

    // full_dataset_sa_ch = MERGE_SA_RESULTS(
    //     sassy_postfilt.chunks.map { it[1] }.collect()
    // )

    emit:
    cleared_chunks = bt_postfilt.chunks // change to sassy_postfilt.chunks
    // merged_final   = full_dataset_sa_ch.merged_csv
}