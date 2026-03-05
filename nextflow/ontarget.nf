#!/usr/bin/env nextflow
// ontarget.nf

process chunkFasta {
    cache true
    debug false

    input:
    path subset_fasta
    path python_script
    val n_chunks
    
    output:
    path "gencode.${params.gencode_version}.transcripts.subset.chunk_*"

    script:
    """
    python ${python_script} ${subset_fasta} --n_chunks ${n_chunks}
    """
}

process testSmall {
    cache true
    debug false

    input:
    path chunk_file

    output:
    path "${chunk_file.name}"

    script:
    """
    seqkit sort -l ${chunk_file} | seqkit head -n 50 > tmp_small && mv tmp_small ${chunk_file.name}
    """
}

process setupTIGER {
    // Publish to a versioned directory to avoid conflicts
    publishDir "tiger_${params.tiger_version}", mode: 'copy', overwrite: true
    cache true
    debug false

    input:
    // Use version string as cache key instead of boolean trigger
    val tiger_version
    path tiger_pipeline_src

    output:
    path "tiger", emit: tiger_dir
    val true, emit: setup_complete

    script:
    """
    # Remove existing tiger directory if it exists
    rm -rf tiger

    # Clone the repository
    git clone ${params.tiger_hf} tiger

    # Download Git LFS files from Hugging Face (git lfs not available)
    cd tiger
    curl -L "https://huggingface.co/spaces/Knowles-Lab/tiger/resolve/main/model/saved_model.pb?download=true" -o model/saved_model.pb
    curl -L "https://huggingface.co/spaces/Knowles-Lab/tiger/resolve/main/model/variables/variables.data-00000-of-00001?download=true" -o model/variables/variables.data-00000-of-00001
    curl -L "https://huggingface.co/spaces/Knowles-Lab/tiger/resolve/main/model/variables/variables.index?download=true" -o model/variables/variables.index
    curl -L "https://huggingface.co/spaces/Knowles-Lab/tiger/resolve/main/model/fingerprint.pb?download=true" -o model/fingerprint.pb
    curl -L "https://huggingface.co/spaces/Knowles-Lab/tiger/resolve/main/model/keras_metadata.pb?download=true" -o model/keras_metadata.pb
    cd ..

    # Copy modified pipeline script (now tracked as input)
    cp ${tiger_pipeline_src} tiger/tiger_pipeline.py
    echo "TIGER ${tiger_version} is setup"
    """
}

process runTIGER {
    debug false
    cache 'deep'

    input:
    tuple val(chunk_id), path(chunk_file)
    path tiger_dir
    val setup_complete
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}_on_target.csv"), emit: chunk_csv
    
    script:
    """
    echo "[\$(date +%r)] Running TIGER on ${chunk_file}"

    # Prevent Python bytecode
    export PYTHONDONTWRITEBYTECODE=1

    # Prevent TensorFlow from writing to tiger_dir
    export TF_CPP_MIN_LOG_LEVEL=2
    export CUDA_CACHE_DISABLE=1

    # Keep any TF/Keras cache in current work dir, not tiger_dir
    export TFHUB_CACHE_DIR=\$PWD/.tfhub_cache
    export HF_HOME=\$PWD/.hf_cache

    # Symlink ALL files and directories from tiger_dir
    for item in ${tiger_dir}/*; do
        ln -s "\$item" . 2>/dev/null || true
    done

    python3 tiger_pipeline.py \
        --fasta_path "${chunk_file}" \
        --output "chunk_${chunk_id}_on_target.csv"

    echo "[\$(date +%r)] TIGER completed successfully for chunk ${chunk_id}"
    """
}

process filtHomopolymers {
    cache true

    input:
    tuple val(chunk_id), path(on_target_csv)
    path prefilter_script

    output:
    tuple val(chunk_id), path("*.pt.csv"), emit: chunks
    tuple val(chunk_id), path("*.pt.failed.csv"), emit: failed_csv

    script:
    basename = on_target_csv.baseName
    """
    echo "Starting prefilter (remove TTTT, AAAAA, CCCCC, GGGGG)"
    wc -l ${on_target_csv} | awk '{print \$1-1, "total guides"}'

    awk -f ${prefilter_script} \
        -v fails="${basename}.pt.failed.csv" \
        < ${on_target_csv} \
        > ${basename}.pt.csv

    wc -l ${basename}.pt.csv | awk '{print \$1-1, "lines remain after removing guides containing TTTT/V*5"}'
    awk -F',' 'NR > 1 {genes[\$6]++} END {print length(genes), "unique genes remain"}' ${basename}.pt.csv
    """
}


workflow ONTARGET {
    take:
    subset_fasta
    
    main:
    tiger_pipeline_src = file(params.tiger_src)

    ontarget_chunk_fasta_script = file(params.ontarget_chunk_fasta_script)
    chunks = chunkFasta(subset_fasta, ontarget_chunk_fasta_script, params.n_chunks)

    if (params.small_test) {
        chunks = chunks.flatten() | testSmall
    }

    // Pass version string and source file for proper caching
    // Change params.tiger_version to invalidate cache (e.g., "v1.0", "v1.1")
    tiger_setup = setupTIGER(params.tiger_version, tiger_pipeline_src)

    // Flatten to process each chunk individually
    chunk_channel = chunks
        .flatten()
        .map { file -> 
            def chunk_id = file.name.replaceAll(/.*chunk_/, '').replaceAll(/\..*/, '')
            tuple(chunk_id, file)
        }
    
    on_target_chunk_csvs = runTIGER(chunk_channel, tiger_setup.tiger_dir, tiger_setup.setup_complete)
    
    full_dataset_1 = on_target_chunk_csvs.chunk_csv
        .map { chunk_id, csv -> csv }
        .collectFile(
            name: "gencode.${params.gencode_version}.ot.tsv",
            storeDir: "datasets",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\d+/) { it.padLeft(10, '0') } }
        )

    // Filter away homopolymers and low TIGER guides
    ontarget_prefilter_poly_script = file(params.ontarget_prefilter_poly_script)
    pt = filtHomopolymers(on_target_chunk_csvs.chunk_csv, ontarget_prefilter_poly_script)

    full_dataset_2 = pt.chunks
        .map { chunk_id, csv -> csv }
        .collectFile(
            name: "gencode.${params.gencode_version}.ot.pt.tsv",
            storeDir: "datasets",
            keepHeader: true,
            sort: { file -> file.name.replaceAll(/\d+/) { it.padLeft(10, '0') } }
        )

    emit:
    ontarget_chunks = pt.chunks
}
