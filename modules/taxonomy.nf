/*
========================================================================================
    Taxonomy Module
========================================================================================
*/

process DOWNLOAD_KRAKEN2_DB {
    storeDir "${params.db_cache}/kraken2"
    
    output:
    path("babykraken"), emit: db
    
    script:
    """
    # Check if babykraken database already exists
    if [ -d "babykraken" ] && [ -f "babykraken/hash.k2d" ]; then
        echo "Kraken2 (BabyKraken) database already exists. Skipping download." >&2
        echo "Database location: \$(pwd)/babykraken" >&2
        exit 0
    fi
    
    echo "Downloading BabyKraken database..." >&2
    mkdir -p babykraken
    curl -L "https://github.com/MDU-PHL/babykraken/raw/master/dist/babykraken.tar.gz" | tar xz -C babykraken --strip-components=1
    echo "BabyKraken download complete." >&2
    """
}

process KRAKEN2 {
    tag "$sample"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    stageInMode 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    path(db)
    
    output:
    tuple val(sample), path("${sample}_kraken2.report"), emit: report
    tuple val(sample), path("${sample}_kraken2.output"), emit: output
    
    script:
    """
    echo "[Taxonomy] Running Kraken2 for ${sample}..."
    
    # 1. Set the DB path variable to current dir as a safety net
    export KRAKEN2_DB_PATH=.

    # 2. Robustly find the actual database directory containing the hash.k2d file
    #    Broaden search to catch any file ending in hash.k2d
    TARGET_FILE=\$(find -L ${db} -name "*hash.k2d" | head -n 1)
    
    if [ -z "\$TARGET_FILE" ]; then
        echo "WARNING: Could not find *hash.k2d in ${db}. Listing contents for debug:"
        ls -laR ${db}
        # Fallback to the staged directory path
        DB_DIR=${db}
    else
        DB_DIR=\$(dirname \$TARGET_FILE)
    fi
    
    echo "Using Kraken2 DB at: \$DB_DIR"
    
    kraken2 \\
        --db \$DB_DIR \\
        --threads $task.cpus \\
        --paired \\
        --report ${sample}_kraken2.report \\
        --output ${sample}_kraken2.output \\
        ${read1} ${read2}
    """
}
