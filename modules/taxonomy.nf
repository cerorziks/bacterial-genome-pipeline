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
    if [ ! -d "babykraken" ]; then
        mkdir babykraken
        curl -L "https://github.com/MDU-PHL/babykraken/raw/master/dist/babykraken.tar.gz" | tar xz -C babykraken --strip-components=1
    fi
    """
}

process KRAKEN2 {
    tag "$sample"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    path(db)
    
    output:
    tuple val(sample), path("${sample}_kraken2.report"), emit: report
    tuple val(sample), path("${sample}_kraken2.output"), emit: output
    
    script:
    """
    echo "[Taxonomy] Running Kraken2 for ${sample}..."
    
    # Robustly find the actual database directory containing the hash.k2d file
    # If find returns nothing, fallback to the staged dir itself
    TARGET_FILE=\$(find -L ${db} -name "hash.k2d" | head -n 1)
    
    if [ -z "\$TARGET_FILE" ]; then
        echo "WARNING: Could not find hash.k2d in ${db}. Listing contents:"
        ls -R ${db}
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
