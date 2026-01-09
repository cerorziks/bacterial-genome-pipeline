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
    kraken2 \\
        --db ${db} \\
        --threads $task.cpus \\
        --paired \\
        --report ${sample}_kraken2.report \\
        --output ${sample}_kraken2.output \\
        ${read1} ${read2}
    """
}
