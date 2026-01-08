/*
========================================================================================
    AMR Detection Module
========================================================================================
*/

process DOWNLOAD_AMR_DB {
    publishDir "${params.outdir}/databases", mode: 'copy'
    storeDir "${params.amr_db}"
    
    output:
    path "database", emit: db
    
    script:
    """
    mkdir -p database
    echo "Starting AMR database download..."
    
    # Run update in background and capture PID
    amrfinder_update -d database > amrfinder_update.log 2>&1 &
    UPDATE_PID=\$!
    
    # Watcher loop to report percentage
    while kill -0 \$UPDATE_PID 2>/dev/null; do
        SIZE=\$(du -sk database | cut -f1)
        # 450000 KB is the approximate size of the full database
        PERCENT=\$((SIZE * 100 / 450000))
        if [ \$PERCENT -gt 100 ]; then PERCENT=100; fi
        
        echo "[Download Progress: \$PERCENT%]"
        sleep 5
    done
    
    wait \$UPDATE_PID
    echo "Download completed successfully."
    """
}

process AMRFINDERPLUS {
    tag "$sample"
    publishDir "${params.outdir}/amr", mode: 'copy'
    
    input:
    tuple val(sample), path(protein_fasta)
    path amr_db
    
    output:
    tuple val(sample), path("${sample}_amr.tsv"), emit: results
    tuple val(sample), path("${sample}_amr_summary.txt"), emit: summary
    
    script:
    """
    # Use the provided database path
    # Handle the case where 'latest' symlink might not be created or is relative
    DB_PATH="${amr_db}/latest"
    if [ ! -d "\$DB_PATH" ]; then
        # Fallback: Find the actual dated directory
        DB_PATH=\$(ls -d ${amr_db}/2* 2>/dev/null | head -n 1)
    fi

    if [ -n "\$DB_PATH" ] && [ -d "\$DB_PATH" ]; then
        amrfinder \\
            --protein ${protein_fasta} \\
            --threads $task.cpus \\
            --database \$DB_PATH \\
            --output ${sample}_amr.tsv \\
            --plus
    else
        echo "Error: No valid AMR database found in ${amr_db}"
        exit 1
    fi
    
    # Create summary
    if [ -s ${sample}_amr.tsv ]; then
        echo "AMR Genes Found for ${sample}:" > ${sample}_amr_summary.txt
        echo "================================" >> ${sample}_amr_summary.txt
        tail -n +2 ${sample}_amr.tsv | cut -f6,11 | sort -u >> ${sample}_amr_summary.txt
        echo "" >> ${sample}_amr_summary.txt
        echo "Total AMR genes: \$(tail -n +2 ${sample}_amr.tsv | wc -l)" >> ${sample}_amr_summary.txt
    else
        echo "No AMR genes detected for ${sample}" > ${sample}_amr_summary.txt
    fi
    """
}
