/*
========================================================================================
    AMR Detection Module
========================================================================================
*/

process DOWNLOAD_AMR_DB {
    publishDir "${params.outdir}/databases", mode: 'copy'
    storeDir "${params.amr_db}"
    debug true
    
    output:
    path "latest", emit: db
    
    script:
    """
    echo "Starting/Resuming AMR database download from NCBI..."
    echo "Note: This can take 5-10 minutes on slower connections."
    
    # Run update directly in the storeDir
    # We pipe to amrfinder_update.log for our watcher to read
    amrfinder_update -d . > amrfinder_update.log 2>&1 &
    UPDATE_PID=\$!
    
    # Watcher loop to report percentage and size
    echo "[Download Progress: 0%] (Connecting to NCBI...)"
    while kill -0 \$UPDATE_PID 2>/dev/null; do
        # Count all files in dated directories
        SIZE=\$(du -sk 2* 2>/dev/null | cut -f1 | awk '{s+=\$1} END {print s+0}')
        # 450000 KB is the approximate size of the full database
        PERCENT=\$((SIZE * 100 / 450000))
        if [ \$PERCENT -gt 100 ]; then PERCENT=100; fi
        
        echo "[Download Progress: \$PERCENT%] (Downloaded: \${SIZE}KB / ~450,000KB)"
        sleep 10
    done
    
    wait \$UPDATE_PID
    EXIT_CODE=\$?
    
    if [ \$EXIT_CODE -eq 0 ]; then
        echo "Download completed successfully."
    else
        echo "--------------------------------------------------------"
        echo "ERROR: Database download failed (Exit Code: \$EXIT_CODE)"
        echo "Last lines of the error log:"
        tail -n 15 amrfinder_update.log
        echo "--------------------------------------------------------"
        exit \$EXIT_CODE
    fi
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
    echo "[AMR] Starting AMRFinderPlus for sample: ${sample}"
    # Use the provided database path
    # We check if the input is the base directory or the actual DB directory
    if [ -d "${amr_db}/latest" ]; then
        DB_PATH="${amr_db}/latest"
    elif [ -f "${amr_db}/AMR.LIB" ]; then
        # The input itself is the database directory
        DB_PATH="${amr_db}"
    else
        # Fallback to current directory if it's staged as a flat list, 
        # or look for dated directories
        DB_PATH=\$(ls -d ${amr_db}/2* 2>/dev/null | head -n 1)
        if [ -z "\$DB_PATH" ]; then
            DB_PATH="${amr_db}"
        fi
    fi

    if [ -n "\$DB_PATH" ] && [ -d "\$DB_PATH" ]; then
        amrfinder \\
            --protein ${protein_fasta} \\
            --threads $task.cpus \\
            --database \$DB_PATH \\
            --output ${sample}_amr.tsv \\
            --plus
    else
        echo "Error: No valid AMR database found for ${amr_db}"
        echo "Current directory contents:"
        ls -F
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
