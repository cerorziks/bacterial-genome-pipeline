/*
========================================================================================
    AMR Detection Module
========================================================================================
*/

process DOWNLOAD_AMR_DB {
    publishDir "${params.outdir}/databases", mode: 'copy'
    storeDir "${params.db_cache}/amrfinder"
    debug true
    
    output:
    path "latest", emit: db
    
    script:
    """
    # Check if database already exists
    if [ -d "latest" ] && [ -f "latest/AMR.LIB" ]; then
        echo "AMRFinder database already exists in current directory. Skipping download." >&2
        echo "Database location: \$(pwd)/latest" >&2
        exit 0
    fi
    
    # Check for dated directories (alternative structure)
    EXISTING_DB=\$(ls -d 2* 2>/dev/null | head -n 1)
    if [ -n "\$EXISTING_DB" ] && [ -f "\$EXISTING_DB/AMR.LIB" ]; then
        echo "AMRFinder database already exists: \$EXISTING_DB" >&2
        echo "Creating 'latest' symlink to existing database..." >&2
        ln -sf "\$EXISTING_DB" latest
        exit 0
    fi
    
    echo "Starting AMR database download from NCBI..." >&2
    echo "Note: This can take 5-10 minutes on slower connections." >&2
    
    # Run update directly in the storeDir
    # We pipe to amrfinder_update.log for our watcher to read
    amrfinder_update -d . > amrfinder_update.log 2>&1 &
    UPDATE_PID=\$!
    
    # Watcher loop to report percentage and size
    echo "[Download Progress: 0%] (Connecting to NCBI...)" >&2
    LAST_PERCENT=0
    while kill -0 \$UPDATE_PID 2>/dev/null; do
        # Count all files in dated directories
        SIZE=\$(du -sk 2* 2>/dev/null | cut -f1 | awk '{s+=\$1} END {print s+0}')
        # 450000 KB is the approximate size of the full database
        PERCENT=\$((SIZE * 100 / 450000))
        if [ \$PERCENT -gt 100 ]; then PERCENT=100; fi
        
        # Only print if percentage changed (reduce spam)
        if [ \$PERCENT -ne \$LAST_PERCENT ]; then
            echo "[Download Progress: \$PERCENT%] (Downloaded: \${SIZE}KB / ~450,000KB)" >&2
            LAST_PERCENT=\$PERCENT
        fi
        sleep 5
    done
    
    wait \$UPDATE_PID
    EXIT_CODE=\$?
    
    if [ \$EXIT_CODE -eq 0 ]; then
        echo "Download completed successfully." >&2
    else
        echo "--------------------------------------------------------" >&2
        echo "ERROR: Database download failed (Exit Code: \$EXIT_CODE)" >&2
        echo "Last lines of the error log:" >&2
        tail -n 15 amrfinder_update.log >&2
        echo "--------------------------------------------------------" >&2
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
        # Build organism flag if provided
        ORG_FLAG=""
        if [ "${params.organism}" != "null" ] && [ -n "${params.organism}" ]; then
            ORG_FLAG="--organism ${params.organism}"
        fi

        amrfinder \
            --protein ${protein_fasta} \
            --threads $task.cpus \
            --database \$DB_PATH \
            --output ${sample}_amr.tsv \
            --plus \
            \$ORG_FLAG
    else
        echo "Error: No valid AMR database found for ${amr_db}"
        echo "Current directory contents:"
        ls -F
        exit 1
    fi
    
        echo "Resistome Profile for ${sample}:" > ${sample}_amr_summary.txt
        echo "======================================" >> ${sample}_amr_summary.txt
        
        for CAT in AMR VIRULENCE METAL BIOCIDE ACID HEAT; do
            count=$(awk -v cat="$CAT" -F'\t' '$9 ~ cat {count++} END {print count+0}' ${sample}_amr.tsv)
            if [ "$count" -gt 0 ]; then
                echo "" >> ${sample}_amr_summary.txt
                echo "[$CAT Resistance: $count genes]" >> ${sample}_amr_summary.txt
                awk -v cat="$CAT" -F'\t' '$9 ~ cat {print $6 " (" $11 ")"}' ${sample}_amr.tsv | sort -u | sed 's/^/  - /' >> ${sample}_amr_summary.txt
            fi
        done
        
        echo "" >> ${sample}_amr_summary.txt
        echo "Total genes detected: $(tail -n +2 ${sample}_amr.tsv | wc -l)" >> ${sample}_amr_summary.txt
    else
        echo "No AMR genes detected for ${sample}" > ${sample}_amr_summary.txt
    fi
    """
}
