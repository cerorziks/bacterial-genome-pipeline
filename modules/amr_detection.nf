/*
========================================================================================
    AMR Detection Module
========================================================================================
*/

process AMRFINDERPLUS {
    tag "$sample"
    publishDir "${params.outdir}/amr", mode: 'copy'
    
    input:
    tuple val(sample), path(protein_fasta)
    
    output:
    tuple val(sample), path("${sample}_amr.tsv"), emit: results
    tuple val(sample), path("${sample}_amr_summary.txt"), emit: summary
    
    script:
    """
    # Create local database directory and update
    mkdir -p amr_db
    amrfinder_update -d amr_db
    
    # Run AMRFinderPlus using the local database
    # Handle the case where 'latest' symlink might not be created
    DB_PATH="amr_db/latest"
    if [ ! -d "\$DB_PATH" ]; then
        # Fallback: Find the actual dated directory
        DB_PATH=\$(ls -d amr_db/2* 2>/dev/null | head -n 1)
    fi

    if [ -n "\$DB_PATH" ] && [ -d "\$DB_PATH" ]; then
        amrfinder \\
            --protein ${protein_fasta} \\
            --threads $task.cpus \\
            --database \$DB_PATH \\
            --output ${sample}_amr.tsv \\
            --plus
    else
        echo "Error: No valid AMR database found in amr_db/"
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
