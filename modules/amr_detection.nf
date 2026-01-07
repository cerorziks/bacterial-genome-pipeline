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
    # Create local database directory
    mkdir amr_db
    
    # Update AMRFinderPlus database in the local directory
    amrfinder_update -d amr_db || true
    
    # Run AMRFinderPlus using the local database
    # If update failed, it might still find a database if the container has one, 
    # but usually we need to point to the local one
    DB_OPTS=""
    if [ -d "amr_db/latest" ]; then
        DB_OPTS="-d amr_db/latest"
    fi
    
    amrfinder \\
        --protein ${protein_fasta} \\
        --threads $task.cpus \\
        --output ${sample}_amr.tsv \\
        \$DB_OPTS \\
        --plus
    
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
