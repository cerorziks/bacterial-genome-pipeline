/*
========================================================================================
    MLST Module
========================================================================================
*/

process MLST {
    tag "$sample"
    publishDir "${params.outdir}/mlst", mode: 'copy'
    
    input:
    tuple val(sample), path(assembly)
    
    output:
    tuple val(sample), path("${sample}_mlst.tsv"), emit: results
    
    script:
    """
    mlst \\
        --threads $task.cpus \\
        ${assembly} > ${sample}_mlst.tsv
    """
}
