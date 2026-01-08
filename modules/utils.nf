/*
========================================================================================
    Utility Module (SeqKit)
========================================================================================
*/

process SEQKIT_STATS {
    tag "$sample"
    publishDir "${params.outdir}/seqkit", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_stats.txt"), emit: stats
    
    script:
    """
    seqkit stats ${read1} ${read2} -T -a > ${sample}_stats.txt
    """
}

process SEQKIT_STATS_ASSEMBLY {
    tag "$sample"
    publishDir "${params.outdir}/seqkit", mode: 'copy'
    
    input:
    tuple val(sample), path(assembly)
    
    output:
    tuple val(sample), path("${sample}_assembly_stats.txt"), emit: stats
    
    script:
    """
    seqkit stats ${assembly} -T -a > ${sample}_assembly_stats.txt
    """
}
