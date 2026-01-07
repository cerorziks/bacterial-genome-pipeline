/*
========================================================================================
    Annotation Module
========================================================================================
*/

process PROKKA {
    tag "$sample"
    publishDir "${params.outdir}/annotation", mode: 'copy'
    
    input:
    tuple val(sample), path(assembly)
    
    output:
    tuple val(sample), path("${sample}"), emit: results
    tuple val(sample), path("${sample}/${sample}.faa"), emit: faa
    tuple val(sample), path("${sample}/${sample}.ffn"), emit: ffn
    tuple val(sample), path("${sample}/${sample}.gff"), emit: gff
    tuple val(sample), path("${sample}/${sample}.gbk"), emit: gbk
    tuple val(sample), path("${sample}/${sample}.txt"), emit: txt
    
    script:
    """
    prokka \\
        --outdir ${sample} \\
        --prefix ${sample} \\
        --cpus $task.cpus \\
        --force \\
        --mincontiglen 200 \\
        --centre XXX \\
        --compliant \\
        ${assembly}
    """
}
