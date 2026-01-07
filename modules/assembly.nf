/*
========================================================================================
    Assembly Module
========================================================================================
*/

process SPADES {
    tag "$sample"
    publishDir "${params.outdir}/assembly/spades", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}/contigs.fasta"), emit: assembly
    tuple val(sample), path("${sample}/scaffolds.fasta"), emit: scaffolds
    tuple val(sample), path("${sample}/spades.log"), emit: log
    
    script:
    def kmers = params.spades_kmers == 'auto' ? '' : "-k ${params.spades_kmers}"
    """
    spades.py \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o ${sample} \\
        --threads $task.cpus \\
        --memory ${task.memory.toGiga()} \\
        ${kmers} \\
        --careful
    """
}

process SHOVILL {
    tag "$sample"
    publishDir "${params.outdir}/assembly/shovill", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}/contigs.fasta"), emit: assembly
    tuple val(sample), path("${sample}/shovill.log"), emit: log
    
    script:
    """
    shovill \\
        --R1 ${read1} \\
        --R2 ${read2} \\
        --outdir ${sample} \\
        --cpus $task.cpus \\
        --ram ${task.memory.toGiga()} \\
        --force
        
    # Standardize output name
    mv ${sample}/contigs.fa ${sample}/contigs.fasta
    """
}

process QUAST {
    publishDir "${params.outdir}/quast", mode: 'copy'
    
    input:
    path(assemblies)
    
    output:
    path("quast_results"), emit: results
    path("quast_results/report.html"), emit: html
    path("quast_results/report.tsv"), emit: tsv
    
    script:
    """
    quast.py \\
        ${assemblies} \\
        -o quast_results \\
        --threads $task.cpus \\
        --min-contig 500
    """
}
