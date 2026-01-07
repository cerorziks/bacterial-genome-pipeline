/*
========================================================================================
    Quality Control Module
========================================================================================
*/

process FASTQC {
    tag "$sample"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("*.html"), emit: html
    tuple val(sample), path("*.zip"), emit: zip
    
    script:
    """
    fastqc -q -t $task.cpus ${read1} ${read2}
    """
}

process FASTP {
    tag "$sample"
    publishDir "${params.outdir}/fastp", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    
    output:
    tuple val(sample), path("${sample}_R1_trimmed.fastq.gz"), path("${sample}_R2_trimmed.fastq.gz"), emit: reads
    tuple val(sample), path("${sample}.fastp.json"), emit: json
    tuple val(sample), path("${sample}.fastp.html"), emit: html
    
    script:
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample}_R1_trimmed.fastq.gz \\
        -O ${sample}_R2_trimmed.fastq.gz \\
        --thread $task.cpus \\
        --qualified_quality_phred ${params.fastp_qualified_quality_phred} \\
        --length_required ${params.fastp_min_length} \\
        --detect_adapter_for_pe \\
        --correction \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html
    """
}
