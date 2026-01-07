/*
========================================================================================
    Reporting Module
========================================================================================
*/

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path(files)
    
    output:
    path("multiqc_report.html"), emit: html
    path("*_data"), emit: data
    
    script:
    """
    multiqc \\
        --title "${params.multiqc_title}" \\
        --filename multiqc_report.html \\
        --force \\
        .
    """
}
