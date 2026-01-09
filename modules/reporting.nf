/*
========================================================================================
    Reporting Module
========================================================================================
*/

process COMPILED_HTML_REPORT {
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path(seqkit_stats)
    path(seqkit_assembly_stats)
    path(quast_stats)
    path(amr_stats)
    path(mlst_stats)
    path(virulence_stats)
    path(kraken_reports)
    path(pan_stats)
    path(tree)

    output:
    path("summary_report.html"), emit: html

    script:
    """
    generate_summary_report.py
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(files)

    output:
    path("multiqc_report.html"), emit: html
    path("*_data"), emit: data

    script:
    """
    multiqc --title "${params.multiqc_title}" --filename multiqc_report.html --force .
    """
}
