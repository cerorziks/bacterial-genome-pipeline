/*
========================================================================================
    Reporting Module
========================================================================================
*/

process BOHRA_STYLE_SUMMARY {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path(seqkit_stats)
    path(quast_stats)
    path(amr_stats)
    path(mlst_stats)
    
    output:
    path("summary_mqc.tsv"), emit: tsv
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import csv

    summary = []
    
    # Header for MDU/Bohra style summary
    header = ['Sample', 'Total Reads', 'Total Bases (Mb)', 'N50 (bp)', 'Genome Size (Mb)', 'MLST/ST', 'AMR Genes Count']
    
    # This is a simplified parser - in a real scenario we'd be more robust
    # But for our pipeline structure, we can expect certain files
    
    # We'll collect data sample by sample if multiple... 
    # For now, let's assume one sample test or handle basic aggregation
    
    # [Bohra reports usually have a main table at the top]
    # We create a TSV for MultiQC custom content
    
    with open('summary_mqc.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\\t')
        writer.writerow(['# id: "bohra_summary"'])
        writer.writerow(['# section_name: "Executive Summary (MDU-Style)"'])
        writer.writerow(['# description: "Consolidated key metrics for all samples."'])
        writer.writerow(['# format: "tsv"'])
        writer.writerow(['# plot_type: "table"'])
        writer.writerow(header)
        
        # Logic to match samples across files would go here
        # For this pilot, we'll write a placeholder or attempt a basic join if possible
        # Since nextflow passes all files as a list, we have to find them
        
        # Placeholder/Simplified logic:
        # writer.writerow(['Example_Sample', '1.2M', '300', '120,000', '5.1', 'ST11', '12'])
        
        # If we have real files, we'd parse them here.
        # Since this is a template for the workflow:
        writer.writerow(['Summary Metrics', 'Check individual sections below', '-', '-', '-', '-', '-'])

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
    multiqc \\
        --title "${params.multiqc_title}" \\
        --filename multiqc_report.html \\
        --force \\
        .
    """
}
