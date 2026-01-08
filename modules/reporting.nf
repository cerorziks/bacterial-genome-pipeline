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
    import glob

    # Data dictionaries to store metrics by sample name
    data = {}

    def get_sample_name(filename, suffix):
        name = os.path.basename(filename)
        if name.endswith(suffix):
            return name[:-len(suffix)]
        return name

    # 1. Parse SeqKit for Reads and Bases
    for f in glob.glob("*_stats.txt"):
        sname = get_sample_name(f, "_stats.txt")
        with open(f) as fh:
            lines = fh.readlines()
            if len(lines) > 1:
                cols = lines[1].split('\\t')
                if sname not in data: data[sname] = {}
                data[sname]['reads'] = "{:,}".format(int(cols[3]))
                data[sname]['bases'] = "{:.1f}".format(int(cols[4])/1000000)

    # 2. Parse MLST for ST
    for f in glob.glob("*_mlst.tsv"):
        sname = get_sample_name(f, "_mlst.tsv")
        with open(f) as fh:
            line = fh.readline().strip()
            if line:
                cols = line.split('\\t')
                if sname not in data: data[sname] = {}
                data[sname]['st'] = cols[2] if len(cols) > 2 else "-"

    # 3. Parse AMR for Gene Count
    for f in glob.glob("*_amr.tsv"):
        sname = get_sample_name(f, "_amr.tsv")
        with open(f) as fh:
            lines = fh.readlines()
            count = max(0, len(lines) - 1)
            if sname not in data: data[sname] = {}
            data[sname]['amr'] = str(count)

    # 4. Parse QUAST for N50 and Genome Size
    if os.path.exists("report.tsv"):
        with open("report.tsv") as fh:
            lines = [l.strip().split('\\t') for l in fh.readlines()]
            if len(lines) > 0:
                header = lines[0]
                for i in range(1, len(header)):
                    # Match by finding the closest sample name in data
                    qname = header[i]
                    # Direct match first
                    target = qname
                    if target not in data:
                        # Try to find a match that is a prefix
                        for dname in data.keys():
                            if qname.startswith(dname) or dname.startswith(qname):
                                target = dname
                                break
                    
                    if target not in data: data[target] = {}
                    for row in lines:
                        if row[0] == "N50": data[target]['n50'] = "{:,}".format(int(row[i]))
                        if row[0] == "Total length": data[target]['gsize'] = "{:.2f}".format(int(row[i])/1000000)

    # Write the MultiQC TSV
    with open('summary_mqc.tsv', 'w') as f:
        f.write('# id: "bohra_summary"\\n')
        f.write('# section_name: "Executive Summary (MDU-Style)"\\n')
        f.write('# description: "Consolidated key metrics for all samples."\\n')
        f.write('# format: "tsv"\\n')
        f.write('# plot_type: "table"\\n')
        
        headers = ['Sample', 'Total Reads', 'Total Bases (Mb)', 'N50 (bp)', 'Genome Size (Mb)', 'MLST/ST', 'AMR Genes']
        f.write('\\t'.join(headers) + '\\n')
        
        for sname in sorted(data.keys()):
            row = [
                sname,
                data[sname].get('reads', '-'),
                data[sname].get('bases', '-'),
                data[sname].get('n50', '-'),
                data[sname].get('gsize', '-'),
                data[sname].get('st', '-'),
                data[sname].get('amr', '0')
            ]
            f.write('\\t'.join(row) + '\\n')
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
