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
    '''
    #!/usr/bin/env python3
    import os
    import glob
    import datetime
    import sys

    # Helper for logging
    def log(msg):
        print(f"[Report Gen] {msg}", file=sys.stderr)

    log("Starting report generation...")

    # Data dictionaries to store metrics by sample name
    data = {}
    pan_data = {}

    def get_sample_name(filename, suffix):
        name = os.path.basename(filename)
        if name.endswith(suffix):
            return name[:-len(suffix)]
        return name

    # 1. Parse SeqKit for Raw Reads
    stats_files = glob.glob("*_stats.txt")
    log(f"Found {len(stats_files)} seqkit stats files")
    for f in stats_files:
        if "_assembly_stats.txt" in f: continue
        sname = get_sample_name(f, "_stats.txt")
        with open(f) as fh:
            lines = fh.readlines()
            if len(lines) > 1:
                cols = lines[1].split('\t')
                if sname not in data: data[sname] = {}
                data[sname]['reads'] = "{:,}".format(int(cols[3]))
                data[sname]['bases'] = "{:.1f}".format(int(cols[4])/1000000)

    # 2. Parse MLST
    mlst_files = glob.glob("*_mlst.tsv")
    for f in mlst_files:
        sname = get_sample_name(f, "_mlst.tsv")
        with open(f) as fh:
            line = fh.readline().strip()
            if line:
                cols = line.split('\t')
                if sname not in data: data[sname] = {}
                data[sname]['st'] = cols[2] if len(cols) > 2 else "-"
                data[sname]['scheme'] = cols[1] if len(cols) > 1 else "-"

    # 3. Parse AMR (Categorized)
    amr_files = glob.glob("*_amr.tsv")
    for f in amr_files:
        sname = get_sample_name(f, "_amr.tsv")
        if sname not in data: data[sname] = {}
        counts = {'AMR': 0, 'VIRULENCE': 0, 'METAL': 0, 'BIOCIDE': 0, 'ACID': 0, 'HEAT': 0, 'OTHER': 0}
        categorized_genes = {k: [] for k in counts.keys()}
        with open(f) as fh:
            header = fh.readline()
            for l in fh:
                c = l.strip().split('\t')
                if len(c) > 10:
                    gene, etype, subclass = c[5], c[8].upper(), c[10]
                    cat = 'OTHER'
                    if 'AMR' in etype: cat = 'AMR'
                    elif 'VIRULENCE' in etype: cat = 'VIRULENCE'
                    elif 'METAL' in etype: cat = 'METAL'
                    elif 'BIOCIDE' in etype: cat = 'BIOCIDE'
                    elif 'ACID' in etype: cat = 'ACID'
                    elif 'HEAT' in etype: cat = 'HEAT'
                    categorized_genes[cat].append(f"{gene} ({subclass})")
                    counts[cat] += 1
        data[sname]['amr_counts'] = counts
        data[sname]['amr_categorized'] = {k: ", ".join(sorted(set(v))) for k, v in categorized_genes.items() if v}
        data[sname]['amr_total'] = sum(counts.values())

    # 4. Parse Virulence
    vir_files = glob.glob("*_virulence_summary.txt")
    for f in vir_files:
        sname = get_sample_name(f, "_virulence_summary.txt")
        with open(f) as fh:
            content = fh.read()
            if "Total virulence factors" in content:
                count = content.split("Total virulence factors")[-1].split(":")[-1].strip()
                if sname not in data: data[sname] = {}
                data[sname]['vf_count'] = count
            else:
                if sname not in data: data[sname] = {}
                data[sname]['vf_count'] = "0"

    # 5. Parse Kraken2
    kraken_files = glob.glob("*_kraken2.report")
    for f in kraken_files:
        sname = get_sample_name(f, "_kraken2.report")
        with open(f) as fh:
            top_species = "Unknown"
            top_pct = 0
            for line in fh:
                cols = line.strip().split('\t')
                if len(cols) >= 6:
                    rank, pct, name = cols[3], float(cols[0]), cols[5].strip()
                    if rank == 'S' and pct > top_pct:
                        top_pct = pct
                        top_species = name
            if sname not in data: data[sname] = {}
            data[sname]['taxonomy'] = f"{top_species} ({top_pct}%)"

    # 6. Parse Panaroo
    if os.path.exists("summary_statistics.txt"):
        with open("summary_statistics.txt") as fh:
            for line in fh:
                if "Core genes" in line: pan_data['core'] = line.split()[-1]
                if "Total genes" in line: pan_data['total'] = line.split()[-1]

    # 7. Parse QUAST
    if os.path.exists("report.tsv"):
        with open("report.tsv") as fh:
            lines = [l.strip().split('\t') for l in fh.readlines()]
            if len(lines) > 0:
                header = lines[0]
                for i in range(1, len(header)):
                    qname = header[i]
                    target = qname
                    if target not in data:
                        for dname in data.keys():
                            if qname.startswith(dname) or dname.startswith(qname):
                                target = dname
                                break
                    if target not in data: data[target] = {}
                    for row in lines:
                        if row[0] == "N50": data[target]['n50'] = "{:,}".format(int(row[i]))
                        if row[0] == "Total length": data[target]['gsize'] = "{:.2f}".format(int(row[i])/1000000)
                        if row[0] == "# contigs (>= 0 bp)": data[target]['contigs'] = row[i]

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8"><title>Summary Report</title>
        <style>
            body {{ font-family: sans-serif; background: #f4f6f9; padding: 20px; }}
            .container {{ max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
            th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #eee; }}
            th {{ background: #2c3e50; color: white; }}
            .badge {{ padding: 2px 8px; border-radius: 10px; font-size: 0.8em; background: #eef2f7; }}
            .amr-cat {{ margin-bottom: 8px; border-left: 3px solid #3498db; padding-left: 10px; }}
            .amr-label {{ font-weight: bold; font-size: 0.8em; color: #2c3e50; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Analysis Summary <small style="font-size:0.4em; color:grey;">{now}</small></h1>
            {f'''<div style="display:flex; gap:20px; margin:20px 0;">
                <div style="flex:1; background:#f8fafc; padding:15px; border-radius:5px; text-align:center;"><b>{pan_data.get('core','-')}</b><br>Core Genes</div>
                <div style="flex:1; background:#f8fafc; padding:15px; border-radius:5px; text-align:center;"><b>{pan_data.get('total','-')}</b><br>Total Genes</div>
            </div>''' if pan_data else ''}
            <table>
                <thead><tr><th>Sample</th><th>Taxonomy</th><th>Reads</th><th>Yield(Mb)</th><th>Contigs</th><th>N50</th></tr></thead>
                <tbody>
    """
    for sname in sorted(data.keys()):
        html += f"""
                    <tr>
                        <td><b>{sname}</b></td>
                        <td><span class="badge">{data[sname].get('taxonomy','Unknown')}</span></td>
                        <td>{data[sname].get('reads','-')}</td>
                        <td>{data[sname].get('bases','-')}</td>
                        <td>{data[sname].get('contigs','-')}</td>
                        <td>{data[sname].get('n50','-')}</td>
                    </tr>
        """
    html += """</tbody></table>
            <h2>Resistome & MLST</h2>
            <table><thead><tr><th>Sample</th><th>MLST</th><th>AMR</th><th>Detailed Profile</th></tr></thead><tbody>"""
    for sname in sorted(data.keys()):
        st = f"{data[sname].get('scheme','-')}: {data[sname].get('st','-')}" if data[sname].get('st') else "-"
        amr_cat, amr_counts = data[sname].get('amr_categorized',{}), data[sname].get('amr_counts',{})
        res_html = ""
        for cat in ['AMR','VIRULENCE','METAL','BIOCIDE','ACID','HEAT']:
            if cat in amr_cat:
                res_html += f'<div class="amr-cat"><span class="amr-label">{cat} ({amr_counts[cat]})</span><div style="font-size:0.85em;">{amr_cat[cat]}</div></div>'
        html += f"""<tr><td><b>{sname}</b></td><td>{st}</td><td>{data[sname].get('amr_total','0')}</td><td>{res_html or 'None'}</td></tr>"""
    html += "</tbody></table></div></body></html>"
    with open('summary_report.html', 'w') as f: f.write(html)
    '''
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
