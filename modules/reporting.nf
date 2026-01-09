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
                cols = lines[1].split('\\t')
                if sname not in data: data[sname] = {}
                data[sname]['reads'] = "{:,}".format(int(cols[3]))
                data[sname]['bases'] = "{:.1f}".format(int(cols[4])/1000000)

    # 2. Parse MLST for ST
    mlst_files = glob.glob("*_mlst.tsv")
    log(f"Found {len(mlst_files)} MLST files")
    for f in mlst_files:
        sname = get_sample_name(f, "_mlst.tsv")
        with open(f) as fh:
            line = fh.readline().strip()
            if line:
                cols = line.split('\\t')
                if sname not in data: data[sname] = {}
                data[sname]['st'] = cols[2] if len(cols) > 2 else "-"
                data[sname]['scheme'] = cols[1] if len(cols) > 1 else "-"

    # 3. Parse AMR for Gene Count and Detailed List (Categorized)
    amr_files = glob.glob("*_amr.tsv")
    log(f"Found {len(amr_files)} AMR files")
    for f in amr_files:
        sname = get_sample_name(f, "_amr.tsv")
        if sname not in data: data[sname] = {}
        
        counts = {'AMR': 0, 'VIRULENCE': 0, 'METAL': 0, 'BIOCIDE': 0, 'ACID': 0, 'HEAT': 0, 'OTHER': 0}
        categorized_genes = {k: [] for k in counts.keys()}
        
        with open(f) as fh:
            header = fh.readline()
            for l in fh:
                c = l.strip().split('\\t')
                if len(c) > 10:
                    gene = c[5]
                    element_type = c[8].upper()
                    subclass = c[10]
                    
                    cat = 'OTHER'
                    if 'AMR' in element_type: cat = 'AMR'
                    elif 'VIRULENCE' in element_type: cat = 'VIRULENCE'
                    elif 'METAL' in element_type: cat = 'METAL'
                    elif 'BIOCIDE' in element_type: cat = 'BIOCIDE'
                    elif 'ACID' in element_type: cat = 'ACID'
                    elif 'HEAT' in element_type: cat = 'HEAT'
                    
                    categorized_genes[cat].append(f"{gene} ({subclass})")
                    counts[cat] += 1
            
            data[sname]['amr_counts'] = counts
            data[sname]['amr_categorized'] = {k: ", ".join(sorted(set(v))) for k, v in categorized_genes.items() if v}
            data[sname]['amr_total'] = sum(counts.values())

    # 4. Parse Virulence
    vir_files = glob.glob("*_virulence_summary.txt")
    log(f"Found {len(vir_files)} virulence summary files")
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

    # 5. Parse Kraken2 for Taxonomy
    kraken_files = glob.glob("*_kraken2.report")
    log(f"Found {len(kraken_files)} Kraken2 report files")
    for f in kraken_files:
        sname = get_sample_name(f, "_kraken2.report")
        with open(f) as fh:
            top_species = "Unknown"
            top_pct = 0
            for line in fh:
                cols = line.strip().split('\\t')
                if len(cols) >= 6:
                    rank = cols[3]
                    pct = float(cols[0])
                    name = cols[5].strip()
                    if rank == 'S' and pct > top_pct:
                        top_pct = pct
                        top_species = name
            if sname not in data: data[sname] = {}
            data[sname]['taxonomy'] = f"{top_species} ({top_pct}%)"

    # 6. Parse Panaroo Stats
    if os.path.exists("summary_statistics.txt"):
        log("Found Panaroo summary statistics")
        with open("summary_statistics.txt") as fh:
            for line in fh:
                if "Core genes" in line: pan_data['core'] = line.split()[-1]
                if "Total genes" in line: pan_data['total'] = line.split()[-1]
    else:
        log("No Panaroo summary statistics found")

    # 7. Parse QUAST
    if os.path.exists("report.tsv"):
        log("Found QUAST report")
        with open("report.tsv") as fh:
            lines = [l.strip().split('\\t') for l in fh.readlines()]
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
    else:
        log("No QUAST report found")

    # Generate HTML
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    
    html = f\"\"\"
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Bacterial Genome Analysis Summary</title>
        <link rel="icon" href="data:image/svg+xml,<svg xmlns=%22http://www.w3.org/2000/svg%22 viewBox=%220 0 100 100%22><text y=%22.9em%22 font-size=%2290%22>🧬</text></svg>">
        <style>
            :root {{
                --primary: #2c3e50;
                --secondary: #34495e;
                --accent: #3498db;
                --success: #27ae60;
                --bg: #f4f7f6;
                --text: #333;
                --white: #ffffff;
            }}
            body {{
                font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
                background-color: var(--bg);
                color: var(--text);
                margin: 0;
                padding: 20px;
                line-height: 1.6;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background: var(--white);
                padding: 40px;
                box-shadow: 0 4px 6px rgba(0,0,0,0.1);
                border-radius: 8px;
            }}
            h1 {{
                color: var(--primary);
                border-bottom: 3px solid var(--accent);
                padding-bottom: 10px;
                margin-top: 0;
                display: flex;
                justify-content: space-between;
                align-items: center;
            }}
            .timestamp {{
                font-size: 0.4em;
                color: #777;
                font-weight: normal;
            }}
            h2 {{
                color: var(--secondary);
                background: #edf2f7;
                padding: 10px 15px;
                border-radius: 4px;
                margin-top: 30px;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                table-layout: auto;
            }}
            th, td {{
                padding: 12px 15px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }}
            th {{
                background-color: var(--primary);
                color: var(--white);
                font-weight: 600;
                text-transform: uppercase;
                font-size: 0.85em;
                letter-spacing: 1px;
            }}
            tr:hover {{
                background-color: #f9f9f9;
            }}
            .badge {{
                display: inline-block;
                padding: 2px 8px;
                border-radius: 12px;
                font-size: 0.85em;
                font-weight: 600;
            }}
            .badge-primary {{ background: #ebf8ff; color: #2b6cb0; }}
            .badge-success {{ background: #f0fff4; color: #276749; }}
            .badge-info {{ background: #e6fffa; color: #2c7a7b; }}
            .badge-warn {{ background: #fffaf0; color: #744210; }}
            .amr-list {{
                font-size: 0.85em;
                color: #555;
                font-style: italic;
                max-width: 400px;
                word-wrap: break-word;
            }}
            .resistome-cat {{
                margin-bottom: 5px;
            }}
            .resistome-name {{
                font-weight: bold;
                font-size: 0.8em;
                color: var(--primary);
                display: block;
            }}
            .resistome-genes {{
                font-size: 0.9em;
            }}
            .tech-note {{
                font-size: 0.8em;
                color: #718096;
                background: #fdf2f2;
                padding: 10px;
                border-left: 4px solid #feb2b2;
                margin-bottom: 20px;
                border-radius: 0 4px 4px 0;
            }}
            footer {{
                margin-top: 50px;
                text-align: center;
                font-size: 0.9em;
                color: #999;
                border-top: 1px solid #eee;
                padding-top: 20px;
            }}
            .summary-box {{
                display: flex;
                gap: 20px;
                margin-bottom: 20px;
            }}
            .metric {{
                flex: 1;
                background: #f8fafc;
                padding: 15px;
                border-radius: 8px;
                border: 1px solid #e2e8f0;
                text-align: center;
                transition: transform 0.2s;
            }}
            .metric:hover {{
                transform: translateY(-2px);
                box-shadow: 0 4px 12px rgba(0,0,0,0.05);
            }}
            .metric-val {{
                font-size: 1.5em;
                font-weight: bold;
                color: var(--accent);
                display: block;
            }}
            .metric-label {{
                font-size: 0.8em;
                color: #64748b;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }}
            .tree-section {{
                padding: 20px;
                background: #f8fafc;
                border: 1px dashed #cbd5e0;
                border-radius: 4px;
                text-align: center;
                color: #718096;
            }}
        </style>
    </head>
    <body>
        <div class=\"container\">
            <h1>
                Bacterial Genome Analysis Summary
                <span class=\"timestamp\">Generated on: {now}</span>
            </h1>

            <div class=\"tech-note\">
                <strong>Note:</strong> Taxonomic verification performed using the MDU-PHL BabyKraken database (10MB). Results are optimized for common public health pathogens; low percentages or \"Unknown\" may occur for rare species.
            </div>
    \"\"\"

    if pan_data:
        html += f\"\"\"
            <section>
                <h2>0. Pangenome Overview (Population Analysis)</h2>
                <div class=\"summary-box\">
                    <div class=\"metric\">
                        <span class=\"metric-val\">{pan_data.get('core', '-')}</span>
                        <span class=\"metric-label\">Core Genes (Shared by all)</span>
                    </div>
                    <div class=\"metric\">
                        <span class=\"metric-val\">{pan_data.get('total', '-')}</span>
                        <span class=\"metric-label\">Total Pangenome Genes</span>
                    </div>
                    <div class=\"metric\">
                        <span class=\"metric-val\">{len(data.keys())}</span>
                        <span class=\"metric-label\">Total Samples Analyzed</span>
                    </div>
                </div>
            </section>
        \"\"\"

    html += \"\"\"
            <section>
                <h2>1. Sequencing and Assembly Statistics</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>Taxonomy (Best Match)</th>
                            <th>Total Reads</th>
                            <th>Yield (Mb)</th>
                            <th>Contigs</th>
                            <th>N50 (bp)</th>
                            <th>Genome Size (Mb)</th>
                        </tr>
                    </thead>
                    <tbody>
    \"\"\"

    for sname in sorted(data.keys()):
        html += f\"\"\"
                        <tr>
                            <td><strong>{sname}</strong></td>
                            <td><span class=\"badge badge-info\">{data[sname].get('taxonomy', 'Unknown')}</span></td>
                            <td>{data[sname].get('reads', '-')}</td>
                            <td>{data[sname].get('bases', '-')}</td>
                            <td>{data[sname].get('contigs', '-')}</td>
                            <td>{data[sname].get('n50', '-')}</td>
                            <td>{data[sname].get('gsize', '-')}</td>
                        </tr>
        \"\"\"

    html += \"\"\"
                    </tbody>
                </table>
            </section>

            <section>
                <h2>2. MLST, AMR & Virulence Profiles</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>MLST / ST</th>
                            <th>Total Genes</th>
                            <th>Virulence Factors</th>
                            <th>Detailed Resistome Profile</th>
                        </tr>
                    </thead>
                    <tbody>
    \"\"\"

    for sname in sorted(data.keys()):
        st = f\"{data[sname].get('scheme', '-')}: {data[sname].get('st', '-')}\" if data[sname].get('st') else \"-\"
        amr_counts = data[sname].get('amr_counts', {})
        amr_cat = data[sname].get('amr_categorized', {})
        
        resistome_html = \"\"
        for cat in ['AMR', 'VIRULENCE', 'METAL', 'BIOCIDE', 'ACID', 'HEAT']:
            if cat in amr_cat:
                count = amr_counts.get(cat, 0)
                genes = amr_cat.get(cat, '')
                resistome_html += f'<div class=\"resistome-cat\"><span class=\"resistome-name\">{cat} ({count})</span><span class=\"resistome-genes\">{genes}</span></div>'
        
        if not resistome_html: resistome_html = \"None detected\"

        html += f\"\"\"
                        <tr>
                            <td><strong>{sname}</strong></td>
                            <td><span class=\"badge badge-primary\">{st}</span></td>
                            <td><span class=\"badge badge-success\">{data[sname].get('amr_total', '0')}</span></td>
                            <td>{data[sname].get('vf_count', '0')}</td>
                            <td class=\"amr-list\">{resistome_html}</td>
                        </tr>
        \"\"\"

    html += \"\"\"
                    </tbody>
                </table>
            </section>
    \"\"\"

    if os.path.exists(\"*.treefile\") or glob.glob(\"*.treefile\"):
        html += \"\"\"
            <section>
                <h2>3. Phylogenetic Analysis</h2>
                <div class=\"tree-section\">
                    <p>Phylogenetic tree generation completed. See results directory for NEWICK files (.treefile).</p>
                </div>
            </section>
        \"\"\"

    html += f\"\"\"
            <footer>
                Generated by Bacterial Genome Analysis Pipeline v1.1.0-polished
            </footer>
        </div>
    </body>
    </html>
    \"\"\"

    with open('summary_report.html', 'w') as f:
        f.write(html)
    
    log(\"Report generation completed successfully.\")
    \"\"\"
}

process MULTIQC {
    publishDir \"${params.outdir}/multiqc\", mode: 'copy'
    
    input:
    path(files)
    
    output:
    path(\"multiqc_report.html\"), emit: html
    path(\"*_data\"), emit: data
    
    script:
    \"\"\"
    multiqc \\\\
        --title \"${params.multiqc_title}\" \\\\
        --filename multiqc_report.html \\\\
        --force \\\\
        .
    \"\"\"
}
