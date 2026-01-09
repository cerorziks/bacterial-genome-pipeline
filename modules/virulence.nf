/*
========================================================================================
    Virulence Factor Detection Module
========================================================================================
*/

process DOWNLOAD_VFDB {
    publishDir "${params.outdir}/databases", mode: 'copy'
    storeDir "${params.db_cache}/vfdb"
    
    output:
    path "VFDB_setB_pro.fas", emit: db
    
    script:
    """
    # Check if VFDB database already exists
    if [ -f "VFDB_setB_pro.fas" ]; then
        echo "VFDB database already exists. Skipping download." >&2
        echo "Database location: \$(pwd)/VFDB_setB_pro.fas" >&2
        exit 0
    fi
    
    echo "[DB] Starting VFDB database download..." >&2
    if curl -L -o VFDB_setB_pro.fas.gz http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz; then
        echo "Downloaded from primary source." >&2
    else
        echo "Primary source failed, trying backup..." >&2
        curl -L -o VFDB_setB_pro.fas.gz https://github.com/arpcard/VFDB/raw/master/VFDB_setB_pro.fas.gz
    fi
    
    gunzip -f VFDB_setB_pro.fas.gz
    echo "[DB] VFDB download and extraction complete." >&2
    """
}

process VFDB_BLAST {
    tag "$sample"
    publishDir "${params.outdir}/virulence", mode: 'copy'
    
    input:
    tuple val(sample), path(protein_fasta)
    path vfdb_file
    
    output:
    tuple val(sample), path("${sample}_virulence.tsv"), emit: results
    tuple val(sample), path("${sample}_virulence_summary.txt"), emit: summary
    
    script:
    """
    # Create BLAST database
    makeblastdb -in ${vfdb_file} -dbtype prot -out vfdb_db
    
    # Run BLAST
    blastp \\
        -query ${protein_fasta} \\
        -db vfdb_db \\
        -out ${sample}_virulence.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -evalue 1e-5 \\
        -num_threads $task.cpus \\
        -max_target_seqs 1
    
    # Create summary
    if [ -s ${sample}_virulence.tsv ]; then
        echo "Virulence Factors Found for ${sample}:" > ${sample}_virulence_summary.txt
        echo "=======================================" >> ${sample}_virulence_summary.txt
        awk '\$3 >= 70 && \$4 >= 50' ${sample}_virulence.tsv | \\
            cut -f13 | sed 's/(.*//g' | sort -u >> ${sample}_virulence_summary.txt
        echo "" >> ${sample}_virulence_summary.txt
        echo "Total virulence factors (>70% identity, >50 aa): \$(awk '\$3 >= 70 && \$4 >= 50' ${sample}_virulence.tsv | wc -l)" >> ${sample}_virulence_summary.txt
    else
        echo "No virulence factors detected for ${sample}" > ${sample}_virulence_summary.txt
    fi
    """
}
