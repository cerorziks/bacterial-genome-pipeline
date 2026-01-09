/*
========================================================================================
    Pangenomics Module (Panaroo)
========================================================================================
*/

process PANAROO {
    label 'process_high'
    publishDir "${params.outdir}/pangenomics", mode: 'copy'

    input:
    path gffs

    output:
    path "results/*", emit: results
    path "results/gene_presence_absence.csv", emit: csv
    path "results/pan_genome_reference.fa", emit: reference
    path "results/summary_statistics.txt", emit: stats

    script:
    """
    num_files=\$(echo ${gffs} | wc -w)
    if [ "\$num_files" -lt 2 ]; then
        echo "Pangenomics requires at least 2 samples. Skipping Panaroo."
        mkdir results
        echo "Core genes: N/A" > results/summary_statistics.txt
        echo "Total genes: N/A" >> results/summary_statistics.txt
        touch results/gene_presence_absence.csv
        touch results/pan_genome_reference.fa
    else
        echo "[Pangenomics] Running Panaroo on $num_files input files..."
        panaroo \
            -i ${gffs} \
            -o results \
            --clean-mode strict \
            -t ${task.cpus} \
            --remove-invalid-genes
    fi
    """
}
