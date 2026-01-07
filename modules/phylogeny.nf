/*
========================================================================================
    Phylogeny Module
========================================================================================
*/

process SNIPPY {
    tag "$sample"
    publishDir "${params.outdir}/phylogeny/snippy", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    path reference
    
    output:
    tuple val(sample), path("${sample}"), emit: snippy_dir
    tuple val(sample), path("${sample}/${sample}.vcf"), emit: vcf
    
    script:
    """
    snippy \\
        --outdir ${sample} \\
        --ref ${reference} \\
        --R1 ${read1} \\
        --R2 ${read2} \\
        --prefix ${sample} \\
        --cpus $task.cpus \\
        --force
    """
}

process SNIPPY_CORE {
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    
    input:
    path(snippy_dirs)
    path(reference)
    
    output:
    path("core.aln"), emit: alignment
    path("core.tab"), emit: tab
    path("core.txt"), emit: txt
    path("core.vcf"), emit: vcf
    
    script:
    """
    snippy-core \\
        --ref ${reference} \\
        --prefix core \\
        ${snippy_dirs}
    """
}

process IQTREE {
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    
    input:
    path(alignment)
    
    output:
    path("core.aln.treefile"), emit: tree
    path("core.aln.iqtree"), emit: iqtree
    path("core.aln.log"), emit: log
    
    script:
    """
    iqtree \\
        -s ${alignment} \\
        -m GTR+G \\
        -bb 1000 \\
        -nt $task.cpus \\
        -redo
    """
}

process SNP_DISTS {
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    
    input:
    path(alignment)
    
    output:
    path("snp_distances.tsv"), emit: distances
    
    script:
    """
    snp-dists ${alignment} > snp_distances.tsv
    """
}
