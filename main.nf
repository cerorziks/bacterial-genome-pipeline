#!/usr/bin/env nextflow

/*
========================================================================================
    Bacterial Genome Analysis Pipeline
========================================================================================
    Author: Nextflow Pipeline
    Version: 1.0.0
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FASTQC as FASTQC_RAW    } from './modules/quality_control'
include { FASTP                   } from './modules/quality_control'
include { FASTQC as FASTQC_TRIMMED} from './modules/quality_control'
include { SPADES                  } from './modules/assembly'
include { SHOVILL                 } from './modules/assembly'
include { QUAST                   } from './modules/assembly'
include { PROKKA                  } from './modules/annotation'
include { AMRFINDERPLUS           } from './modules/amr_detection'
include { DOWNLOAD_AMR_DB         } from './modules/amr_detection'
include { VFDB_BLAST              } from './modules/virulence'
include { DOWNLOAD_VFDB           } from './modules/virulence'
include { MLST                    } from './modules/mlst'
include { SNIPPY                  } from './modules/phylogeny'
include { SNIPPY_CORE             } from './modules/phylogeny'
include { IQTREE                  } from './modules/phylogeny'
include { SNP_DISTS               } from './modules/phylogeny'
include { SEQKIT_STATS            } from './modules/utils'
include { SEQKIT_STATS_ASSEMBLY   } from './modules/utils'
include { PANAROO                 } from './modules/pangenomics'
include { KRAKEN2                 } from './modules/taxonomy'
include { DOWNLOAD_KRAKEN2_DB     } from './modules/taxonomy'
include { MULTIQC                  } from './modules/reporting'
include { COMPILED_HTML_REPORT     } from './modules/reporting'

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================
    Bacterial Genome Analysis Pipeline v1.0.0
    =========================================
    
    Usage:
        nextflow run main.nf --input samplesheet.csv --outdir results [options]
    
    Required Arguments:
        --input         Path to samplesheet CSV file with columns: sample,read1,read2
        --outdir        Path to output directory
    
    Optional Arguments:
        --reference     Path to reference genome for phylogeny (FASTA format)
        --vfdb          Path to VFDB database (default: auto-download)
        --skip_phylogeny  Skip phylogenetic analysis (default: false)
        --skip_mlst     Skip MLST typing (default: false)
        
    Profile Options:
        -profile docker     Run with Docker containers
        -profile singularity Run with Singularity containers
        -profile local      Run locally without containers
    
    Example:
        nextflow run main.nf \\
            --input samples.csv \\
            --outdir results \\
            --reference ref.fasta \\
            -profile docker
    
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input) {
    log.error "ERROR: --input parameter is required"
    helpMessage()
    exit 1
}

if (!params.outdir) {
    log.error "ERROR: --outdir parameter is required"
    helpMessage()
    exit 1
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Create input channel from samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            tuple(row.sample, file(row.read1), file(row.read2))
        }
        .set { reads_ch }
    
    // 1. Quality Control - Raw reads
    FASTQC_RAW(reads_ch)
    
    // 2. Read trimming and quality filtering
    FASTP(reads_ch)

    // 2.0 Taxonomy verification
    if (params.kraken_db) {
        kraken_db = Channel.fromPath(params.kraken_db)
    } else {
        DOWNLOAD_KRAKEN2_DB()
        kraken_db = DOWNLOAD_KRAKEN2_DB.out.db
    }
    KRAKEN2(FASTP.out.reads, kraken_db)
    
    // 2.1 SeqKit stats on trimmed reads
    SEQKIT_STATS(FASTP.out.reads)
    
    // 3. Quality Control - Trimmed reads
    FASTQC_TRIMMED(FASTP.out.reads)
    
    // 4. Genome assembly
    if (params.assembler == 'shovill') {
        SHOVILL(FASTP.out.reads)
        assembly_ch = SHOVILL.out.assembly
    } else {
        SPADES(FASTP.out.reads)
        assembly_ch = SPADES.out.assembly
    }
    
    // 4.1 SeqKit stats on assembly
    SEQKIT_STATS_ASSEMBLY(assembly_ch)
    
    // 5. Assembly quality assessment
    QUAST(
        assembly_ch.map { it[1] }.collect(),
        assembly_ch.map { it[0] }.collect()
    )
    
    // 6. Genome annotation
    PROKKA(assembly_ch)
    
    // 7. AMR detection
    if (params.amr_db) {
        amr_db_ch = Channel.fromPath(params.amr_db)
    } else {
        DOWNLOAD_AMR_DB()
        amr_db_ch = DOWNLOAD_AMR_DB.out.db
    }
    AMRFINDERPLUS(PROKKA.out.faa, amr_db_ch)
    
    // 8. Virulence factor detection
    if (params.vfdb) {
        vfdb_ch = Channel.fromPath(params.vfdb)
    } else {
        DOWNLOAD_VFDB()
        vfdb_ch = DOWNLOAD_VFDB.out.db
    }
    VFDB_BLAST(PROKKA.out.faa, vfdb_ch)
    
    // 9. MLST typing
    if (!params.skip_mlst) {
        MLST(assembly_ch)
    }
    
    // 10. Phylogenetic analysis (if reference provided)
    if (params.reference && !params.skip_phylogeny) {
        reference_file = file(params.reference)
        
        // Variant calling for each sample (using reads for better accuracy)
        SNIPPY(FASTP.out.reads, reference_file)
        
        // Core SNP alignment
        SNIPPY_CORE(
            SNIPPY.out.snippy_dir.map { it[1] }.collect(),
            reference_file
        )
        
        // Build phylogenetic tree
        IQTREE(SNIPPY_CORE.out.alignment)
        
        // Calculate SNP distances
        SNP_DISTS(SNIPPY_CORE.out.alignment)
    }
    
    // 10. Pangenomics (if enabled)
    pangenome_stats_ch = Channel.value([])
    if (!params.skip_pangenomics) {
        PANAROO(PROKKA.out.gff.map{it[1]}.collect())
        pangenome_stats_ch = PANAROO.out.stats
    }

    // 11. Final Consolidated HTML Report (Standalone)
    tree_report_ch = (params.reference && !params.skip_phylogeny) ? IQTREE.out.tree : Channel.value([])
    
    COMPILED_HTML_REPORT(
        SEQKIT_STATS.out.stats.map{it[1]}.collect().ifEmpty([]),
        SEQKIT_STATS_ASSEMBLY.out.stats.map{it[1]}.collect().ifEmpty([]),
        QUAST.out.tsv.ifEmpty([]),
        AMRFINDERPLUS.out.results.map{it[1]}.collect().ifEmpty([]),
        MLST.out.results.map{it[1]}.collect().ifEmpty([]),
        VFDB_BLAST.out.summary.map{it[1]}.collect().ifEmpty([]),
        KRAKEN2.out.report.map{it[1]}.collect().ifEmpty([]),
        pangenome_stats_ch.ifEmpty([]),
        tree_report_ch
    )

    // 12. MultiQC reporting (Preprocessing QC only)
    multiqc_files = Channel.empty()
        .mix(FASTQC_RAW.out.zip.map{it[1]}.collect().ifEmpty([]))
        .mix(FASTQC_TRIMMED.out.zip.map{it[1]}.collect().ifEmpty([]))
        .mix(FASTP.out.json.map{it[1]}.collect().ifEmpty([]))
        .mix(SEQKIT_STATS.out.stats.map{it[1]}.collect().ifEmpty([]))
    
    MULTIQC(multiqc_files.collect())
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "========================================="
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Duration: $workflow.duration"
    log.info "Output directory: ${params.outdir}"
    log.info "========================================="
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
