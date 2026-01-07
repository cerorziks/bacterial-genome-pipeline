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
include { QUAST                   } from './modules/assembly'
include { PROKKA                  } from './modules/annotation'
include { AMRFINDERPLUS           } from './modules/amr_detection'
include { VFDB_BLAST              } from './modules/virulence'
include { MLST                    } from './modules/mlst'
include { SNIPPY                  } from './modules/phylogeny'
include { SNIPPY_CORE             } from './modules/phylogeny'
include { IQTREE                  } from './modules/phylogeny'
include { SNP_DISTS               } from './modules/phylogeny'
include { MULTIQC                 } from './modules/reporting'

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
    
    // 3. Quality Control - Trimmed reads
    FASTQC_TRIMMED(FASTP.out.reads)
    
    // 4. Genome assembly
    SPADES(FASTP.out.reads)
    
    // 5. Assembly quality assessment
    QUAST(SPADES.out.assembly.map { it[1] }.collect())
    
    // 6. Genome annotation
    PROKKA(SPADES.out.assembly)
    
    // 7. AMR detection
    AMRFINDERPLUS(PROKKA.out.faa)
    
    // 8. Virulence factor detection
    VFDB_BLAST(PROKKA.out.faa)
    
    // 9. MLST typing
    if (!params.skip_mlst) {
        MLST(SPADES.out.assembly)
    }
    
    // 10. Phylogenetic analysis (if reference provided)
    if (params.reference && !params.skip_phylogeny) {
        reference_ch = Channel.fromPath(params.reference)
        
        // Variant calling for each sample
        SNIPPY(SPADES.out.assembly, reference_ch)
        
        // Core SNP alignment
        SNIPPY_CORE(
            SNIPPY.out.snippy_dir.collect(),
            reference_ch
        )
        
        // Build phylogenetic tree
        IQTREE(SNIPPY_CORE.out.alignment)
        
        // Calculate SNP distances
        SNP_DISTS(SNIPPY_CORE.out.alignment)
    }
    
    // 11. MultiQC reporting
    multiqc_files = Channel.empty()
        .mix(FASTQC_RAW.out.zip.map{it[1]}.collect().ifEmpty([]))
        .mix(FASTQC_TRIMMED.out.zip.map{it[1]}.collect().ifEmpty([]))
        .mix(FASTP.out.json.map{it[1]}.collect().ifEmpty([]))
        .mix(QUAST.out.results.ifEmpty([]))
        .mix(PROKKA.out.txt.map{it[1]}.collect().ifEmpty([]))
    
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
