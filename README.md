# Bacterial Genome Analysis Pipeline

A comprehensive Nextflow pipeline for analyzing bacterial genomes from Illumina paired-end reads, performing quality control, assembly, AMR detection, phylogenetic analysis, and virulence factor identification.

## Features

- **Quality Control**: FastQC and fastp for read quality assessment and trimming
- **Statistics**: **SeqKit** for rapid read and assembly metrics
- **Assembly**: de novo assembly using **SPAdes** (default) or **Shovill** (faster)
- **Assembly QC**: QUAST quality metrics for assembly assessment
- **Annotation**: Prokka genome annotation
- **AMR Detection**: AMRFinderPlus with automated database management and persistent storage
- **Virulence**: BLAST against VFDB for virulence factor detection
- **MLST**: Multi-locus sequence typing
- **Phylogeny**: Core genome SNP-based phylogenetic tree construction
- **Reporting**: **Consolidated HTML report** summarizing results from all analysis stages.

## Pipeline Workflow

The pipeline processes bacterial genome data through the following stages:

```mermaid
graph TD
    A[Paired-End FASTQ Files] --> B[FastQC - Raw Reads]
    A --> C[fastp Trimming & QC]
    
    C --> D[FastQC - Trimmed Reads]
    C --> S1[SeqKit Read Stats]
    C --> E{Assembly Selection}
    
    E -- Default --> E1[SPAdes Assembly]
    E -- "--assembler shovill" --> E2[Shovill Assembly]
    
    E1 --> F[Assembly Results]
    E2 --> F
    
    F --> S2[SeqKit Assembly Stats]
    F --> G[QUAST Assembly QC]
    F --> H[Prokka Annotation]
    F --> I[MLST Typing]
    
    H --> J[AMRFinderPlus]
    H --> K[VFDB BLAST Virulence]
    
    C -.Reference Genome.-> L[Snippy Variant Calling]
    L --> M[snippy-core Alignment]
    M --> N[IQ-TREE Phylogeny]
    M --> P[SNP Distance Matrix]
    
    B --> MQC[MultiQC Dashboard]
    D --> MQC
    C --> MQC
    S1 --> MQC
    
    BS[Consolidated Report Generator]
    
    S1 --> BS
    S2 --> BS
    G --> BS
    J --> BS
    I --> BS
    K --> BS
    N --> BS
    
    BS --> Q[Final HTML Summary]
    
    style Q fill:#c8e6c9,stroke:#388e3c,stroke-width:2px
    style N fill:#fff9c4,stroke:#f57c00,stroke-width:2px
    style J fill:#ffccbc,stroke:#d84315,stroke-width:2px
    style K fill:#f8bbd0,stroke:#c2185b,stroke-width:2px
    style BS fill:#b2ebf2,stroke:#00acc1,stroke-width:2px
    style MQC fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style Q fill:#c8e6c9,stroke:#388e3c,stroke-width:2px
    style N fill:#fff9c4,stroke:#f57c00,stroke-width:2px
    style J fill:#ffccbc,stroke:#d84315,stroke-width:2px
    style K fill:#f8bbd0,stroke:#c2185b,stroke-width:2px
    style BS fill:#b2ebf2,stroke:#00acc1,stroke-width:2px
```

**Legend:**
- 🔵 **Blue** - Input data
- 🟢 **Green** - Final outputs and reports
- 🟡 **Yellow** - Phylogenetic analysis
- 🟠 **Orange** - AMR detection
- 🔴 **Pink** - Virulence detection
- **Dashed lines** - Optional analysis (requires reference genome)

## Requirements

### Software
- **Nextflow** >= 22.10.0 (Tested with v25.10.2)
- **Java** 11 or 17 (Tested with OpenJDK 17)
- **Docker** (Recommended) or **Singularity**
- **Git** (for cloning the repository)

### Hardware Recommendations
- **Minimum**: 8 CPUs, 32 GB RAM, 100 GB Free Disk Space
- **Recommended**: 16 CPUs, 128 GB RAM, 500 GB Free Disk Space (for large datasets)
- *Note: SPAdes assembly is memory-intensive. For large genomes or high depth, ensure at least 32GB RAM is available.*

## Environment Setup

### 1. Install Java and Docker
Ensure your system has Java and Docker installed. On Ubuntu/Debian:
```bash
sudo apt update
sudo apt install openjdk-17-jdk docker.io
# Add your user to the docker group
sudo usermod -aG docker $USER
# Log out and back in for changes to take effect
```

### 2. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
sudo chmod +x /usr/local/bin/nextflow
```

## Quick Start

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### 2. Prepare Input Samplesheet

Create a CSV file with your samples:

```csv
sample,read1,read2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz
```

### 3. Run the Pipeline

**Basic usage (without phylogeny):**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  -profile docker
```

**Full analysis with phylogeny:**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --reference reference_genome.fasta \
  -profile docker
```

## Parameters

### Required
| Parameter | Description |
|-----------|-------------|
| `--input` | Path to CSV samplesheet with sample information |
| `--outdir` | Path to output directory |

### Optional
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reference` | null | Reference genome for phylogeny (FASTA) |
| `--vfdb` | null | Path to VFDB database (auto-downloaded if not provided) |
| `--amr_db` | './amrfinder_db' | Path to store/load AMRFinder database |
| `--skip_phylogeny` | false | Skip phylogenetic analysis |
| `--skip_mlst` | false | Skip MLST typing |
| `--assembler` | 'spades' | Assembler to use ('spades' or 'shovill') |
| `--fastp_qualified_quality_phred` | 20 | Minimum quality score for fastp filtering |
| `--fastp_min_length` | 50 | Minimum read length after trimming |
| `--spades_kmers` | auto | K-mer sizes for SPAdes assembly |
| `--multiqc_title` | 'Bacterial Genome Analysis Report' | Title for MultiQC report |

### Resource Limits
| Parameter | Default |
|-----------|---------|
| `--max_cpus` | 16 |
| `--max_memory` | 128.GB |
| `--max_time` | 240.h |

## Output Structure

```
results/
├── fastqc/              # FastQC reports (raw and trimmed)
├── fastp/               # fastp trimming reports
├── seqkit/              # Read and assembly statistics
├── assembly/            # SPAdes/Shovill assemblies
├── quast/               # Assembly quality metrics
├── annotation/          # Prokka annotations
├── amr/                 # AMR resistance genes
├── virulence/           # Virulence factors
├── mlst/                # MLST typing results
├── phylogeny/           # Phylogenetic tree and SNP distances
│   ├── snippy/          # Per-sample variant calls
│   ├── core.aln         # Core genome alignment
│   ├── core.aln.treefile # Phylogenetic tree (Newick format)
│   └── snp_distances.tsv # Pairwise SNP distances
├── multiqc/             # MultiQC dashboard (all technical reports)
├── final_report/        # Consolidated standalone summary report (HTML)
└── pipeline_info/       # Pipeline execution reports
```

## Profiles

The pipeline includes several execution profiles:

- **`docker`**: Run with Docker containers (recommended)
- **`singularity`**: Run with Singularity containers
- **`local`**: Run without containers (requires all tools installed)
- **`slurm`**: Submit jobs to SLURM cluster
- **`test`**: Small-scale testing with reduced resources

Example:
```bash
nextflow run main.nf --input samples.csv --outdir results -profile docker
```

## Understanding Results

### AMR Detection
- **Location**: `results/amr/`
- Each sample has a `*_amr.tsv` file with detailed AMR gene predictions
- `*_amr_summary.txt` provides a quick overview of detected resistance genes

### Virulence Factors
- **Location**: `results/virulence/`
- `*_virulence.tsv` contains BLAST hits against VFDB
- `*_virulence_summary.txt` lists virulence factors with >70% identity

### Phylogeny
- **Location**: `results/phylogeny/`
- `core.aln.treefile` is the phylogenetic tree (view with FigTree, iTOL, etc.)
- `snp_distances.tsv` shows pairwise SNP distances between samples

### MultiQC Report
- **Location**: `results/multiqc/multiqc_report.html`
- Open this HTML file in a web browser for an interactive dashboard
- Aggregates QC metrics from all pipeline steps

## Troubleshooting

### Memory Issues
If SPAdes runs out of memory:
```bash
nextflow run main.nf --input samples.csv --outdir results --max_memory 64.GB
```

### Container Issues
If containers fail to pull, try:
```bash
# For Docker
docker pull staphb/spades:3.15.5

# For Singularity
singularity pull docker://staphb/spades:3.15.5
```

### Pipeline Resumption
If the pipeline fails, you can resume from where it stopped:
```bash
nextflow run main.nf --input samples.csv --outdir results -resume
```

## Citation

If you use this pipeline, please cite the tools:

- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology ([GitHub](https://github.com/nextflow-io/nextflow))
- **FastQC**: Andrews S. (2010) Babraham Bioinformatics ([GitHub](https://github.com/s-andrews/FastQC))
- **fastp**: Chen et al. (2018) Bioinformatics ([GitHub](https://github.com/OpenGene/fastp))
- **SPAdes**: Bankevich et al. (2012) Journal of Computational Biology ([GitHub](https://github.com/ablab/spades))
- **Prokka**: Seemann (2014) Bioinformatics ([GitHub](https://github.com/tseemann/prokka))
- **AMRFinderPlus**: Feldgarden et al. (2021) Scientific Reports ([GitHub](https://github.com/ncbi/amr))
- **VFDB**: Chen et al. (2005) Nucleic Acids Research ([Home](http://www.mgc.ac.cn/VFs/))
- **mlst**: Seemann T. ([GitHub: tseemann/mlst](https://github.com/tseemann/mlst))
- **Snippy**: Seemann T. ([GitHub: tseemann/snippy](https://github.com/tseemann/snippy))
- **IQ-TREE**: Nguyen et al. (2015) Molecular Biology and Evolution ([GitHub](https://github.com/iqtree/iqtree2))
- **MultiQC**: Ewels et al. (2016) Bioinformatics ([GitHub](https://github.com/MultiQC/MultiQC))

## License

This pipeline is distributed under the MIT License.

## Support

For issues and questions, please open an issue on the repository.
