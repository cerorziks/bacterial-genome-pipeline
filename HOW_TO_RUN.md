# 🚀 HOW TO RUN THE PIPELINE

## Prerequisites
Before running, ensure you have the following installed:
- **Nextflow** (v22.10.0 or later)
- **Docker** or **Singularity**
- **Java** (v11 or v17)

## Quick Start in 3 Steps

### STEP 1: Prepare Your Data

Create a CSV file listing your samples. Example:

```bash
# Create your samplesheet
cat > my_samples.csv << 'EOF'
sample,read1,read2
sample1,/full/path/to/sample1_R1.fastq.gz,/full/path/to/sample1_R2.fastq.gz
sample2,/full/path/to/sample2_R1.fastq.gz,/full/path/to/sample2_R2.fastq.gz
EOF
```

**Important:** Use **absolute paths** (full paths starting with `/`) for your FASTQ files.

### STEP 2: Run the Pipeline

**Basic run (without phylogeny):**
```bash
# Navigate to the pipeline directory
cd bacterial-genome-pipeline

# Run the pipeline using Docker
nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  -profile docker
```

**Full analysis with phylogeny:**
```bash
nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  --reference /path/to/reference.fasta \
  -profile docker
```

### STEP 3: Check Results

```bash
# View the interactive report in your browser
open results/multiqc/multiqc_report.html

# Or check specific results
ls -lh results/amr/          # AMR genes
ls -lh results/virulence/    # Virulence factors
ls -lh results/phylogeny/    # Trees and SNP distances
```

## Advanced Usage

### Resource Management
The pipeline uses these resources by default:
- **Max CPUs**: 16 (adjust with `--max_cpus 8`)
- **Max Memory**: 128 GB (adjust with `--max_memory 32.GB`)
- **Max Time**: 240 hours

Example with limited resources:
```bash
nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  --max_cpus 4 \
  --max_memory 16.GB \
  -profile docker
```

### Resume Failed Runs
If the pipeline stops due to an error, you can resume from the last successful step using `-resume`:
```bash
nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  -profile docker \
  -resume
```

## What You'll Get

After completion, your `results/` folder contains:

```
results/
├── multiqc/
│   └── multiqc_report.html    ← Interactive summary report
├── amr/
│   ├── sample1_amr.tsv        ← AMR genes detected
│   └── sample1_amr_summary.txt
├── virulence/
│   ├── sample1_virulence.tsv  ← Virulence factors
│   └── sample1_virulence_summary.txt
├── assembly/
│   └── sample1/contigs.fasta  ← Assembled genomes
├── phylogeny/                  ← If reference provided
│   ├── core.aln.treefile      ← Phylogenetic tree
│   └── snp_distances.tsv      ← Pairwise SNP matrix
└── mlst/
    └── sample1_mlst.tsv       ← Sequence types (ST)
```

## Troubleshooting

### Docker Permissions
If you get a "Permission denied" error with Docker, ensure your user is in the `docker` group:
```bash
sudo usermod -aG docker $USER
# Then LOG OUT and LOG BACK IN
```

---

## Ready to Run! 🎉

For more detailed information on parameters and output structure, please refer to the main [README.md](README.md).
