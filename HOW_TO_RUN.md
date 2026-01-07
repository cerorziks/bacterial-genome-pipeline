# 🚀 HOW TO RUN THE PIPELINE

## Installation Complete! ✅

Your system is ready:
- ✅ Nextflow 25.10.2 installed at `/usr/local/bin/nextflow`
- ✅ Docker installed and running
- ✅ Pipeline validated and working

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

**Important:** Use **absolute paths** (full paths starting with `/`)

### STEP 2: Run the Pipeline

**Basic run (without phylogeny):**
```bash
cd /home/genomic/google_antigravity/AMR

# Note: Use sudo for Docker
sudo nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  -profile docker
```

**Full analysis with phylogeny:**
```bash
sudo nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  --reference /path/to/reference.fasta \
  -profile docker
```

### STEP 3: Check Results

```bash
# View the report in your browser
firefox results/multiqc/multiqc_report.html

# Or check specific results
ls -lh results/amr/          # AMR genes
ls -lh results/virulence/    # Virulence factors
ls -lh results/phylogeny/    # Trees and SNP distances
```

## Example with Real Data

Let's say your FASTQ files are in `/data/genomes/`:

```bash
# 1. Create samplesheet
cat > bacteria_samples.csv << 'EOF'
sample,read1,read2
ecoli_01,/data/genomes/ecoli_01_R1.fastq.gz,/data/genomes/ecoli_01_R2.fastq.gz
ecoli_02,/data/genomes/ecoli_02_R1.fastq.gz,/data/genomes/ecoli_02_R2.fastq.gz
salmonella,/data/genomes/salm_R1.fastq.gz,/data/genomes/salm_R2.fastq.gz
EOF

# 2. Run the pipeline
cd /home/genomic/google_antigravity/AMR
sudo nextflow run main.nf \
  --input bacteria_samples.csv \
  --outdir results_bacteria \
  -profile docker

# 3. Wait for completion (check progress in terminal)
# The pipeline shows: [█████████] 100% 15 of 15 ✔

# 4. View results
firefox results_bacteria/multiqc/multiqc_report.html
```

## Important Notes

### Docker Permissions
You need to use `sudo` with Docker commands, or:
```bash
# Add yourself to docker group (one-time setup)
sudo usermod -aG docker $USER
# Then LOG OUT and LOG BACK IN
# After that, you can run without sudo
```

### Resource Management
The pipeline uses these resources by default:
- **Max CPUs**: 16 (adjust with `--max_cpus 8`)
- **Max Memory**: 128 GB (adjust with `--max_memory 32.GB`)
- **Max Time**: 240 hours

Example with limited resources:
```bash
sudo nextflow run main.nf \
  --input my_samples.csv \
  --outdir results \
  --max_cpus 4 \
  --max_memory 16.GB \
  -profile docker
```

### Resume Failed Runs
If the pipeline stops, resume with `-resume`:
```bash
sudo nextflow run main.nf \
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
│   └── multiqc_report.html    ← Open this in browser!
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
│   └── snp_distances.tsv      ← SNP matrix
└── mlst/
    └── sample1_mlst.tsv       ← Sequence types
```

## Test the Installation

Run the test script:
```bash
cd /home/genomic/google_antigravity/AMR
./test_installation.sh
```

## Need Help?

- **Full documentation**: See `README.md`
- **Detailed guide**: See `QUICKSTART.md`
- **Troubleshooting**: Check `.nextflow.log`
- **GitHub**: https://github.com/cerorziks/bacterial-genome-pipeline

---

## Ready to Run! 🎉

Just prepare your samplesheet and execute:

```bash
cd /home/genomic/google_antigravity/AMR
sudo nextflow run main.nf \
  --input YOUR_SAMPLES.csv \
  --outdir results \
  -profile docker
```

**That's it!** The pipeline will handle everything automatically. ✨
