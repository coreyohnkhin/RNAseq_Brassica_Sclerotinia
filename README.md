# *Brassica-Sclerotinia* interaction RNAseq

RNA-Seq + DEG analysis of [*B. oleracea, B. villosa* with Sclerotinia infection](https://pubmed.ncbi.nlm.nih.gov/36966424/). To conserve storage on the server, we will run each species separately. A lot of the commands will be run with the `nohup` command. It's tedious to keep a live terminal for these long processes (esp. the RNA-seq download and read mapping). Everything should be run in sequence, however, so use `jobs` to check if each process is finished before continuing. Use `watch tail [logfile.log]` to view the live progress and check for any errors.

## Usage
```bash
git clone https://github.com/coreyohnkhin/RNAseq_Brassica_Sclerotinia
cd RNAseq_Brassica_Sclerotinia
```

### Install software
```bash
./scripts/install_sratoolkit.sh
./scripts/install_bowtie2.sh
./scripts/install_samtools.sh
./scripts/install_subread.sh
./scripts/install_gffread.sh
```

`FastQC` and `DESeq2` should already be installed.

## *Brassica villosa* (BRA1896)
### 1. Data Retrieval

Run the below command to download the RNAseq data for *B. villosa* in the background.
```bash
nohup ./scripts/download_rnaseq.sh data/N1896*.txt > data/download_rnaseq.log 2>&1 &
```
### 2. QC
Run the below command to perform QC with `FastQC` on the .fastq files.
```bash
nohup ./scripts/run_fastqc.sh N1896* > data/fastqc.log 2>&1 &
```
All reads for this data pass QC, but requires aligner clipping (strong base (G/C) bias in first 10 bp)
### 3. Mapping
Download the reference genomes for *B. oleracea* and *Sclerotinia sclerotiorum*.
```bash
nohup ./scripts/download_genomes.sh > data/genomes.log 2>&1 &
```
Build `bowtie2` indices.
```bash
mkdir -p data/indices

# Brassica oleracea reference genome index
nohup ./software/bowtie2-2.5.4-linux-x86_64/bowtie2-build \
data/reference_genomes/Brassica_oleracea.BOL.dna.toplevel.fa \
data/indices/brassica_oleracea \
> data/brassica_oleracea_index.log 2>&1 &

# Sclerotinia sclerotiorum reference genome index
nohup ./software/bowtie2-2.5.4-linux-x86_64/bowtie2-build \
data/reference_genomes/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa \
data/indices/sclerotinia_sclerotiorum \
> data/sclerotinia_sclerotiorum_index.log 2>&1 &
```

Map the reads to the reference indices. This removes the downloaded RNA-seq data after each replicate is mapped and converts SAM to BAM.
```bash
nohup ./scripts/map_b_villosa.sh > data/b_villosa_mapping.log 2>&1 &
```

### 4. Generate gene count matrices

Download the reference genome (*B. oleracea* and *S. sclerotiorum*) annotations (GFF3 format).
```bash
nohup ./scripts/download_annotations.sh > data/download_annotation.log 2>&1 &
```
Convert annotations to GTF

```bash
./software/gffread/gffread \
    data/annotations/Brassica_oleracea.BOL.61.gff3 \
    -T -o data/annotations/Brassica_oleracea.BOL.61.gtf
./software/gffread/gffread \
    data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gff3 \
    -T -o data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gtf
```
Generate gene count matrices.
```bash
nohup ./scripts/b_villosa_featurecounts.sh > data/featurecounts.log 2>&1 &
```

### 5. Differentially expressed genes analysis
Run the below **locally** (will need the download the gene count matrices off of the server) and then check the output.
```bash
Rscript deg_b_villosaBOL.r
Rscript deg_b_villosaSsc.r
```
Output:

 - **Significant genes count** - prints in console n = # DEGs

 - `results/results_[BOL/Ssc]infected_vs_control.csv` - log₂FC, pval, padj for each DEG

 - `figures/heatmap_top50_[BOL/Ssc]_.pdf` - clustered expression map

### 6. Gene set enrichment analysis

Also run locally.

```bash
Rscript scripts/gsea_b_villosa.r results/results_BOL_infected_vs_control.csv results/results_Ssc_infected_vs_control.csv
```
## *Brassica villosa* (BRA1896)
### 1. Data retrieval
```bash
nohup ./scripts/download_rnaseq.sh data/N1909*.txt \
      > data/download_rnaseq_N1909.log 2>&1 &
```
### 2. QC
```bash
nohup ./scripts/run_fastqc.sh N1909* \
      > data/fastqc_N1909.log 2>&1 &
```

### 3. Mapping
```bash
# use the same reference genomes / indices built earlier
nohup ./scripts/map_b_oleracea.sh > data/b_oleracea_mapping.log 2>&1 &
```

### 4. Gene counts
```bash
nohup ./scripts/b_oleracea_featurecounts.sh > data/featurecounts_N1909.log 2>&1 &
```

### 5. DEG analysis
```bash
Rscript deg_b_oleraceaBOL.r
Rscript deg_b_oleraceaSsc.r
```

### 6. GSEA
```bash
Rscript gsea_b_oleracea.r results/results_BOL_infected_vs_control.csv \
                           results/results_Ssc_infected_vs_control.csv
```
