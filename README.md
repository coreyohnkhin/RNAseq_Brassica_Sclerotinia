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
./scripts/install_subread.sh
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

Map the reads to the reference indices.
```bash
nohup ./scripts/map_b_villosa.sh > data/b_villosa_mapping.log 2>&1 &
```

### 4. Generate gene count matrices

Download the reference genome (*B. oleracea* and *S. sclerotiorum*) annotations (GFF3 format).
```bash
nohup ./scripts/download_annotations.sh > data/download_annotation.log 2>&1 &
```

`featureCounts` analysis.
```bash
# WIP
```

### DEG - `DESeq2`
WIP

 - check metadata for experimental bias etc.

 - functional analysis with *A. thaliana* reference, focus on well-annotated genes
