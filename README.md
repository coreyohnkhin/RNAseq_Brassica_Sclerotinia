# *Brassica-Sclerotinia* interaction RNAseq

RNA-Seq + DEG analysis of *B. oleracea, B. villosa* with Sclerotinia infection. To conserve storage on the server, we will run teach species separately.

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
All reads for this data pass QC.
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

Map the reads to the reference genomes.

Gene count matrices - `featureCounts`

DEG - `DESeq2`

 - check metadata for experimental bias etc.

 - functional analysis with *A. thaliana* reference, focus on well-annotated genes
