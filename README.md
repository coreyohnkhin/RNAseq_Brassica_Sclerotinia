# *Brassica-Sclerotinia* interaction RNAseq

RNA-Seq + DEG analysis of *B. oleracea, B. villosa* with Sclerotinia infection.

## Usage
```bash
git clone https://github.com/coreyohnkhin/RNAseq_Brassica_Sclerotinia
cd RNAseq_Brassica_Sclerotinia
```

### Install software
```bash
./scripts/install_sratoolkit.sh
./scripts/install_bowtie2.sh
```

`FastQC` and `DESeq2` should already be installed.

### 1. Data Retrieval

Run the below command to download the RNAseq data in the background.
```bash
nohup ./scripts/download_rnaseq.sh PRJNA706136.txt > data/download_rnaseq.log 2>&1 &
```
### 2. QC
Run the below command to perform QC with `FastQC` on the .fastq files for each run.
```bash
nohup ./scripts/fastqc > data/fastqc.log 2>&1 &
```

Mapping - `Bowtie2`

DEG - `DESeq2`

 - check metadata for experimental bias etc.

 - functional analysis with *A. thaliana* reference, focus on well-annotated genes
