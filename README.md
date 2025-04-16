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
./scripts/install_fastqc.sh
```

### 1. Data Retrieval

Run the below command to download the RNAseq data in the background.
```bash
nohup ./scripts/download-data.sh -a PRJNA706136 -o ./data > data/download.log 2>&1 &
```
QC - FastQC

Mapping - Bowtie2

DEG - DESeq2

 - check metadata for experimental bias etc.

 - read original paper

 - functional analysis with *A. thaliana* reference, focus on well-annotated genes
