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
Map the reads to the reference genomes.
```bash
mkdir -p mapping_logs

# Map B. villosa control samples to B. oleracea reference
nohup ./scripts/run_bowtie2.sh data/N1896Pet_Control \
     data/reference_genomes/Brassica_oleracea.BOL.dna.toplevel.fa \
     data/N1896Pet_Control_BOL_mapped 8 \
     > mapping_logs/N1896Pet_Control.BOL.bowtie2.log 2>&1 &

# Map B. villosa control samples to S. sclerotiorum reference
nohup ./scripts/run_bowtie2.sh data/N1896Pet_Control \
     data/reference_genomes/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa \
     data/N1896Pet_Control_Ssc_mapped 8 \
     > mapping_logs/N1896Pet_Control.Ssc.bowtie2.log 2>&1 &

# Map B. villosa infected samples to B. oleracea reference
nohup ./scripts/run_bowtie2.sh data/N1896Pet_Infected \
     data/reference_genomes/Brassica_oleracea.BOL.dna.toplevel.fa \
     data/N1896Pet_Infected_BOL_mapped 8 \
     > mapping_logs/N1896Pet_Infected.BOL.bowtie2.log 2>&1 &

# Map B. villosa infected samples to S. sclerotiorum reference
nohup ./scripts/run_bowtie2.sh data/N1896Pet_Infected \
     data/reference_genomes/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa \
     data/N1896Pet_Infected_Ssc_mapped 8 \
     > mapping_logs/N1896Pet_Infected.Ssc.bowtie2.log 2>&1 &

```

DEG - `DESeq2`

 - check metadata for experimental bias etc.

 - functional analysis with *A. thaliana* reference, focus on well-annotated genes
