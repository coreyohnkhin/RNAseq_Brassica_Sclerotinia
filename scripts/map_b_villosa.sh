#!/bin/bash

mkdir -p data/bowtie2
mkdir -p data/bowtie2/N1896Pet_ControlBOL
mkdir -p data/bowtie2/N1896Pet_ControlSsc
mkdir -p data/bowtie2/N1896Pet_InfectedBOL
mkdir -p data/bowtie2/N1896Pet_InfectedSsc

# --local allows for softclipping (fastqc results show strong base (G/C) bias in first 10 bp)
# map control to BOL reference index
for sample in 75 78 79; do
    ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR138380${sample}_1.fastq -2 data/N1896Pet_Control/SRR138380${sample}_2.fastq -S data/bowtie2/N1896Pet_ControlBOL/${sample}BOL.sam
    rm data/N1896Pet_Control/SRR138380${sample}_1.fastq data/N1896Pet_Control/SRR138380${sample}_2.fastq
done

# map infected to BOL reference index
for sample in 72 73 74; do
    ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR138380${sample}_1.fastq -2 data/N1896Pet_Infected/SRR138380${sample}_2.fastq -S data/bowtie2/N1896Pet_InfectedBOL/${sample}BOL.sam
    rm data/N1896Pet_Infected/SRR138380${sample}_1.fastq data/N1896Pet_Infected/SRR138380${sample}_2.fastq
done

# map control to Ssc reference index
for sample in 75 78 79; do
    ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR138380${sample}_1.fastq -2 data/N1896Pet_Control/SRR138380${sample}_2.fastq -S data/bowtie2/N1896Pet_ControlSsc/${sample}Ssc.sam
    rm data/N1896Pet_Control/SRR138380${sample}_1.fastq data/N1896Pet_Control/SRR138380${sample}_2.fastq
done

# map infected to Ssc reference index
for sample in 72 73 74; do
    ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR138380${sample}_1.fastq -2 data/N1896Pet_Infected/SRR138380${sample}_2.fastq -S data/bowtie2/N1896Pet_InfectedSsc/${sample}Ssc.sam
    rm data/N1896Pet_Infected/SRR138380${sample}_1.fastq data/N1896Pet_Infected/SRR138380${sample}_2.fastq
done
