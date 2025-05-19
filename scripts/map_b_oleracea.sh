#!/bin/bash

set -euo pipefail
mkdir -p data/bowtie2 \
         data/bowtie2/N1909Pet_{Control,Infected}{BOL,Ssc}

for sample in 69 70 71; do
  ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -p 8 \
      -x data/indices/brassica_oleracea \
      -1 data/N1909Pet_Control/SRR138380${sample}_1.fastq \
      -2 data/N1909Pet_Control/SRR138380${sample}_2.fastq \
  | software/samtools-1.21/samtools view -bS -o \
      data/bowtie2/N1909Pet_ControlBOL/${sample}BOL.bam

  ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -p 8 \
      -x data/indices/sclerotinia_sclerotiorum \
      -1 data/N1909Pet_Control/SRR138380${sample}_1.fastq \
      -2 data/N1909Pet_Control/SRR138380${sample}_2.fastq \
  | software/samtools-1.21/samtools view -bS -o \
      data/bowtie2/N1909Pet_ControlSsc/${sample}Ssc.bam
  rm data/N1909Pet_Control/SRR138380${sample}_*.fastq
done

for sample in 68 76 77; do
  ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -p 8 \
      -x data/indices/brassica_oleracea \
      -1 data/N1909Pet_Infected/SRR138380${sample}_1.fastq \
      -2 data/N1909Pet_Infected/SRR138380${sample}_2.fastq \
  | software/samtools-1.21/samtools view -bS -o \
      data/bowtie2/N1909Pet_InfectedBOL/${sample}BOL.bam

  ./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -p 8 \
      -x data/indices/sclerotinia_sclerotiorum \
      -1 data/N1909Pet_Infected/SRR138380${sample}_1.fastq \
      -2 data/N1909Pet_Infected/SRR138380${sample}_2.fastq \
  | software/samtools-1.21/samtools view -bS -o \
      data/bowtie2/N1909Pet_InfectedSsc/${sample}Ssc.bam
  rm data/N1909Pet_Infected/SRR138380${sample}_*.fastq
done
