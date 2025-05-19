#!/bin/bash

mkdir -p data/bowtie2
mkdir -p data/bowtie2/N1896Pet_ControlBOL
mkdir -p data/bowtie2/N1896Pet_ControlSsc
mkdir -p data/bowtie2/N1896Pet_InfectedBOL
mkdir -p data/bowtie2/N1896Pet_InfectedSsc

# map control to BOL reference index
# --local allows for softclipping (fastqc results show strong base (G/C) bias in first 10 bp)
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR13838075_1.fastq -2 data/N1896Pet_Control/SRR13838075_2.fastq -S data/bowtie2/N1896Pet_ControlBOL/75BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR13838078_1.fastq -2 data/N1896Pet_Control/SRR13838078_2.fastq -S data/bowtie2/N1896Pet_ControlBOL/78BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR13838079_1.fastq -2 data/N1896Pet_Control/SRR13838079_2.fastq -S data/bowtie2/N1896Pet_ControlBOL/79BOL.sam

# map infected to BOL reference index
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR13838072_1.fastq -2 data/N1896Pet_Infected/SRR13838072_2.fastq -S data/bowtie2/N1896Pet_InfectedBOL/72BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR13838073_1.fastq -2 data/N1896Pet_Infected/SRR13838073_2.fastq -S data/bowtie2/N1896Pet_InfectedBOL/73BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR13838074_1.fastq -2 data/N1896Pet_Infected/SRR13838074_2.fastq -S data/bowtie2/N1896Pet_InfectedBOL/74BOL.sam

# map control to Ssc reference index
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR13838075_1.fastq -2 data/N1896Pet_Control/SRR13838075_2.fastq -S data/bowtie2/N1896Pet_ControlSsc/75Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR13838078_1.fastq -2 data/N1896Pet_Control/SRR13838078_2.fastq -S data/bowtie2/N1896Pet_ControlSsc/78Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR13838079_1.fastq -2 data/N1896Pet_Control/SRR13838079_2.fastq -S data/bowtie2/N1896Pet_ControlSsc/79Ssc.sam

# map infected to Ssc reference index
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR13838072_1.fastq -2 data/N1896Pet_Infected/SRR13838072_2.fastq -S data/bowtie2/N1896Pet_InfectedSsc/72Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR13838073_1.fastq -2 data/N1896Pet_Infected/SRR13838073_2.fastq -S data/bowtie2/N1896Pet_InfectedSsc/73Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR13838074_1.fastq -2 data/N1896Pet_Infected/SRR13838074_2.fastq -S data/bowtie2/N1896Pet_InfectedSsc/74Ssc.sam
