#!/bin/bash

mkdir -p data/N1896Pet_Control_bowtie2
mkdir -p data/N1896Pet_Control_bowtie2/BOL
mkdir -p data/N1896Pet_Control_bowtie2/Ssc
mkdir -p data/N1896Pet_Infected_bowtie2
mkdir -p data/N1896Pet_Infected_bowtie2/BOL
mkdir -p data/N1896Pet_Infected_bowtie2/SSc

# map control to BOL reference index
# --local allows for softclipping (fastqc results show strong base (G/C) bias in first 10 bp)
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR13838075_1.fastq -2 data/N1896Pet_Control/SRR13838075_2.fastq -S data/N1896Pet_Control_bowtie2/BOL/75BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR13838078_1.fastq -2 data/N1896Pet_Control/SRR13838078_2.fastq -S data/N1896Pet_Control_bowtie2/BOL/78BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Control/SRR13838079_1.fastq -2 data/N1896Pet_Control/SRR13838079_2.fastq -S data/N1896Pet_Control_bowtie2/BOL/79BOL.sam

# map infected to BOL reference index
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR13838072_1.fastq -2 data/N1896Pet_Infected/SRR13838072_2.fastq -S data/N1896Pet_Infected_bowtie2/BOL/72BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR13838073_1.fastq -2 data/N1896Pet_Infected/SRR13838073_2.fastq -S data/N1896Pet_Infected_bowtie2/BOL/73BOL.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/brassica_oleracea -1 data/N1896Pet_Infected/SRR13838074_1.fastq -2 data/N1896Pet_Infected/SRR13838074_2.fastq -S data/N1896Pet_Infected_bowtie2/BOL/74BOL.sam

# map control to Ssc reference index
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR13838075_1.fastq -2 data/N1896Pet_Control/SRR13838075_2.fastq -S data/N1896Pet_Control_bowtie2/Ssc/75Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR13838078_1.fastq -2 data/N1896Pet_Control/SRR13838078_2.fastq -S data/N1896Pet_Control_bowtie2/Ssc/78Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Control/SRR13838079_1.fastq -2 data/N1896Pet_Control/SRR13838079_2.fastq -S data/N1896Pet_Control_bowtie2/Ssc/79Ssc.sam

# map infected to Ssc reference index
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR13838072_1.fastq -2 data/N1896Pet_Infected/SRR13838072_2.fastq -S data/N1896Pet_Infected_bowtie2/Ssc/72Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR13838073_1.fastq -2 data/N1896Pet_Infected/SRR13838073_2.fastq -S data/N1896Pet_Infected_bowtie2/Ssc/73Ssc.sam
./software/bowtie2-2.5.4-linux-x86_64/bowtie2 --local -x data/indices/sclerotinia_sclerotiorum -1 data/N1896Pet_Infected/SRR13838074_1.fastq -2 data/N1896Pet_Infected/SRR13838074_2.fastq -S data/N1896Pet_Infected_bowtie2/Ssc/74Ssc.sam
