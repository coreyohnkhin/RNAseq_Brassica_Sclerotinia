#!/bin/bash

mkdir -p data/featurecountmatrices

# GCM for Control BOL
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene id \
    -a data/annotations/Brassica_oleracea.BOL.61.gff3 \
    -o data/featurecountmatrices/N1896Pet_ControlBOL.txt \
    data/bowtie2/N1896Pet_ControlBOL/75BOL.sam data/bowtie2/N1896Pet_ControlBOL/78BOL.sam data/bowtie2/N1896Pet_ControlBOL/79BOL.sam

# GCM for Infected BOL
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene id \
    -a data/annotations/Brassica_oleracea.BOL.61.gff3 \
    -o data/featurecountmatrices/N1896Pet_InfectedBOL.txt \
    data/bowtie2/N1896Pet_InfectedBOL/72BOL.sam data/bowtie2/N1896Pet_InfectedBOL/73BOL.sam data/bowtie2/N1896Pet_InfectedBOL/74BOL.sam

# GCM for Control Ssc
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene id \
    -a data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gff3 \
    -o data/featurecountmatrices/N1896Pet_ControlSsc.txt \
    data/bowtie2/N1896Pet_ControlSsc/75Ssc.sam data/bowtie2/N1896Pet_ControlSsc/78Ssc.sam data/bowtie2/N1896Pet_ControlSsc/79Ssc.sam

# GCM for Infected Ssc
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene id \
    -a data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gff3 \
    -o data/featurecountmatrices/N1896Pet_InfectedSsc.txt \
    data/bowtie2/N1896Pet_InfectedSsc/72Ssc.sam data/bowtie2/N1896Pet_InfectedSsc/73Ssc.sam data/bowtie2/N1896Pet_InfectedSsc/74Ssc.sam
