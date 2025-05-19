#!/bin/bash

mkdir -p data/featurecountmatrices

# GCM for Control BOL
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene_id \
    -a data/annotations/Brassica_oleracea.BOL.61.gtf \
    -o data/featurecountmatrices/N1896Pet_ControlBOL.txt \
    data/bowtie2/N1896Pet_ControlBOL/75BOL.bam \
    data/bowtie2/N1896Pet_ControlBOL/78BOL.bam \
    data/bowtie2/N1896Pet_ControlBOL/79BOL.bam

# GCM for Infected BOL
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene_id \
    -a data/annotations/Brassica_oleracea.BOL.61.gtf \
    -o data/featurecountmatrices/N1896Pet_InfectedBOL.txt \
    data/bowtie2/N1896Pet_InfectedBOL/72BOL.bam \
    data/bowtie2/N1896Pet_InfectedBOL/73BOL.bam \
    data/bowtie2/N1896Pet_InfectedBOL/74BOL.bam

# GCM for Control Ssc
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene_id \
    -a data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gtf \
    -o data/featurecountmatrices/N1896Pet_ControlSsc.txt \
    data/bowtie2/N1896Pet_ControlSsc/75Ssc.bam \
    data/bowtie2/N1896Pet_ControlSsc/78Ssc.bam \
    data/bowtie2/N1896Pet_ControlSsc/79Ssc.bam

# GCM for Infected Ssc
./software/subread-2.1.1-Linux-x86_64/bin/featureCounts \
    -p --countReadPairs -t exon -g gene_id \
    -a data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gtf \
    -o data/featurecountmatrices/N1896Pet_InfectedSsc.txt \
    data/bowtie2/N1896Pet_InfectedSsc/72Ssc.bam \
    data/bowtie2/N1896Pet_InfectedSsc/73Ssc.bam \
    data/bowtie2/N1896Pet_InfectedSsc/74Ssc.bam
