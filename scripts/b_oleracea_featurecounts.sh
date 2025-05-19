#!/bin/bash

set -euo pipefail
mkdir -p data/featurecountmatrices

FC=./software/subread-2.1.1-Linux-x86_64/bin/featureCounts
GTF_BO=data/annotations/Brassica_oleracea.BOL.61.gtf
GTF_SSC=data/annotations/Sclerotinia_sclerotiorum.ASM14694v1.61.gtf

$FC -p --countReadPairs -t exon -g gene_id -a $GTF_BO \
    -o data/featurecountmatrices/N1909Pet_ControlBOL.txt \
    data/bowtie2/N1909Pet_ControlBOL/{69,70,71}BOL.bam

$FC -p --countReadPairs -t exon -g gene_id -a $GTF_BO \
    -o data/featurecountmatrices/N1909Pet_InfectedBOL.txt \
    data/bowtie2/N1909Pet_InfectedBOL/{68,76,77}BOL.bam

$FC -p --countReadPairs -t exon -g gene_id -a $GTF_SSC \
    -o data/featurecountmatrices/N1909Pet_ControlSsc.txt \
    data/bowtie2/N1909Pet_ControlSsc/{69,70,71}Ssc.bam

$FC -p --countReadPairs -t exon -g gene_id -a $GTF_SSC \
    -o data/featurecountmatrices/N1909Pet_InfectedSsc.txt \
    data/bowtie2/N1909Pet_InfectedSsc/{68,76,77}Ssc.bam
