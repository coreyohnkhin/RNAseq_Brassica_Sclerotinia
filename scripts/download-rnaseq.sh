#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <SRA_ACCESSION>"
    exit 1
fi

SRA_ACCESSION=$1

mkdir -p data

software/sratoolkit.3.2.1-ubuntu64/bin/prefetch $SRA_ACCESSION

for sra_file in $(find . -name "SRR*" -type d); do
    software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files --threads 4 --outdir data $sra_file
done
