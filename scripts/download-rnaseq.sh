#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <SRA_ACCESSION>"
    echo "Example: $0 PRJNA706136"
    exit 1
fi

SRA_ACCESSION=$1
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
DATA_DIR="$PROJECT_ROOT/data"

mkdir -p "$DATA_DIR"

echo "$(date) - Starting download of $SRA_ACCESSION"

echo "$(date) - Running prefetch for $SRA_ACCESSION"
"$PROJECT_ROOT/software/sratoolkit.3.2.1-ubuntu64/bin/prefetch" "$SRA_ACCESSION"

echo "$(date) - Converting SRA files to FASTQ"
for sra_file in $(find "$PROJECT_ROOT" -name "SRR*" -type d); do
    echo "$(date) - Processing $sra_file"
    "$PROJECT_ROOT/software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump" --split-files --threads 4 --outdir "$DATA_DIR" "$sra_file"
done

echo "$(date) - Download and conversion complete"
