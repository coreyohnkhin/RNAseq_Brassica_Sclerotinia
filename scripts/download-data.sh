#!/bin/bash

SRA_TOOLKIT_PATH="../sratoolkit.3.2.1-ubuntu64/bin"
OUTPUT_DIR="../data"
ACCESSION="PRJNA706136"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

$SRA_TOOLKIT_PATH/vdb-config --root -o cache

echo "Downloading metadata for $ACCESSION..."
$SRA_TOOLKIT_PATH/prefetch "$ACCESSION" --output-directory .

echo "Getting run accessions for $ACCESSION..."
SRA_RUNS=$($SRA_TOOLKIT_PATH/sra-stat --meta --quick "$ACCESSION" | grep "Accession:" | grep "SRR" | awk '{print $2}')

echo "Downloading individual runs..."
for run in $SRA_RUNS; do
  echo "Processing $run..."

  $SRA_TOOLKIT_PATH/prefetch "$run" --output-directory .

  echo "Converting $run to FASTQ..."
  $SRA_TOOLKIT_PATH/fasterq-dump "$run" --outdir . --split-files

  echo "Compressing FASTQ files for $run..."
  gzip "${run}"*.fastq
done

echo "Download complete. Data saved to $OUTPUT_DIR"
