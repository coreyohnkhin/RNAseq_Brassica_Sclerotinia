#!/bin/bash

SRA_TOOLKIT="software/.3.2.1-ubuntu64/bin"

ACCESSION=""
OUTPUT_DIR="results"
TEMP_DIR="results/tmp"

while getopts "a:" opt; do
  case $opt in
    a)
      ACCESSION=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "Usage: ./scripts/download-rnaseq.sh -a <accession>" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "Usage: ./scripts/download-rnaseq.sh -a <accession>" >&2
      exit 1
      ;;
  esac
done

if [ -z "$ACCESSION" ]; then
  echo "Error: No accession provided."
  echo "Usage: ./scripts/download-rnaseq.sh -a <accession>"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

echo "Starting download for accession: $ACCESSION"
echo "Output will be stored in: $OUTPUT_DIR"

echo "Retrieving run information for $ACCESSION..."
"$SRA_TOOLKIT/efetch" -db sra -format runinfo -id "$ACCESSION" > "$OUTPUT_DIR/${ACCESSION}_runinfo.csv"

if [ ! -s "$OUTPUT_DIR/${ACCESSION}_runinfo.csv" ]; then
  echo "Error: Failed to retrieve run information for $ACCESSION."
  exit 1
fi

SRRS=$(tail -n +2 "$OUTPUT_DIR/${ACCESSION}_runinfo.csv" | cut -d "," -f 1)

if [ -z "$SRRS" ]; then
  echo "Error: No SRR numbers found for accession $ACCESSION."
  exit 1
fi

SRR_COUNT=$(echo "$SRRS" | wc -l)
echo "Found $SRR_COUNT run(s) for accession $ACCESSION."

COUNTER=1
for SRR in $SRRS; do
  echo "[$COUNTER/$SRR_COUNT] Processing $SRR..."

  echo "Downloading $SRR..."
  "$SRA_TOOLKIT/prefetch" "$SRR" -O "$TEMP_DIR"

  echo "Converting $SRR to FASTQ format..."
  "$SRA_TOOLKIT/fasterq-dump" "$TEMP_DIR/$SRR/$SRR.sra" \
    --outdir "$OUTPUT_DIR" \
    --temp "$TEMP_DIR" \
    --threads 4 \
    --split-files

  if [ $? -ne 0 ]; then
    echo "Error: Failed to convert $SRR to FASTQ format."
  else
    echo "Successfully processed $SRR."
    echo "Compressing FASTQ files..."
    find "$OUTPUT_DIR" -name "${SRR}*.fastq" -exec gzip -f {} \;
  fi

  COUNTER=$((COUNTER+1))
done

echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "Download and processing complete for accession $ACCESSION."
