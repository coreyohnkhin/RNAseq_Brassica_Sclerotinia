#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Error: No directory pattern provided"
    echo "Usage: $0 pattern"
    echo "Example: $0 N1896*"
    exit 1
fi

PATTERN="$1"

DIRS=$(find "data/" -type d -name "$PATTERN" -not -path "*/\.*")

if [ -z "$DIRS" ]; then
    echo "No directories matching pattern '$PATTERN' found in data/"
    exit 1
fi

echo "Found directories matching pattern '$PATTERN':"
echo "$DIRS"
echo ""

for DIR in $DIRS; do
    BASE_DIR=$(basename "$DIR")

    OUTPUT_DIR="data/${BASE_DIR}_fastqc"
    mkdir -p "$OUTPUT_DIR"

    echo "Processing directory: $DIR"
    echo "Output directory: $OUTPUT_DIR"

    FASTQ_COUNT=$(find "$DIR" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | wc -l)

    if [ "$FASTQ_COUNT" -eq 0 ]; then
        echo "Warning: No fastq files found in $DIR"
        continue
    fi

    echo "Running FastQC on $FASTQ_COUNT files..."
    fastqc "$DIR"/*.fastq "$DIR"/*.fastq.gz -o "$OUTPUT_DIR" 2>/dev/null || fastqc "$DIR"/*.fastq.gz -o "$OUTPUT_DIR" 2>/dev/null || fastqc "$DIR"/*.fastq -o "$OUTPUT_DIR" 2>/dev/null

    echo "FastQC processing completed for $DIR"
    echo ""
done

echo "All FastQC processing completed!"
