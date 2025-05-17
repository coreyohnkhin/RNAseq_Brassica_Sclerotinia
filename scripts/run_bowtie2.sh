#!/bin/bash

set -e

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BOWTIE2_PATH="$PROJECT_ROOT/software/bowtie2-2.5.4-linux-x86_64"
BOWTIE2="$BOWTIE2_PATH/bowtie2"
BOWTIE2_BUILD="$BOWTIE2_PATH/bowtie2-build"
SAMTOOLS_PATH="$PROJECT_ROOT/software/samtools-1.21"
SAMTOOLS="$SAMTOOLS_PATH/samtools"
THREADS=4

usage() {
    echo "Usage: $0 <fastq_directory> <reference_genome_fasta> [output_directory] [threads]"
    echo "  <fastq_directory>: Directory containing paired-end FASTQ files (_1.fastq and _2.fastq)"
    echo "  <reference_genome_fasta>: Path to the reference genome FASTA file"
    echo "  [output_directory]: Optional. Output directory (default: <fastq_directory>_bowtie2_results)"
    echo "  [threads]: Optional. Number of threads to use (default: 4)"
    exit 1
}

check_dependencies() {
    if [ ! -f "$BOWTIE2" ]; then
        echo "ERROR: Bowtie2 not found at $BOWTIE2"
        echo "Please run ./scripts/install_bowtie2.sh first"
        exit 1
    fi

    if [ ! -f "$BOWTIE2_BUILD" ]; then
        echo "ERROR: bowtie2-build not found at $BOWTIE2_BUILD"
        echo "Please run ./scripts/install_bowtie2.sh first"
        exit 1
    fi

    if [ ! -f "$SAMTOOLS" ]; then
        echo "ERROR: samtools not found at $SAMTOOLS"
        echo "Please run ./scripts/install_samtools.sh first"
        exit 1
    fi
}

if [ "$#" -lt 2 ]; then
    usage
fi

FASTQ_DIR="$1"
REF_GENOME="$2"
OUTPUT_DIR="${3:-${FASTQ_DIR}_bowtie2_results}"
[ "$4" ] && THREADS="$4"

INDEX_BASE="${REF_GENOME%.*}"

if [ ! -d "$FASTQ_DIR" ]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

if [ ! -f "$REF_GENOME" ]; then
    echo "ERROR: Reference genome file not found: $REF_GENOME"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
echo "Output will be stored in: $OUTPUT_DIR"

check_dependencies

if [ ! -f "${INDEX_BASE}.1.bt2" ]; then
    echo "Building Bowtie2 index for $REF_GENOME..."
    "$BOWTIE2_BUILD" "$REF_GENOME" "$INDEX_BASE" || { echo "ERROR: Failed to build Bowtie2 index"; exit 1; }
else
    echo "Using existing Bowtie2 index for $REF_GENOME"
fi

SAMPLES=$(find "$FASTQ_DIR" -name "*_1.fastq" | sed 's/_1\.fastq$//')

if [ -z "$SAMPLES" ]; then
    echo "ERROR: No paired-end FASTQ files found in $FASTQ_DIR"
    echo "Expected files should end with _1.fastq and _2.fastq"
    exit 1
fi

echo "Found $(echo "$SAMPLES" | wc -l) samples to process using $THREADS threads"

for SAMPLE in $SAMPLES; do
    SAMPLE_NAME=$(basename "$SAMPLE")
    echo "Processing sample: $SAMPLE_NAME"

    READ1="${SAMPLE}_1.fastq"
    READ2="${SAMPLE}_2.fastq"

    if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
        echo "WARNING: Missing one of the paired files for $SAMPLE_NAME, skipping"
        continue
    fi

    SAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.sam"
    BAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.bam"
    SORTED_BAM="$OUTPUT_DIR/${SAMPLE_NAME}.sorted.bam"

    echo "Mapping reads to reference genome..."
    "$BOWTIE2" \
        --local \
        --threads "$THREADS" \
        -x "$INDEX_BASE" \
        -1 "$READ1" \
        -2 "$READ2" \
        -S "$SAM_FILE" \
        2> "$OUTPUT_DIR/${SAMPLE_NAME}.bowtie2.log"

    echo "Converting SAM to BAM..."
    "$SAMTOOLS" view -bS "$SAM_FILE" > "$BAM_FILE"

    echo "Sorting BAM file..."
    "$SAMTOOLS" sort "$BAM_FILE" -o "$SORTED_BAM"

    echo "Indexing sorted BAM file..."
    "$SAMTOOLS" index "$SORTED_BAM"

    echo "Generating alignment statistics..."
    "$SAMTOOLS" flagstat "$SORTED_BAM" > "$OUTPUT_DIR/${SAMPLE_NAME}.flagstat.txt"

    echo "Completed processing $SAMPLE_NAME"
done

echo "All samples processed successfully"
