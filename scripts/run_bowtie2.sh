#!/bin/bash

set -e

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BOWTIE2_PATH="$PROJECT_ROOT/software/bowtie2-2.5.4-linux-x86_64"
BOWTIE2="$BOWTIE2_PATH/bowtie2"
BOWTIE2_BUILD="$BOWTIE2_PATH/bowtie2-build"
SAMTOOLS_PATH="$PROJECT_ROOT/software/samtools-1.21/bin"
SAMTOOLS="$SAMTOOLS_PATH/samtools"

usage() {
    echo "Usage: $0 <fastq_directory> <reference_genome_fasta>"
    echo "  <fastq_directory>: Directory containing paired-end FASTQ files (_1.fastq and _2.fastq)"
    echo "  <reference_genome_fasta>: Path to the reference genome FASTA file"
    exit 1
}

check_dependencies() {
    if [ ! -f "$BOWTIE2" ]; then
        exit 1
    fi

    if [ ! -f "$BOWTIE2_BUILD" ]; then
        exit 1
    fi

    if [ ! -f "$SAMTOOLS" ]; then
        exit 1
    fi
}

if [ "$#" -ne 2 ]; then
    usage
fi

FASTQ_DIR="$1"
REF_GENOME="$2"
INDEX_BASE="${REF_GENOME%.*}"
OUTPUT_DIR="${FASTQ_DIR}_bowtie2_results"
THREADS=4

if [ ! -d "$FASTQ_DIR" ]; then
    exit 1
fi

if [ ! -f "$REF_GENOME" ]; then
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

check_dependencies

if [ ! -f "${INDEX_BASE}.1.bt2" ]; then
    "$BOWTIE2_BUILD" "$REF_GENOME" "$INDEX_BASE"
else
    :
fi

SAMPLES=$(find "$FASTQ_DIR" -name "*_1.fastq" | sed 's/_1\.fastq$//')

if [ -z "$SAMPLES" ]; then
    exit 1
fi

for SAMPLE in $SAMPLES; do
    SAMPLE_NAME=$(basename "$SAMPLE")

    READ1="${SAMPLE}_1.fastq"
    READ2="${SAMPLE}_2.fastq"

    if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
        continue
    fi

    SAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.sam"
    BAM_FILE="$OUTPUT_DIR/${SAMPLE_NAME}.bam"
    SORTED_BAM="$OUTPUT_DIR/${SAMPLE_NAME}.sorted.bam"

    "$BOWTIE2" \
        --local \
        --threads "$THREADS" \
        -x "$INDEX_BASE" \
        -1 "$READ1" \
        -2 "$READ2" \
        -S "$SAM_FILE" \
        2> "$OUTPUT_DIR/${SAMPLE_NAME}.bowtie2.log"

    "$SAMTOOLS" view -bS "$SAM_FILE" > "$BAM_FILE"
    "$SAMTOOLS" sort "$BAM_FILE" -o "$SORTED_BAM"
    "$SAMTOOLS" index "$SORTED_BAM"
    "$SAMTOOLS" flagstat "$SORTED_BAM" > "$OUTPUT_DIR/${SAMPLE_NAME}.flagstat.txt"
done
