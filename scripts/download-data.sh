#!/bin/bash

DEFAULT_SRA_TOOLKIT_PATH="./software/sratoolkit.3.2.1-ubuntu64/bin"
DEFAULT_OUTPUT_DIR="./data"
DEFAULT_ACCESSION="PRJNA706136"

function print_usage {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  -a, --accession ACCESSION   SRA accession number (default: $DEFAULT_ACCESSION)"
  echo "  -o, --output DIR            Output directory (default: $DEFAULT_OUTPUT_DIR)"
  echo "  -t, --toolkit PATH          Path to SRA toolkit bin directory (default: $DEFAULT_SRA_TOOLKIT_PATH)"
  echo "  -h, --help                  Show this help message"
  exit 1
}

SRA_TOOLKIT_PATH=$DEFAULT_SRA_TOOLKIT_PATH
OUTPUT_DIR=$DEFAULT_OUTPUT_DIR
ACCESSION=$DEFAULT_ACCESSION

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -a|--accession)
      ACCESSION="$2"
      shift
      shift
      ;;
    -o|--output)
      OUTPUT_DIR="$2"
      shift
      shift
      ;;
    -t|--toolkit)
      SRA_TOOLKIT_PATH="$2"
      shift
      shift
      ;;
    -h|--help)
      print_usage
      ;;
    *)
      echo "Unknown option: $1"
      print_usage
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

if [[ ! "$SRA_TOOLKIT_PATH" = /* ]]; then
  SRA_TOOLKIT_PATH="$BASE_DIR/$SRA_TOOLKIT_PATH"
fi

if [[ ! "$OUTPUT_DIR" = /* ]]; then
  OUTPUT_DIR="$BASE_DIR/$OUTPUT_DIR"
fi

echo "Using SRA toolkit at: $SRA_TOOLKIT_PATH"
echo "Saving data to: $OUTPUT_DIR"
echo "Processing accession: $ACCESSION"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

mkdir -p "$HOME/.ncbi"
echo "Configuring SRA toolkit..."
$SRA_TOOLKIT_PATH/vdb-config --restore-defaults

echo "Downloading metadata for $ACCESSION..."
$SRA_TOOLKIT_PATH/prefetch "$ACCESSION" --output-directory .

echo "Getting run accessions for $ACCESSION..."
SRA_RUNS=$($SRA_TOOLKIT_PATH/sra-stat --meta --quick "$ACCESSION" | grep "Accession:" | grep "SRR" | awk '{print $2}')

if [ -z "$SRA_RUNS" ]; then
  echo "Could not get run accessions via sra-stat. Checking for downloaded SRA files..."

  SRA_RUNS=$(find . -maxdepth 1 -type d -name "SRR*" | sed 's/\.\///')

  if [ -z "$SRA_RUNS" ]; then
    echo "Using hardcoded list of known runs for $ACCESSION"
    SRA_RUNS="SRR13838068 SRR13838069 SRR13838070 SRR13838071 SRR13838072 SRR13838073 SRR13838074 SRR13838075 SRR13838076 SRR13838077 SRR13838078 SRR13838079"
  fi
fi

echo "Found the following runs: $SRA_RUNS"

echo "Downloading and processing individual runs..."
for run in $SRA_RUNS; do
  echo "-------------------------------------------------------------------------"
  echo "Processing $run..."

  if [ ! -d "$run" ]; then
    echo "Downloading $run..."
    $SRA_TOOLKIT_PATH/prefetch "$run" --output-directory .
  else
    echo "Found existing directory for $run, skipping download"
  fi

  if ls "${run}"*.fastq.gz 1> /dev/null 2>&1; then
    echo "FASTQ files for $run already exist, skipping conversion"
    continue
  fi

  TEMP_DIR=$(mktemp -d -p .)
  echo "Created temporary directory: $TEMP_DIR"

  echo "Checking available disk space..."
  AVAIL_SPACE=$(df -k . | awk 'NR==2 {print $4}')
  echo "Available space: $(($AVAIL_SPACE / 1024)) MB"

  echo "Converting $run to FASTQ..."
  $SRA_TOOLKIT_PATH/fasterq-dump "$run" \
    --outdir . \
    --temp "$TEMP_DIR" \
    --split-files \
    --buffer-size 100MB \
    --threads 4 || {
      echo "Error converting $run to FASTQ. Trying alternate method..."
      $SRA_TOOLKIT_PATH/fastq-dump "$run" --outdir . --split-files
    }

  rm -rf "$TEMP_DIR"

  if ls "${run}"*.fastq 1> /dev/null 2>&1; then
    echo "Compressing FASTQ files for $run..."
    gzip "${run}"*.fastq
    echo "$run processed successfully"
  else
    echo "WARNING: No FASTQ files were created for $run"
  fi
done

echo "Download complete. Data saved to $OUTPUT_DIR"
