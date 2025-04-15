#!/bin/bash

# Script to download SRA data using sratoolkit
# Can be run from the RNAseq_Brassica_Sclerotinia directory

# Default values
DEFAULT_SRA_TOOLKIT_PATH="./software/sratoolkit.3.2.1-ubuntu64/bin"
DEFAULT_OUTPUT_DIR="./data"
DEFAULT_ACCESSION="PRJNA706136"

# Help function
function print_usage {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  -a, --accession ACCESSION   SRA accession number (default: $DEFAULT_ACCESSION)"
  echo "  -o, --output DIR            Output directory (default: $DEFAULT_OUTPUT_DIR)"
  echo "  -t, --toolkit PATH          Path to SRA toolkit bin directory (default: $DEFAULT_SRA_TOOLKIT_PATH)"
  echo "  -h, --help                  Show this help message"
  exit 1
}

# Parse command line arguments
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

# Ensure the script works correctly from any location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Resolve paths if they're relative
if [[ ! "$SRA_TOOLKIT_PATH" = /* ]]; then
  SRA_TOOLKIT_PATH="$BASE_DIR/$SRA_TOOLKIT_PATH"
fi

if [[ ! "$OUTPUT_DIR" = /* ]]; then
  OUTPUT_DIR="$BASE_DIR/$OUTPUT_DIR"
fi

echo "Using SRA toolkit at: $SRA_TOOLKIT_PATH"
echo "Saving data to: $OUTPUT_DIR"
echo "Processing accession: $ACCESSION"

# Create output directory and navigate to it
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Configure the SRA toolkit
$SRA_TOOLKIT_PATH/vdb-config --root -o cache

echo "Downloading metadata for $ACCESSION..."
$SRA_TOOLKIT_PATH/prefetch "$ACCESSION" --output-directory .

echo "Getting run accessions for $ACCESSION..."
SRA_RUNS=$($SRA_TOOLKIT_PATH/sra-stat --meta --quick "$ACCESSION" | grep "Accession:" | grep "SRR" | awk '{print $2}')

echo "Found the following runs: $SRA_RUNS"

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
