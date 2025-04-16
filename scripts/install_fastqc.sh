#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
SOFTWARE_DIR="${BASE_DIR}/software"
FASTQC_VERSION="0.12.1"
FASTQC_URL="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip"
FASTQC_ZIP="${SOFTWARE_DIR}/fastqc_v${FASTQC_VERSION}.zip"

echo "Installing FastQC v${FASTQC_VERSION} to ${SOFTWARE_DIR}"

if ! command -v java &> /dev/null; then
    echo "Java is not installed. FastQC requires Java to run."
    echo "Please install Java before continuing (OpenJDK 11+ recommended)."
    exit 1
fi

mkdir -p "${SOFTWARE_DIR}"

echo "Downloading FastQC..."
if command -v wget &> /dev/null; then
    wget -q --show-progress -O "${FASTQC_ZIP}" "${FASTQC_URL}"
elif command -v curl &> /dev/null; then
    curl --progress-bar -o "${FASTQC_ZIP}" "${FASTQC_URL}"
else
    echo "Error: Neither wget nor curl is installed. Cannot download FastQC."
    exit 1
fi

if [ ! -f "${FASTQC_ZIP}" ]; then
    echo "Error: Failed to download FastQC."
    exit 1
fi

unzip -q -o "${FASTQC_ZIP}" -d "${SOFTWARE_DIR}"

chmod +x "${SOFTWARE_DIR}/FastQC/fastqc"

rm "${FASTQC_ZIP}"

echo "Testing FastQC installation..."
"${SOFTWARE_DIR}/FastQC/fastqc" --version

echo "Installation successful!"
echo "Full path to FastQC: ${SOFTWARE_DIR}/FastQC/fastqc"
