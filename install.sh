#!/bin/bash
# LowBioPipe Installation Script
# Creates the conda environment, makes the scripts executable and downloads the NCBI taxdump

set -euo pipefail

cd "$(dirname "$0")"

ENV_NAME="lowbiopipe"
TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

fail() {
    echo -e "${RED}ERROR: $1${NC}"
    exit 1
}

download() {
    if command -v curl &> /dev/null; then
        curl -fsSL -o "$2" "$1"
    elif command -v wget &> /dev/null; then
        wget -q -O "$2" "$1"
    else
        return 1
    fi
}

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}  LowBioPipe Installation${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""

# Find conda or mamba
echo -e "${GREEN}[1/4]${NC} Checking for conda/mamba..."
if command -v mamba &> /dev/null; then
    CONDA=mamba
elif command -v conda &> /dev/null; then
    CONDA=conda
else
    fail "conda or mamba is required. Install Miniforge: https://conda-forge.org/download/"
fi
echo -e "  Using ${CONDA} ${GREEN}OK${NC}"

# Create or update the environment (Python dependencies, Recentrifuge and Nextflow)
echo ""
echo -e "${GREEN}[2/4]${NC} Setting up the '${ENV_NAME}' conda environment from environment.yml..."
if "$CONDA" env list --json | grep -q "[/\\\\]envs[/\\\\]${ENV_NAME}\""; then
    "$CONDA" env update -n "$ENV_NAME" -f environment.yml --prune
else
    "$CONDA" env create -n "$ENV_NAME" -f environment.yml
fi
NXF_VERSION=$("$CONDA" run -n "$ENV_NAME" nextflow -version 2>/dev/null | grep -oE 'version [0-9.]+' | awk '{print $2}' || true)
[ -n "$NXF_VERSION" ] || fail "Nextflow not found in the '${ENV_NAME}' environment"
echo -e "  Environment ready (Nextflow ${NXF_VERSION}) ${GREEN}OK${NC}"

# Make scripts executable
echo ""
echo -e "${GREEN}[3/4]${NC} Making scripts executable..."
chmod +x bin/*.py
echo -e "  Scripts are now executable ${GREEN}OK${NC}"

# Download NCBI taxonomy database, verified against NCBI's checksum
echo ""
echo -e "${GREEN}[4/4]${NC} Checking NCBI taxonomy database..."
if [ ! -f "data/taxdump/nodes.dmp" ]; then
    echo -e "${YELLOW}  Taxonomy database not found. Downloading (~70 MB)...${NC}"
    TMP_DIR=$(mktemp -d)
    trap 'rm -rf "$TMP_DIR"' EXIT
    if download "$TAXDUMP_URL" "$TMP_DIR/taxdump.tar.gz" && download "$TAXDUMP_URL.md5" "$TMP_DIR/taxdump.tar.gz.md5"; then
        EXPECTED=$(awk '{print $1}' "$TMP_DIR/taxdump.tar.gz.md5")
        ACTUAL=$( (md5sum "$TMP_DIR/taxdump.tar.gz" 2> /dev/null || md5 -r "$TMP_DIR/taxdump.tar.gz") | awk '{print $1}')
        [ "$EXPECTED" = "$ACTUAL" ] || fail "taxdump.tar.gz checksum mismatch (expected $EXPECTED, got $ACTUAL); please re-run"
        mkdir -p data/taxdump
        tar -xzf "$TMP_DIR/taxdump.tar.gz" -C data/taxdump
        echo -e "  Taxonomy database downloaded and verified ${GREEN}OK${NC}"
    else
        echo -e "${YELLOW}WARNING: Failed to download the taxonomy database (needs curl or wget)${NC}"
        echo -e "${YELLOW}You can download it manually:${NC}"
        echo -e "  mkdir -p data/taxdump"
        echo -e "  curl -O ${TAXDUMP_URL} && tar -xzf taxdump.tar.gz -C data/taxdump"
    fi
else
    echo -e "  Taxonomy database already present ${GREEN}OK${NC}"
fi

# Installation complete
echo ""
echo -e "${BLUE}======================================${NC}"
echo -e "${GREEN}Installation Complete!${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""
echo -e "Next steps:"
echo -e "  1. Activate the environment:  ${BLUE}conda activate ${ENV_NAME}${NC}"
echo -e "  2. Check the installation:    ${BLUE}nextflow run main.nf -profile test${NC}"
echo -e "  3. Customize config/contaminants_example.txt for your study"
echo -e "  4. See all parameters:        ${BLUE}nextflow run main.nf --help${NC}"
echo ""
echo -e "For questions and troubleshooting, see ${BLUE}README.md${NC}"
echo ""
