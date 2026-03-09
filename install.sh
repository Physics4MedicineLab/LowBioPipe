#!/bin/bash
# LowBioPipe Installation Script
# Installs dependencies and sets up the pipeline

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}======================================${NC}"
echo -e "${BLUE}  LowBioPipe Installation${NC}"
echo -e "${BLUE}======================================${NC}"
echo ""

# Check Python version
echo -e "${GREEN}[1/4]${NC} Checking Python version..."
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: Python 3 is not installed${NC}"
    exit 1
fi

PYTHON_VERSION=$(python3 --version | cut -d' ' -f2 | cut -d'.' -f1,2)
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d'.' -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d'.' -f2)

if [ "$PYTHON_MAJOR" -lt 3 ] || ([ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 8 ]); then
    echo -e "${RED}ERROR: Python >= 3.8 required, found $PYTHON_VERSION${NC}"
    exit 1
fi

echo -e "  Found Python $PYTHON_VERSION ${GREEN}OK${NC}"

# Install Python dependencies
echo ""
echo -e "${GREEN}[2/4]${NC} Installing Python dependencies..."
if pip3 install -r requirements.txt; then
    echo -e "  Python dependencies installed ${GREEN}OK${NC}"
else
    echo -e "${YELLOW}WARNING: Some dependencies failed to install${NC}"
    echo -e "${YELLOW}You may need to install them manually${NC}"
fi

# Make scripts executable
echo ""
echo -e "${GREEN}[3/4]${NC} Making scripts executable..."
chmod +x bin/*.py
echo -e "  Scripts are now executable ${GREEN}OK${NC}"

# Download NCBI taxonomy database
echo ""
echo -e "${GREEN}[4/4]${NC} Checking NCBI taxonomy database..."
if [ ! -f "data/taxdump/nodes.dmp" ]; then
    echo -e "${YELLOW}  Taxonomy database not found. Downloading...${NC}"
    mkdir -p data/taxdump
    cd data/taxdump

    echo -e "  Downloading taxdump.tar.gz (~50 MB)..."
    if wget -q --show-progress https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz; then
        echo -e "  Extracting..."
        tar -xzf taxdump.tar.gz
        rm taxdump.tar.gz
        cd ../..
        echo -e "  Taxonomy database downloaded ${GREEN}OK${NC}"
    else
        cd ../..
        echo -e "${YELLOW}WARNING: Failed to download taxonomy database${NC}"
        echo -e "${YELLOW}You can download it manually:${NC}"
        echo -e "  mkdir -p data/taxdump"
        echo -e "  cd data/taxdump"
        echo -e "  wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        echo -e "  tar -xzf taxdump.tar.gz"
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
echo -e "  1. Review the README.md for usage instructions"
echo -e "  2. Customize config/contaminants_example.txt for your study"
echo -e "  3. Run the pipeline on your data:"
echo -e ""
echo -e "     ${BLUE}bin/filter_taxa.py --help${NC}"
echo -e "     ${BLUE}bin/rcf_to_abundance.py --help${NC}"
echo -e ""
echo -e "For questions and troubleshooting, see ${BLUE}README.md${NC}"
echo ""
