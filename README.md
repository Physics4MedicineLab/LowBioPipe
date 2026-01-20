# LowBioPipe

<div align="center">

<img src="logo_lowbiopipe.png" alt="LowBioPipe" width="400"/>

**A Nextflow wrapper for statistical analysis of low biomass microbiome data**

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](VERSION)

</div>

## Overview

LowBioPipe is a Nextflow-based computational framework that serves as a wrapper for processing outputs from the [nf-core/taxprofiler](https://nf-co.re/taxprofiler) pipeline, specifically designed for downstream analysis of low-biomass microbiome samples. In addition to facilitating data processing, LowBioPipe implements comprehensive statistical analyses, data-mining, and data visualization techniques, including contaminant filtering, multi-rank taxonomic abundance profiling, and diversity metrics with graphical representations.

## Pipeline Summary

1. **Contaminant Filtering** - Remove reads assigned to specified taxa from [Kraken2](https://github.com/DerrickWood/kraken2) output
2. **[Recentrifuge](https://github.com/khyox/recentrifuge) Analysis** - Refined taxonomic profiling with negative control integration
3. **Abundance Tables** - Generate taxa × samples matrices at species, genus, and phylum levels
4. **Diversity Analysis** - Alpha/beta diversity, PCoA ordination, clustering heatmaps, PERMANOVA

## Requirements

- Nextflow ≥ 23.04.0
- Docker, Singularity, or Conda
- [NCBI taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)

### Input: nf-core/taxprofiler Output

LowBioPipe requires output from nf-core/taxprofiler with Kraken2 per-read classifications:

```bash
nextflow run nf-core/taxprofiler \
    -profile docker \
    --input samplesheet.csv \
    --databases databases.csv \
    --outdir taxprofiler_results \
    --run_kraken2 \
    --kraken2_save_readclassifications
```

The `--kraken2_save_readclassifications` flag is required.

## Installation

Clone the repository and run the installation script:

```bash
git clone https://github.com/Physics4MedicineLab/LowBioPipe.git
cd LowBioPipe
./install.sh
```

## Usage

```bash
nextflow run main.nf \
    --taxprofiler_results /path/to/taxprofiler/output \
    --kraken2_db_name <database_name> \
    --samples S1,S2,S3,S4 \
    --controls BLANK1,BLANK2 \
    -profile docker
```

For PERMANOVA statistical testing, provide a groups file:

```bash
nextflow run main.nf \
    --taxprofiler_results /path/to/taxprofiler/output \
    --kraken2_db_name <database_name> \
    --samples S1,S2,S3,S4 \
    --controls BLANK1,BLANK2 \
    --groups_file groups.tsv \
    -profile docker
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--taxprofiler_results` | Path to taxprofiler output directory | *required* |
| `--kraken2_db_name` | Kraken2 database name used in taxprofiler | *required* |
| `--samples` | Sample IDs (comma-separated) | `[]` |
| `--controls` | Negative control IDs (comma-separated) | `[]` |
| `--groups_file` | Sample-to-group mapping for PERMANOVA | `null` |
| `--taxdump` | NCBI taxdump directory | `data/taxdump` |
| `--exclude_taxa` | Contaminant TaxID file | `config/contaminants_example.txt` |
| `--abundance_ranks` | Taxonomic ranks | `['species','genus','phylum']` |
| `--abundance_aggregate` | Aggregate taxa to rank | `true` |
| `--abundance_min_samples` | Minimum samples threshold | `2` |
| `--outdir` | Output directory | `results` |

## Output

```
results/
├── filtered_reads/          # Filtered Kraken2 classifications
├── recentrifuge/            # Recentrifuge reports (.html, .xlsx)
├── abundance/               # Count matrices and metadata
│   ├── *_counts_*.tsv
│   ├── *_relative_*.tsv
│   └── *_taxa_metadata_*.tsv
└── diversity/
    └── {species,genus,phylum}/
        ├── alpha_metrics_*.tsv
        ├── beta_*.tsv
        ├── pcoa_*.png
        ├── heatmap_*.png
        └── permanova_*.txt
```

## Input File Formats

**Contaminant file** - One TaxID per line:
```
9606        # Homo sapiens
1751056     # Flavobacterium ammonificans
```

**Groups file** - Tab-separated, no header:
```
sample1    group_A
sample2    group_A
sample3    group_B
sample4    group_B
```

## Citation

If you use LowBioPipe, please cite:

```bibtex

```

Please also cite the underlying tools:

- **nf-core/taxprofiler**: Stamouli et al. (2023) [doi:10.1101/2023.10.20.563221](https://doi.org/10.1101/2023.10.20.563221)
- **Recentrifuge**: Martí (2019) [doi:10.1371/journal.pcbi.1006967](https://doi.org/10.1371/journal.pcbi.1006967)
- **Kraken2**: Wood et al. (2019) [doi:10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)

## License

MIT License. See [LICENSE](LICENSE) for details.
