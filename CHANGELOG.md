# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-10-12

### Added

- Parameter validation and `--help`.
- `--abundance_source` to build the abundance tables from raw or control-filtered counts.
- PERMDISP next to PERMANOVA, and plots coloured by group.
- `-profile test` to check an installation on a small bundled dataset.

### Changed

- **Breaking:** Nextflow ≥ 25.10.0 is required.
- **Breaking:** `--samples`, `--controls` and `--abundance_ranks` take comma-separated strings, and each ID must match exactly one Kraken2 file.
- Descendants of contaminant taxa are removed only with `--filter_include_descendants` (Recentrifuge used to remove them always).
- Recentrifuge tables are written as TSV instead of Excel, so studies with more than ~500 samples work.

### Removed

- **Breaking:** `--max_cpus`, `--max_memory` and `--max_time`; set `process.resourceLimits` in a custom config instead.

### Fixed

- Abundance tables could merge samples (e.g. `patient_01` and `patient_02`), read the wrong Recentrifuge counts and include controls as samples. Results from 0.1.x should be regenerated.
- The documented `--samples S1,S2,...` usage failed, and `--groups_file` always crashed.
- PERMANOVA treated ungrouped samples as an extra group and was not reproducible with a fixed `--seed`.
- `-profile docker` and `-profile singularity` failed because of missing container images.

## [0.1.1] - 2026-03-09

### Fixed

- Input validation.

## [0.1.0] - 2026-01-21

### Added

- First release: contaminant filtering of Kraken2 per-read output from nf-core/taxprofiler, Recentrifuge analysis with negative controls, multi-rank abundance tables, and diversity analysis (alpha and beta diversity, PCoA, heatmaps, PERMANOVA).

[unreleased]: https://github.com/Physics4MedicineLab/LowBioPipe/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Physics4MedicineLab/LowBioPipe/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/Physics4MedicineLab/LowBioPipe/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Physics4MedicineLab/LowBioPipe/releases/tag/v0.1.0
