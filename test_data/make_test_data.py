#!/usr/bin/env python3
"""
make_test_data.py - Regenerate the synthetic test dataset used by -profile test and CI

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Writes, next to this script:
    taxdump/                   mini NCBI taxdump (nodes.dmp, names.dmp) with real TaxIDs
    taxprofiler_output/        Kraken2 per-read files in nf-core/taxprofiler's layout
    contaminants.txt           exclusion list (with an indented comment and a blank line)
    groups.tsv                 sample-to-group mapping

Samples include IDs with underscores (patient_01, patient_02) and one negative
control (BLANK1), the cases that used to break sample matching.

Usage:
    python test_data/make_test_data.py
"""

import random
from pathlib import Path

OUT = Path(__file__).resolve().parent

# taxid, parent, rank, scientific name
TAXA = [
    (1, 1, "no rank", "root"),
    (131567, 1, "no rank", "cellular organisms"),
    (2, 131567, "superkingdom", "Bacteria"),
    (1224, 2, "phylum", "Pseudomonadota"),
    (1236, 1224, "class", "Gammaproteobacteria"),
    (72274, 1236, "order", "Pseudomonadales"),
    (135621, 72274, "family", "Pseudomonadaceae"),
    (286, 135621, "genus", "Pseudomonas"),
    (287, 286, "species", "Pseudomonas aeruginosa"),
    (316, 286, "species", "Pseudomonas stutzeri"),
    (91347, 1236, "order", "Enterobacterales"),
    (543, 91347, "family", "Enterobacteriaceae"),
    (561, 543, "genus", "Escherichia"),
    (562, 561, "species", "Escherichia coli"),
    (1239, 2, "phylum", "Bacillota"),
    (91061, 1239, "class", "Bacilli"),
    (1385, 91061, "order", "Bacillales"),
    (90964, 1385, "family", "Staphylococcaceae"),
    (1279, 90964, "genus", "Staphylococcus"),
    (1280, 1279, "species", "Staphylococcus aureus"),
    (46170, 1280, "subspecies", "Staphylococcus aureus subsp. aureus"),
    (1282, 1279, "species", "Staphylococcus epidermidis"),
    (201174, 2, "phylum", "Actinomycetota"),
    (1760, 201174, "class", "Actinomycetes"),
    (85007, 1760, "order", "Mycobacteriales"),
    (1653, 85007, "family", "Corynebacteriaceae"),
    (1716, 1653, "genus", "Corynebacterium"),
    (1719, 1716, "species", "Corynebacterium pseudotuberculosis"),
    (2759, 131567, "superkingdom", "Eukaryota"),
    (33208, 2759, "kingdom", "Metazoa"),
    (7711, 33208, "phylum", "Chordata"),
    (40674, 7711, "class", "Mammalia"),
    (9443, 40674, "order", "Primates"),
    (9604, 9443, "family", "Hominidae"),
    (9605, 9604, "genus", "Homo"),
    (9606, 9605, "species", "Homo sapiens"),
]

# sample -> {taxid: reads}; taxid 0 = unclassified
PROFILES = {
    "S1": {287: 400, 1280: 150, 46170: 50, 286: 30, 2: 20, 9606: 100, 316: 60, 0: 40},
    "S2": {562: 350, 1282: 120, 1719: 80, 9606: 80, 316: 40, 0: 30},
    "S3": {287: 200, 562: 200, 1282: 60, 1719: 20, 316: 30, 0: 10},
    "patient_01": {287: 500, 1280: 20, 316: 10, 0: 20},
    "patient_02": {562: 30, 1719: 600, 1282: 15, 316: 10, 0: 20},
    "BLANK1": {316: 300, 1282: 40, 9606: 50, 0: 60},
}

GROUPS = {"S1": "A", "S2": "B", "S3": "A", "patient_01": "B", "patient_02": "A"}

CONTAMINANTS = """\
# Test exclusion list
9606        # Homo sapiens
    # indented comment line

316         # Pseudomonas stutzeri
"""


def main() -> None:
    """Write the test dataset next to this script."""
    rng = random.Random(0)

    taxdump = OUT / "taxdump"
    taxdump.mkdir(exist_ok=True)
    with open(taxdump / "nodes.dmp", "w") as f:
        for taxid, parent, rank, _ in TAXA:
            f.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n")
    with open(taxdump / "names.dmp", "w") as f:
        for taxid, _, _, name in TAXA:
            f.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    reads_dir = OUT / "taxprofiler_output" / "kraken2" / "testdb"
    reads_dir.mkdir(parents=True, exist_ok=True)
    for sample, profile in PROFILES.items():
        lines = []
        for taxid, n_reads in profile.items():
            for i in range(n_reads):
                status = "U" if taxid == 0 else "C"
                kmers = rng.randint(5, 60)
                lines.append(
                    f"{status}\t{sample}_r{taxid}_{i}\t{taxid}\t150\t{taxid}:{kmers} 0:{100 - kmers}\n"
                )
        rng.shuffle(lines)
        (reads_dir / f"{sample}_testdb.kraken2.kraken2.classifiedreads.txt").write_text(
            "".join(lines)
        )

    (OUT / "contaminants.txt").write_text(CONTAMINANTS)
    (OUT / "groups.tsv").write_text(
        "".join(f"{sample}\t{group}\n" for sample, group in GROUPS.items())
    )
    print(f"Test dataset written to {OUT}")


if __name__ == "__main__":
    main()
