#!/usr/bin/env python3
"""
rcf_to_abundance.py - Convert Recentrifuge detailed output to abundance tables

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Reads unique assignments (UNA counts) from Recentrifuge's detailed output, either
TSV (rcf -e TSV: <base>.data.tsv + <base>.stat.tsv) or Excel (FULL sheet), and
creates taxa x samples abundance matrices for specified taxonomic ranks. Names,
ranks and lineages come from the same NCBI taxdump given to Recentrifuge.

Usage:
    python rcf_to_abundance.py recentrifuge_output.rcf.data.tsv --taxdump taxdump --min-samples 2 --rank species

    # Aggregate species to phylum level:
    python rcf_to_abundance.py recentrifuge_output.rcf.data.tsv --taxdump taxdump --rank phylum --aggregate

    # Use Recentrifuge's control-filtered counts (requires a run with negative controls):
    python rcf_to_abundance.py recentrifuge_output.rcf.data.tsv --taxdump taxdump --source ctrl --controls BLANK1
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="[%(levelname)s] %(message)s", stream=sys.stderr
)
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed argument namespace
    """
    ap = argparse.ArgumentParser(
        description="Build abundance table (taxa x samples) from Recentrifuge detailed output (UNA counts)."
    )
    ap.add_argument(
        "rcf_file",
        help="Recentrifuge .rcf.data.tsv (with .rcf.stat.tsv next to it) or .rcf.xlsx file",
    )
    ap.add_argument(
        "--taxdump",
        type=Path,
        required=True,
        help="NCBI taxdump directory (nodes.dmp, names.dmp) used for the Recentrifuge run",
    )
    ap.add_argument(
        "--source",
        choices=["raw", "ctrl"],
        default="raw",
        help="Counts to use: 'raw' sample counts, or Recentrifuge's control-filtered "
        "'<sample>_CTRL_SUMMARY' counts (controls are then left out). Default: raw",
    )
    ap.add_argument(
        "--controls",
        nargs="*",
        help="Negative-control IDs; with --source ctrl, any other sample without "
        "control-filtered counts is reported as dropped",
    )
    ap.add_argument(
        "--rank",
        default="species",
        help="Rank to filter/aggregate to (e.g., species, genus, family, phylum, any). Default: species",
    )
    ap.add_argument(
        "--aggregate",
        action="store_true",
        help="Aggregate all taxa UP to specified rank (e.g., sum all species into their phyla). "
        "Without this flag, only taxa directly assigned at --rank are included.",
    )
    ap.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Drop taxa with total count < threshold (after pivot). Default 1",
    )
    ap.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Keep only taxa present in >= this number of samples (default 1)",
    )
    ap.add_argument(
        "--out-prefix", default="abundance", help="Output prefix (default: abundance)"
    )
    ap.add_argument(
        "--relative",
        action="store_true",
        help="Also write relative-abundance table (per-sample proportions)",
    )
    ap.add_argument(
        "--keep-unranked",
        action="store_true",
        help="Keep taxa without resolvable rank/TaxID when --rank != any (default: drop)",
    )
    return ap.parse_args()


def load_taxonomy(
    taxdump: Path,
) -> tuple[dict[int, int], dict[int, str], dict[int, str]]:
    """
    Load parent, rank and scientific-name maps from an NCBI taxdump directory.

    Args:
        taxdump: Directory containing nodes.dmp and names.dmp

    Returns:
        Tuple of (taxid -> parent taxid, taxid -> rank, taxid -> scientific name)
    """
    parent: dict[int, int] = {}
    rank: dict[int, str] = {}
    name: dict[int, str] = {}
    for dmp in ("nodes.dmp", "names.dmp"):
        if not (taxdump / dmp).exists():
            logger.error("%s not found in taxdump directory %s", dmp, taxdump)
            sys.exit(1)

    with open(taxdump / "nodes.dmp") as f:
        for line in f:
            fields = line.split("\t|\t", 3)
            taxid = int(fields[0])
            parent[taxid] = int(fields[1])
            rank[taxid] = fields[2]

    with open(taxdump / "names.dmp") as f:
        for line in f:
            fields = line.rstrip("\t|\n").split("\t|\t")
            if fields[3] == "scientific name":
                name[int(fields[0])] = fields[1]

    logger.info("Loaded %d taxa from %s", len(parent), taxdump)
    return parent, rank, name


def ancestor_at_rank(
    taxid: int, target_rank: str, parent: dict[int, int], rank: dict[int, str]
) -> int | None:
    """
    Walk up the taxonomy from taxid (inclusive) to the first node at target_rank.

    Args:
        taxid: Starting taxid
        target_rank: Rank to find (e.g., 'phylum', 'genus')
        parent: taxid -> parent taxid map
        rank: taxid -> rank map

    Returns:
        Taxid of the ancestor at target_rank, or None if the lineage has none
    """
    seen = set()
    while taxid in parent and taxid not in seen:
        if rank.get(taxid) == target_rank:
            return taxid
        seen.add(taxid)
        taxid = parent[taxid]
    return None


def read_recentrifuge(
    rcf_path: str, source: str = "raw", controls: list[str] | None = None
) -> pd.DataFrame:
    """
    Read UNA counts per sample from Recentrifuge's detailed output.

    Accepts the TSV output (<base>.data.tsv, with <base>.stat.tsv next to it) or
    the FULL sheet of the Excel output. Both have one block of columns (cnt, una,
    sco) per sample under a sample-name header, and list Recentrifuge's derived
    pseudo-samples (<sample>_CTRL_<rank>, <sample>_EXCLUSIVE_*, SHARED_*,
    CONTROL_SHARED_*, ...) next to the raw ones. Raw sample names are taken from
    the sample statistics (.stat.tsv or '_sample_stats' sheet) and columns are
    selected by exact (sample, 'una') header match.

    Args:
        rcf_path: Path to Recentrifuge .data.tsv or .xlsx file
        source: 'raw' for raw sample counts, 'ctrl' for each sample's
            control-filtered '<sample>_CTRL_SUMMARY' counts
        controls: Negative-control IDs, to tell controls apart from samples that
            lost all their reads to control filtering

    Returns:
        DataFrame with 'Taxon' column and one UNA-count column per sample
    """
    try:
        if rcf_path.endswith(".data.tsv"):
            stats_path = rcf_path[: -len(".data.tsv")] + ".stat.tsv"
            raw_samples = [
                str(s)
                for s in pd.read_csv(stats_path, sep="\t", index_col=0, nrows=0).columns
            ]
            sheet = pd.read_csv(rcf_path, sep="\t", header=None, dtype=str)
        else:
            raw_samples = [
                str(s)
                for s in pd.read_excel(
                    rcf_path, sheet_name="_sample_stats", index_col=0, engine="openpyxl"
                ).columns
            ]
            sheet = pd.read_excel(
                rcf_path, sheet_name="FULL", header=None, engine="openpyxl"
            )
    except FileNotFoundError as e:
        logger.error("Recentrifuge output not found: %s", e.filename or rcf_path)
        sys.exit(1)
    except ValueError as e:  # includes a missing Excel sheet
        logger.error("Failed to read Recentrifuge output %s: %s", rcf_path, e)
        sys.exit(1)

    # Row 0: sample names (repeated, or merged in Excel, across their stat columns),
    # row 1: stat names, row 2: index name ('Id'), rows 3+: one taxon per row
    header_samples = sheet.iloc[0].ffill().astype(str)
    header_stats = sheet.iloc[1].astype(str)
    una_column = {
        sample: col
        for col, (sample, stat) in enumerate(zip(header_samples, header_stats))
        if stat == "una"
    }
    data = sheet.iloc[3:]

    result_df = pd.DataFrame({"Taxon": data[0].values})
    skipped = []
    for sample in raw_samples:
        column_name = sample if source == "raw" else f"{sample}_CTRL_SUMMARY"
        if column_name not in una_column:
            skipped.append(sample)
            continue
        result_df[sample] = (
            pd.to_numeric(data[una_column[column_name]], errors="coerce")
            .fillna(0)
            .astype(int)
            .values
        )

    if skipped:
        if source == "raw":
            logger.error("Raw samples missing from %s: %s", rcf_path, skipped)
            sys.exit(1)
        # Controls never get a CTRL summary; a sample loses it only if nothing survived filtering
        if controls is None:
            logger.info("No control-filtered counts (controls or empty samples): %s", skipped)
        else:
            left_out = [s for s in skipped if s in controls]
            dropped = [s for s in skipped if s not in controls]
            if left_out:
                logger.info("Controls left out of the tables: %s", left_out)
            if dropped:
                logger.warning(
                    "Samples dropped: all their reads were removed by control filtering: %s",
                    dropped,
                )

    if result_df.shape[1] == 1:
        logger.error("No '%s' sample columns found in %s", source, rcf_path)
        sys.exit(1)

    logger.info(
        "Read Recentrifuge output (%s counts): %d taxa x %d samples",
        source,
        len(result_df),
        result_df.shape[1] - 1,
    )
    return result_df


def main() -> None:
    """Main entry point for abundance table generation."""
    args = parse_args()

    # Load the same taxonomy Recentrifuge used
    parent, taxid2rank, taxid2name = load_taxonomy(args.taxdump)

    # Read the Recentrifuge output
    file_path = Path(args.rcf_file)
    if not file_path.exists():
        logger.error("File not found: %s", args.rcf_file)
        sys.exit(1)

    # Read clean dataframe: 'Taxon' column + sample columns with UNA counts
    df_clean = read_recentrifuge(args.rcf_file, args.source, args.controls)

    # Get sample names (all columns except 'Taxon')
    sample_names = [col for col in df_clean.columns if col != "Taxon"]

    target_rank = args.rank.lower()
    mode = "AGGREGATE" if args.aggregate else "FILTER"
    logger.info(
        "Processing %d taxa across %d samples...", len(df_clean), len(sample_names)
    )
    logger.info("Mode: %s at rank '%s'", mode, target_rank)

    records: list[dict[str, Any]] = []

    # Process each row (each taxon)
    for idx, row in df_clean.iterrows():
        # Get taxon ID (should be numeric)
        taxid_str = str(row["Taxon"]).strip()

        # Skip if not a valid taxon ID
        if not taxid_str or taxid_str == "nan" or taxid_str == "Id":
            continue

        try:
            taxid = int(taxid_str)
        except (ValueError, TypeError):
            continue

        # Get taxon name and rank from the taxdump
        taxon_name = taxid2name.get(taxid, f"TaxID_{taxid}")
        resolved_rank = taxid2rank.get(taxid)

        # Determine which taxon to use (original or ancestor)
        final_taxid = taxid
        final_name = taxon_name
        final_rank = resolved_rank

        # Handle aggregation mode
        if args.aggregate and target_rank != "any":
            # Find ancestor at target rank
            ancestor_taxid = ancestor_at_rank(taxid, target_rank, parent, taxid2rank)

            if ancestor_taxid is not None:
                # Use ancestor for aggregation
                final_taxid = ancestor_taxid
                final_name = taxid2name.get(ancestor_taxid, f"TaxID_{ancestor_taxid}")
                final_rank = target_rank
            else:
                # No ancestor at target rank found, skip this taxon
                continue
        else:
            # Filter mode: only keep taxa at exact rank
            keep_node = False
            if target_rank == "any":
                keep_node = True
            else:
                if resolved_rank == target_rank or args.keep_unranked and resolved_rank is None:
                    keep_node = True

            if not keep_node:
                continue

        # Extract counts for each sample
        for sample_name in sample_names:
            count = row.get(sample_name, 0)
            if pd.isna(count):
                count = 0
            else:
                try:
                    count = int(float(count))
                except (ValueError, TypeError):
                    count = 0

            if count > 0:
                records.append(
                    {
                        "TaxonName": final_name,
                        "TaxID": final_taxid,
                        "Rank": final_rank if final_rank else "",
                        "Sample": sample_name,
                        "Count": count,
                    }
                )

    if not records:
        logger.error("No taxa passed the filter for rank '%s'", target_rank)
        sys.exit(1)

    long_df = pd.DataFrame(records)

    def tax_label(row: pd.Series) -> str:
        if row["TaxID"] and row["TaxID"] != -1:
            return f"{row['TaxonName']} [{row['TaxID']}]"
        return row["TaxonName"]

    long_df["Taxon"] = long_df.apply(tax_label, axis=1)

    abundance = long_df.pivot_table(
        index="Taxon", columns="Sample", values="Count", aggfunc="sum", fill_value=0
    )

    if args.min_count > 1:
        abundance = abundance.loc[abundance.sum(axis=1) >= args.min_count]

    if args.min_samples > 1:
        presence = (abundance > 0).sum(axis=1)
        abundance = abundance.loc[presence >= args.min_samples]

    if abundance.empty:
        logger.error(
            "No taxa at rank '%s' pass --min-count %d and --min-samples %d",
            target_rank,
            args.min_count,
            args.min_samples,
        )
        sys.exit(1)

    out_counts = f"{args.out_prefix}_counts_{target_rank}.tsv"
    abundance.to_csv(out_counts, sep="\t")

    if args.relative:
        col_sums = abundance.sum(axis=0)
        zero_cols = col_sums[col_sums == 0].index.tolist()
        if zero_cols:
            logger.warning(
                "Dropping %d samples with zero total counts before relative "
                "abundance normalization: %s",
                len(zero_cols),
                zero_cols,
            )
            abundance_nz = abundance.drop(columns=zero_cols)
        else:
            abundance_nz = abundance
        rel = abundance_nz.div(abundance_nz.sum(axis=0), axis=1)
        out_rel = f"{args.out_prefix}_relative_{target_rank}.tsv"
        rel.to_csv(out_rel, sep="\t", float_format="%.6f")

    meta_rows: list[dict[str, Any]] = []
    for taxon_label, sub in long_df.groupby("Taxon"):
        if taxon_label not in abundance.index:
            continue
        r = sub.iloc[0]
        meta_rows.append(
            {
                "Taxon": taxon_label,
                "Name": r["TaxonName"],
                "TaxID": r["TaxID"] if r["TaxID"] != -1 else "",
                "Rank": r["Rank"],
            }
        )

    meta_df = pd.DataFrame(meta_rows).drop_duplicates().sort_values("Taxon")
    out_meta = f"{args.out_prefix}_taxa_metadata_{target_rank}.tsv"
    meta_df.to_csv(out_meta, sep="\t", index=False)

    logger.info("Successfully created:")
    logger.info("  - %s", out_counts)
    if args.relative:
        logger.info("  - %s", out_rel)
    logger.info("  - %s", out_meta)
    logger.info(
        "Filtered to %d taxa across %d samples", len(abundance), len(sample_names)
    )


if __name__ == "__main__":
    main()
