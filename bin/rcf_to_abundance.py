#!/usr/bin/env python3
"""
rcf_to_abundance.py - Convert Recentrifuge Excel output to abundance tables

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Reads the 'FULL' sheet with unique assignments (UNA counts) and creates
taxa x samples abundance matrices for specified taxonomic ranks.

Usage:
    python rcf_to_abundance.py recentrifuge_output.rcf.xlsx --min-samples 2 --rank species

    # Aggregate species to phylum level:
    python rcf_to_abundance.py recentrifuge_output.rcf.xlsx --rank phylum --aggregate
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Any

import pandas as pd

from ete3.ncbi_taxonomy import NCBITaxa

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
        description="Build abundance table (taxa x samples) from Recentrifuge Excel output (FULL sheet, UNA counts)."
    )
    ap.add_argument("rcf_file", help="Path to Recentrifuge .rcf.xlsx file")
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
        "--ncbi-update",
        action="store_true",
        help="Run NCBITaxa.update_taxonomy_database() before mapping",
    )
    ap.add_argument(
        "--keep-unranked",
        action="store_true",
        help="Keep taxa without resolvable rank/TaxID when --rank != any (default: drop)",
    )
    return ap.parse_args()


def get_sample_name_from_path(path_text: str) -> str:
    """
    Extract sample name from file path.

    Args:
        path_text: File path string

    Returns:
        Extracted sample name
    """
    tail = str(path_text).strip().split("/")[-1]
    return tail.split("_")[0] if "_" in tail else tail


def get_ancestor_at_rank(
    ncbi: NCBITaxa, taxid: int, target_rank: str, max_depth: int = 50
) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """
    Walk up the taxonomy tree to find an ancestor at the specified rank.

    Args:
        ncbi: NCBITaxa instance
        taxid: Starting taxid
        target_rank: Rank to find (e.g., 'phylum', 'genus')
        max_depth: Maximum steps to prevent infinite loops

    Returns:
        Tuple of (ancestor_taxid, ancestor_name, ancestor_rank) or (None, None, None) if not found
    """
    current = taxid
    for _ in range(max_depth):
        try:
            # Get rank of current taxid
            rank_dict = ncbi.get_rank([current])
            current_rank = rank_dict.get(current)

            # If we found the target rank, return this taxid
            if current_rank == target_rank:
                name_dict = ncbi.get_taxid_translator([current])
                name = name_dict.get(current, f"TaxID_{current}")
                return (current, name, current_rank)

            # Get parent
            lineage = ncbi.get_lineage(current)
            if not lineage or len(lineage) <= 1:
                break  # Reached root

            # Check if any ancestor in lineage has the target rank
            for ancestor_id in reversed(
                lineage[:-1]
            ):  # Exclude current, walk backwards
                rank_dict = ncbi.get_rank([ancestor_id])
                ancestor_rank = rank_dict.get(ancestor_id)
                if ancestor_rank == target_rank:
                    name_dict = ncbi.get_taxid_translator([ancestor_id])
                    name = name_dict.get(ancestor_id, f"TaxID_{ancestor_id}")
                    return (ancestor_id, name, ancestor_rank)

            break  # No ancestor with target rank found in lineage

        except (KeyError, ValueError, TypeError) as e:
            logger.debug("Error looking up ancestor for taxid %d: %s", taxid, e)
            break

    return (None, None, None)


def read_recentrifuge_excel(excel_path: str, sheet_name: str = "FULL") -> pd.DataFrame:
    """
    Read Recentrifuge Excel file from FULL sheet with unique assignments.

    The FULL sheet has structure:
    - Column 0: "Samples" (taxon IDs)
    - For each sample: 3 columns (sample_name, UNA_count, other)
    - Row 0: metadata ("Stats", "cnt", etc.) - skip this
    - Row 1+: actual data

    Args:
        excel_path: Path to Recentrifuge Excel file
        sheet_name: Sheet name to read (default: 'FULL')

    Returns:
        DataFrame with 'Taxon' column and sample columns with UNA counts
    """
    try:
        # Read the sheet
        df_raw = pd.read_excel(excel_path, sheet_name=sheet_name, engine="openpyxl")

        # Skip first row if it contains "Stats" metadata
        if len(df_raw) > 0 and df_raw.iloc[0, 0] == "Stats":
            df_raw = df_raw.iloc[1:].reset_index(drop=True)

        # Identify sample columns and their UNA count columns
        sample_data: Dict[str, Any] = {}
        taxon_col = df_raw.columns[0]  # Should be "Samples"

        i = 1  # Start from column 1 (after "Samples")
        while i < len(df_raw.columns):
            col = df_raw.columns[i]

            # Check if this is a sample column (not "Unnamed")
            if not str(col).startswith("Unnamed"):
                sample_name = get_sample_name_from_path(col)

                # The UNA count is in the next column (first Unnamed after sample name)
                if i + 1 < len(df_raw.columns):
                    una_col = df_raw.columns[i + 1]
                    sample_data[sample_name] = una_col
                    i += 3  # Skip to next sample (sample + 2 unnamed columns)
                else:
                    i += 1
            else:
                i += 1

        # Validate detected structure
        if not sample_data:
            logger.error(
                "No samples detected in sheet '%s'. "
                "Expected columns with sample file paths followed by 'Unnamed' columns.",
                sheet_name,
            )
            sys.exit(1)

        total_data_cols = len(df_raw.columns) - 1  # Subtract taxon column
        expected_cols = len(sample_data) * 3
        if total_data_cols != expected_cols:
            logger.warning(
                "Column count mismatch: %d data columns for %d samples "
                "(expected %d = %d samples x 3). "
                "Recentrifuge Excel format may have changed.",
                total_data_cols,
                len(sample_data),
                expected_cols,
                len(sample_data),
            )

        # Create clean dataframe with taxon IDs and UNA counts
        result_df = pd.DataFrame()
        result_df["Taxon"] = df_raw[taxon_col]

        for sample_name, una_col in sample_data.items():
            result_df[sample_name] = (
                pd.to_numeric(df_raw[una_col], errors="coerce").fillna(0).astype(int)
            )

        logger.info(
            "Successfully read sheet '%s': %d taxa x %d samples",
            sheet_name,
            len(result_df),
            len(sample_data),
        )
        return result_df

    except FileNotFoundError:
        logger.error("Excel file not found: %s", excel_path)
        sys.exit(1)
    except KeyError as e:
        # Try to list available sheets
        try:
            import openpyxl

            wb = openpyxl.load_workbook(excel_path, read_only=True)
            available_sheets = wb.sheetnames
            logger.error(
                "Sheet '%s' not found. Available sheets: %s",
                sheet_name,
                available_sheets,
            )
        except ImportError:
            logger.error("Sheet '%s' not found: %s", sheet_name, e)
        sys.exit(1)
    except ValueError as e:
        logger.error("Failed to parse Excel file: %s", e)
        sys.exit(1)


def get_taxon_info(
    ncbi: NCBITaxa,
    taxid: int,
    name_cache: Dict[int, str],
    rank_cache: Dict[int, Optional[str]],
) -> Tuple[str, Optional[str]]:
    """
    Get taxon name and rank from NCBI, using caches for efficiency.

    Args:
        ncbi: NCBITaxa instance
        taxid: Taxonomy ID
        name_cache: Cache for taxon names
        rank_cache: Cache for taxon ranks

    Returns:
        Tuple of (taxon_name, taxon_rank)
    """
    # Get name
    if taxid not in name_cache:
        try:
            name_dict = ncbi.get_taxid_translator([taxid])
            name_cache[taxid] = name_dict.get(taxid, f"TaxID_{taxid}")
        except (KeyError, ValueError) as e:
            logger.debug("Failed to get name for taxid %d: %s", taxid, e)
            name_cache[taxid] = f"TaxID_{taxid}"

    # Get rank
    if taxid not in rank_cache:
        try:
            rank_dict = ncbi.get_rank([taxid])
            rank_cache[taxid] = rank_dict.get(taxid, None)
        except (KeyError, ValueError) as e:
            logger.debug("Failed to get rank for taxid %d: %s", taxid, e)
            rank_cache[taxid] = None

    return name_cache[taxid], rank_cache[taxid]


def main() -> None:
    """Main entry point for abundance table generation."""
    args = parse_args()

    # Initialize NCBI taxonomy database
    ncbi = NCBITaxa()
    if args.ncbi_update:
        try:
            ncbi.update_taxonomy_database()
        except Exception as e:
            logger.warning("Warning updating NCBITaxa DB: %s", e)

    # Read the Excel file
    file_path = Path(args.rcf_file)
    if not file_path.exists():
        logger.error("File not found: %s", args.rcf_file)
        sys.exit(1)

    # Read clean dataframe: 'Taxon' column + sample columns with UNA counts
    df_clean = read_recentrifuge_excel(args.rcf_file)

    # Get sample names (all columns except 'Taxon')
    sample_names = [col for col in df_clean.columns if col != "Taxon"]

    target_rank = args.rank.lower()
    mode = "AGGREGATE" if args.aggregate else "FILTER"
    logger.info(
        "Processing %d taxa across %d samples...", len(df_clean), len(sample_names)
    )
    logger.info("Mode: %s at rank '%s'", mode, target_rank)

    taxid2rank_cache: Dict[int, Optional[str]] = {}
    taxid2name_cache: Dict[int, str] = {}

    records: List[Dict[str, Any]] = []

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

        # Get taxon name and rank from NCBI
        taxon_name, resolved_rank = get_taxon_info(
            ncbi, taxid, taxid2name_cache, taxid2rank_cache
        )

        # Determine which taxon to use (original or ancestor)
        final_taxid = taxid
        final_name = taxon_name
        final_rank = resolved_rank

        # Handle aggregation mode
        if args.aggregate and target_rank != "any":
            # Find ancestor at target rank
            ancestor_taxid, ancestor_name, ancestor_rank = get_ancestor_at_rank(
                ncbi, taxid, target_rank
            )

            if ancestor_taxid is not None:
                # Use ancestor for aggregation
                final_taxid = ancestor_taxid
                final_name = ancestor_name
                final_rank = ancestor_rank
            else:
                # No ancestor at target rank found, skip this taxon
                continue
        else:
            # Filter mode: only keep taxa at exact rank
            keep_node = False
            if target_rank == "any":
                keep_node = True
            else:
                if resolved_rank == target_rank:
                    keep_node = True
                elif args.keep_unranked and resolved_rank is None:
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

    meta_rows: List[Dict[str, Any]] = []
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
