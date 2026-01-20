#!/usr/bin/env python3
"""
rcf_to_abundance.py - Convert Recentrifuge Excel output to abundance tables

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Reads the 'FULL' sheet with unique assignments (UNA counts) and creates
taxa × samples abundance matrices for specified taxonomic ranks.

Usage:
    python rcf_to_abundance.py recentrifuge_output.rcf.xlsx --min-samples 2 --rank species

    # Aggregate species to phylum level:
    python rcf_to_abundance.py recentrifuge_output.rcf.xlsx --rank phylum --aggregate
"""
import sys
import argparse
from pathlib import Path

import pandas as pd

from ete3.ncbi_taxonomy import NCBITaxa

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build abundance table (taxa x samples) from Recentrifuge Excel output (FULL sheet, UNA counts)."
    )
    ap.add_argument("rcf_file", help="Path to Recentrifuge .rcf.xlsx file")
    ap.add_argument("--rank", default="species",
                    help="Rank to filter/aggregate to (e.g., species, genus, family, phylum, any). Default: species")
    ap.add_argument("--aggregate", action="store_true",
                    help="Aggregate all taxa UP to specified rank (e.g., sum all species into their phyla). "
                         "Without this flag, only taxa directly assigned at --rank are included.")
    ap.add_argument("--min-count", type=int, default=1,
                    help="Drop taxa with total count < threshold (after pivot). Default 1")
    ap.add_argument("--min-samples", type=int, default=1,
                    help="Keep only taxa present in >= this number of samples (default 1)")
    ap.add_argument("--out-prefix", default="abundance",
                    help="Output prefix (default: abundance)")
    ap.add_argument("--relative", action="store_true",
                    help="Also write relative-abundance table (per-sample proportions)")
    ap.add_argument("--ncbi-update", action="store_true",
                    help="Run NCBITaxa.update_taxonomy_database() before mapping")
    ap.add_argument("--keep-unranked", action="store_true",
                    help="Keep taxa without resolvable rank/TaxID when --rank != any (default: drop)")
    return ap.parse_args()

def get_sample_name_from_path(path_text: str) -> str:
    """Extract sample name from file path."""
    tail = str(path_text).strip().split("/")[-1]
    return tail.split("_")[0] if "_" in tail else tail

def get_ancestor_at_rank(ncbi, taxid: int, target_rank: str, max_depth: int = 50):
    """
    Walk up the taxonomy tree to find an ancestor at the specified rank.

    Args:
        ncbi: NCBITaxa instance
        taxid: Starting taxid
        target_rank: Rank to find (e.g., 'phylum', 'genus')
        max_depth: Maximum steps to prevent infinite loops

    Returns:
        tuple: (ancestor_taxid, ancestor_name, ancestor_rank) or (None, None, None) if not found
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
            for ancestor_id in reversed(lineage[:-1]):  # Exclude current, walk backwards
                rank_dict = ncbi.get_rank([ancestor_id])
                ancestor_rank = rank_dict.get(ancestor_id)
                if ancestor_rank == target_rank:
                    name_dict = ncbi.get_taxid_translator([ancestor_id])
                    name = name_dict.get(ancestor_id, f"TaxID_{ancestor_id}")
                    return (ancestor_id, name, ancestor_rank)

            break  # No ancestor with target rank found in lineage

        except Exception:
            break

    return (None, None, None)

def read_recentrifuge_excel(excel_path, sheet_name='FULL'):
    """
    Read Recentrifuge Excel file from FULL sheet with unique assignments.

    The FULL sheet has structure:
    - Column 0: "Samples" (taxon IDs)
    - For each sample: 3 columns (sample_name, UNA_count, other)
    - Row 0: metadata ("Stats", "cnt", etc.) - skip this
    - Row 1+: actual data
    """
    try:
        # Read the sheet
        df_raw = pd.read_excel(excel_path, sheet_name=sheet_name, engine='openpyxl')

        # Skip first row if it contains "Stats" metadata
        if df_raw.iloc[0, 0] == "Stats":
            df_raw = df_raw.iloc[1:].reset_index(drop=True)

        # Identify sample columns and their UNA count columns
        sample_data = {}
        taxon_col = df_raw.columns[0]  # Should be "Samples"

        i = 1  # Start from column 1 (after "Samples")
        while i < len(df_raw.columns):
            col = df_raw.columns[i]

            # Check if this is a sample column (not "Unnamed")
            if not str(col).startswith('Unnamed'):
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

        # Create clean dataframe with taxon IDs and UNA counts
        result_df = pd.DataFrame()
        result_df['Taxon'] = df_raw[taxon_col]

        for sample_name, una_col in sample_data.items():
            result_df[sample_name] = pd.to_numeric(df_raw[una_col], errors='coerce').fillna(0).astype(int)

        print(f"Successfully read sheet '{sheet_name}': {len(result_df)} taxa × {len(sample_data)} samples", file=sys.stderr)
        return result_df

    except Exception as e:
        # Try to list available sheets
        try:
            import openpyxl
            wb = openpyxl.load_workbook(excel_path, read_only=True)
            available_sheets = wb.sheetnames
            sys.exit(f"ERROR: Failed to read sheet '{sheet_name}' from Excel file.\nAvailable sheets: {available_sheets}\nError: {e}")
        except:
            sys.exit(f"ERROR: Failed to read Excel file: {e}\nMake sure openpyxl is installed: pip install openpyxl")

def main():
    args = parse_args()

    # Initialize NCBI taxonomy database
    ncbi = NCBITaxa()
    if args.ncbi_update:
        try:
            ncbi.update_taxonomy_database()
        except Exception as e:
            print(f"[!] Warning updating NCBITaxa DB: {e}", file=sys.stderr)

    # Read the Excel file
    file_path = Path(args.rcf_file)
    if not file_path.exists():
        sys.exit(f"ERROR: File not found: {args.rcf_file}")

    # Read clean dataframe: 'Taxon' column + sample columns with UNA counts
    df_clean = read_recentrifuge_excel(args.rcf_file)

    # Get sample names (all columns except 'Taxon')
    sample_names = [col for col in df_clean.columns if col != 'Taxon']

    target_rank = args.rank.lower()
    mode = "AGGREGATE" if args.aggregate else "FILTER"
    print(f"Processing {len(df_clean)} taxa across {len(sample_names)} samples...", file=sys.stderr)
    print(f"Mode: {mode} at rank '{target_rank}'", file=sys.stderr)
    taxid2rank_cache = {}
    taxid2name_cache = {}

    records = []

    # Process each row (each taxon)
    for idx, row in df_clean.iterrows():
        # Get taxon ID (should be numeric)
        taxid_str = str(row['Taxon']).strip()

        # Skip if not a valid taxon ID
        if not taxid_str or taxid_str == 'nan' or taxid_str == 'Id':
            continue

        try:
            taxid = int(taxid_str)
        except (ValueError, TypeError):
            continue

        # Get taxon name and rank from NCBI
        if taxid not in taxid2name_cache:
            try:
                name_dict = ncbi.get_taxid_translator([taxid])
                taxid2name_cache[taxid] = name_dict.get(taxid, f"TaxID_{taxid}")
            except Exception:
                taxid2name_cache[taxid] = f"TaxID_{taxid}"

        if taxid not in taxid2rank_cache:
            try:
                rank_dict = ncbi.get_rank([taxid])
                taxid2rank_cache[taxid] = rank_dict.get(taxid, None)
            except Exception:
                taxid2rank_cache[taxid] = None

        taxon_name = taxid2name_cache[taxid]
        resolved_rank = taxid2rank_cache[taxid]

        # Determine which taxon to use (original or ancestor)
        final_taxid = taxid
        final_name = taxon_name
        final_rank = resolved_rank

        # Handle aggregation mode
        if args.aggregate and target_rank != "any":
            # Find ancestor at target rank
            ancestor_taxid, ancestor_name, ancestor_rank = get_ancestor_at_rank(ncbi, taxid, target_rank)

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
                records.append({
                    "TaxonName": final_name,
                    "TaxID": final_taxid,
                    "Rank": final_rank if final_rank else "",
                    "Sample": sample_name,
                    "Count": count
                })

    if not records:
        sys.exit(f"ERROR: No taxa passed the filter for rank '{target_rank}'.")

    long_df = pd.DataFrame(records)

    def tax_label(row):
        if row["TaxID"] and row["TaxID"] != -1:
            return f"{row['TaxonName']} [{row['TaxID']}]"
        return row["TaxonName"]

    long_df["Taxon"] = long_df.apply(tax_label, axis=1)

    abundance = long_df.pivot_table(index="Taxon",
                                    columns="Sample",
                                    values="Count",
                                    aggfunc="sum",
                                    fill_value=0)

    if args.min_count > 1:
        abundance = abundance.loc[abundance.sum(axis=1) >= args.min_count]

    if args.min_samples > 1:
        presence = (abundance > 0).sum(axis=1)
        abundance = abundance.loc[presence >= args.min_samples]

    out_counts = f"{args.out_prefix}_counts_{target_rank}.tsv"
    abundance.to_csv(out_counts, sep="\t")

    if args.relative:
        rel = abundance.div(abundance.sum(axis=0), axis=1)
        out_rel = f"{args.out_prefix}_relative_{target_rank}.tsv"
        rel.to_csv(out_rel, sep="\t", float_format="%.6f")

    meta_rows = []
    for taxon_label, sub in long_df.groupby("Taxon"):
        if taxon_label not in abundance.index:
            continue
        r = sub.iloc[0]
        meta_rows.append({
            "Taxon": taxon_label,
            "Name": r["TaxonName"],
            "TaxID": r["TaxID"] if r["TaxID"] != -1 else "",
            "Rank": r["Rank"]
        })

    meta_df = pd.DataFrame(meta_rows).drop_duplicates().sort_values("Taxon")
    out_meta = f"{args.out_prefix}_taxa_metadata_{target_rank}.tsv"
    meta_df.to_csv(out_meta, sep="\t", index=False)

    print(f"\nSuccessfully created:", file=sys.stderr)
    print(f"  - {out_counts}", file=sys.stderr)
    if args.relative:
        print(f"  - {out_rel}", file=sys.stderr)
    print(f"  - {out_meta}", file=sys.stderr)
    print(f"\nFiltered to {len(abundance)} taxa across {len(sample_names)} samples.", file=sys.stderr)

if __name__ == "__main__":
    main()
