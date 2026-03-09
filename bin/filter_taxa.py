#!/usr/bin/env python3
"""
filter_taxa.py - Filter Kraken2 per-read classification files by excluding taxa

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

This script removes reads assigned to specified taxa. By default, only the exact
taxids are filtered. Optionally, you can include ancestors and/or descendants.

Usage:
    ./filter_taxa.py --taxdump /path/to/taxdump \
                     --exclude exclude_taxids.txt \
                     --indir /path/to/per-read-files \
                     --outdir filtered_reads

    # To also filter ancestors and descendants:
    ./filter_taxa.py --taxdump /path/to/taxdump \
                     --exclude exclude_taxids.txt \
                     --indir /path/to/per-read-files \
                     --outdir filtered_reads \
                     --include-ancestors \
                     --include-descendants
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict, deque
from typing import Set, Dict, List, Tuple
import time


class TaxonomyTree:
    """Efficient NCBI taxonomy tree for ancestor/descendant computation."""

    def __init__(self, taxdump_dir: Path):
        """Initialize taxonomy tree from NCBI taxdump directory."""
        self.parent: Dict[int, int] = {}  # taxid -> parent_taxid
        self.children: Dict[int, List[int]] = defaultdict(
            list
        )  # taxid -> [child_taxids]
        self.rank: Dict[int, str] = {}  # taxid -> rank name
        self._load_nodes(taxdump_dir / "nodes.dmp")

    def _load_nodes(self, nodes_file: Path):
        """Parse nodes.dmp and build parent/child relationships."""
        if not nodes_file.exists():
            sys.exit(f"ERROR: nodes.dmp not found at {nodes_file}")

        print(f"Loading taxonomy from {nodes_file}...", file=sys.stderr)
        count = 0

        with open(nodes_file, "r") as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) < 3:
                    continue

                try:
                    taxid = int(parts[0])
                    parent_taxid = int(parts[1])
                    rank = parts[2]
                except (ValueError, IndexError):
                    continue

                self.parent[taxid] = parent_taxid
                self.rank[taxid] = rank

                # Build reverse mapping (parent -> children)
                if taxid != parent_taxid:  # Skip root (taxid 1)
                    self.children[parent_taxid].append(taxid)

                count += 1

        print(f"  Loaded {count:,} taxa from taxonomy", file=sys.stderr)

    def get_ancestors(self, taxid: int) -> Set[int]:
        """Get all ancestors of a taxid (including itself), walking up to root."""
        ancestors = set()
        current = taxid

        while current in self.parent and current not in ancestors:
            ancestors.add(current)
            parent = self.parent[current]
            if parent == current:  # Root node
                break
            current = parent

        return ancestors

    def get_descendants(self, taxid: int) -> Set[int]:
        """Get all descendants of a taxid (including itself) via BFS."""
        descendants = set()
        queue = deque([taxid])

        while queue:
            current = queue.popleft()
            if current in descendants:
                continue

            descendants.add(current)

            # Add all children to queue
            if current in self.children:
                queue.extend(self.children[current])

        return descendants

    def expand_forbidden_taxa(
        self,
        exclude_taxids: List[int],
        include_ancestors: bool = False,
        include_descendants: bool = False,
    ) -> Set[int]:
        """
        Build forbidden taxa set from exclusion list.

        Args:
            exclude_taxids: List of taxids to exclude
            include_ancestors: If True, also exclude all ancestors (walk up tree)
            include_descendants: If True, also exclude all descendants (walk down tree)

        Returns:
            Set of all forbidden taxids
        """
        forbidden = set()

        mode_desc = "exact taxids"
        if include_ancestors and include_descendants:
            mode_desc = "exact taxids + ancestors + descendants"
        elif include_ancestors:
            mode_desc = "exact taxids + ancestors"
        elif include_descendants:
            mode_desc = "exact taxids + descendants"

        print(
            f"\nProcessing {len(exclude_taxids)} excluded taxids ({mode_desc})...",
            file=sys.stderr,
        )

        for i, taxid in enumerate(exclude_taxids, 1):
            if taxid not in self.parent:
                print(
                    f"  WARNING: TaxID {taxid} not found in taxonomy (skipping)",
                    file=sys.stderr,
                )
                continue

            # Always include the exact taxid
            to_add = {taxid}

            # Optionally include ancestors (walk up)
            if include_ancestors:
                ancestors = self.get_ancestors(taxid)
                to_add.update(ancestors)

            # Optionally include descendants (walk down)
            if include_descendants:
                descendants = self.get_descendants(taxid)
                to_add.update(descendants)

            forbidden.update(to_add)

            rank_name = self.rank.get(taxid, "unknown")

            if include_ancestors or include_descendants:
                anc_count = len(self.get_ancestors(taxid)) if include_ancestors else 0
                desc_count = (
                    len(self.get_descendants(taxid)) if include_descendants else 0
                )
                parts = []
                if include_ancestors:
                    parts.append(f"+{anc_count} ancestors")
                if include_descendants:
                    parts.append(f"+{desc_count} descendants")
                extra = " (" + ", ".join(parts) + ")" if parts else ""
                print(
                    f"  [{i}/{len(exclude_taxids)}] TaxID {taxid} ({rank_name}){extra}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"  [{i}/{len(exclude_taxids)}] TaxID {taxid} ({rank_name})",
                    file=sys.stderr,
                )

        print(f"\nTotal forbidden taxa: {len(forbidden):,}", file=sys.stderr)
        return forbidden


def load_exclude_taxids(exclude_file: Path) -> List[int]:
    """Load taxids from exclusion file (one per line)."""
    if not exclude_file.exists():
        sys.exit(f"ERROR: Exclusion file not found: {exclude_file}")

    taxids = []
    with open(exclude_file, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Extract first field (handle lines with comments after taxid)
            taxid_str = line.split()[0]

            try:
                taxid = int(taxid_str)
                taxids.append(taxid)
            except ValueError:
                print(
                    f"WARNING: Invalid taxid on line {line_num}: {taxid_str}",
                    file=sys.stderr,
                )

    print(f"Loaded {len(taxids)} taxids from {exclude_file}", file=sys.stderr)
    return taxids


def filter_classified_reads(
    input_file: Path, output_file: Path, forbidden_taxa: Set[int]
) -> Tuple[int, int, int]:
    """
    Filter a Kraken2 classified reads file, removing forbidden taxa.

    Returns: (total_reads, removed_reads, retained_reads, malformed_count)
    """
    total_reads = 0
    removed_reads = 0
    retained_reads = 0
    malformed_count = 0

    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            total_reads += 1

            # Parse: C/U\treadID\ttaxid\tlengths\tk-mer_mappings
            parts = line.strip().split("\t")
            if len(parts) < 3:
                malformed_count += 1
                fout.write(line)
                retained_reads += 1
                continue

            try:
                taxid = int(parts[2])
            except (ValueError, IndexError):
                malformed_count += 1
                fout.write(line)
                retained_reads += 1
                continue

            # Filter out forbidden taxa
            if taxid in forbidden_taxa:
                removed_reads += 1
            else:
                fout.write(line)
                retained_reads += 1

    if malformed_count > 0:
        print(
            f"  WARNING: {malformed_count} malformed/invalid lines in {input_file.name}",
            file=sys.stderr,
        )

    return total_reads, removed_reads, retained_reads, malformed_count


def process_directory(
    input_dir: Path,
    output_dir: Path,
    forbidden_taxa: Set[int],
    pattern: str = "*.classifiedreads.txt",
):
    """Process all classified reads files in a directory."""

    # Find all matching files
    input_files = sorted(input_dir.glob(pattern))

    if not input_files:
        sys.exit(f"ERROR: No files matching pattern '{pattern}' found in {input_dir}")

    print(f"\nFound {len(input_files)} files to process", file=sys.stderr)
    print(f"Output directory: {output_dir}\n", file=sys.stderr)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Summary statistics
    grand_total = 0
    grand_removed = 0
    grand_retained = 0

    # Process each file
    for i, input_file in enumerate(input_files, 1):
        start_time = time.time()

        # Generate output filename: sample.classifiedreads.txt -> sample.classifiedreads.filtered.txt
        output_name = input_file.stem + ".filtered" + input_file.suffix
        output_file = output_dir / output_name

        print(
            f"[{i}/{len(input_files)}] Processing {input_file.name}...", file=sys.stderr
        )

        # Filter the file
        total, removed, retained, malformed = filter_classified_reads(
            input_file, output_file, forbidden_taxa
        )

        elapsed = time.time() - start_time
        removal_pct = (removed / total * 100) if total > 0 else 0

        print(f"  Total reads:    {total:>12,}", file=sys.stderr)
        print(
            f"  Removed:        {removed:>12,} ({removal_pct:>5.2f}%)", file=sys.stderr
        )
        print(f"  Retained:       {retained:>12,}", file=sys.stderr)
        print(f"  Time:           {elapsed:>12.2f}s", file=sys.stderr)
        print(f"  Output:         {output_file.name}\n", file=sys.stderr)

        grand_total += total
        grand_removed += removed
        grand_retained += retained

    # Print summary
    print("=" * 70, file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print("=" * 70, file=sys.stderr)
    print(f"Files processed:     {len(input_files)}", file=sys.stderr)
    print(f"Total reads:         {grand_total:,}", file=sys.stderr)
    removal_pct = (grand_removed / grand_total * 100) if grand_total > 0 else 0.0
    print(
        f"Total removed:       {grand_removed:,} ({removal_pct:.2f}%)",
        file=sys.stderr,
    )
    print(f"Total retained:      {grand_retained:,}", file=sys.stderr)
    print(f"Forbidden taxa:      {len(forbidden_taxa):,}", file=sys.stderr)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter Kraken2 classified reads by excluding taxa (exact taxids by default)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Filter only exact taxids (default):
  %(prog)s --taxdump /path/to/taxdump \\
           --exclude exclude_taxids.txt \\
           --indir kraken2_reads \\
           --outdir filtered_reads

  # Also filter ancestors and descendants:
  %(prog)s --taxdump /path/to/taxdump \\
           --exclude exclude_taxids.txt \\
           --indir kraken2_reads \\
           --outdir filtered_reads \\
           --include-ancestors \\
           --include-descendants

The exclusion file should contain one taxid per line:
  9606     # Homo sapiens
  1751056  # Flavobacterium ammonificans
  297      # Hydrogenophilus thermoluteolus

Input files should be Kraken2 per-read format (tab-separated):
  readID<TAB>taxid<TAB>score

Output files will have '.filtered' appended before the extension.
        """,
    )

    parser.add_argument(
        "--taxdump",
        type=Path,
        required=True,
        help="Path to NCBI taxdump directory containing nodes.dmp and names.dmp",
    )

    parser.add_argument(
        "--exclude",
        type=Path,
        required=True,
        help="File containing taxids to exclude (one per line, # for comments)",
    )

    parser.add_argument(
        "--indir",
        type=Path,
        required=True,
        help="Directory containing Kraken2 classified reads files",
    )

    parser.add_argument(
        "--outdir", type=Path, required=True, help="Output directory for filtered files"
    )

    parser.add_argument(
        "--pattern",
        default="*.classifiedreads.txt",
        help="Glob pattern for input files (default: *.classifiedreads.txt)",
    )

    parser.add_argument(
        "--include-ancestors",
        action="store_true",
        help="Also exclude all ancestor taxa (walking up to root). Default: False",
    )

    parser.add_argument(
        "--include-descendants",
        action="store_true",
        help="Also exclude all descendant taxa (walking down the tree). Default: False",
    )

    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()

    print("=" * 70, file=sys.stderr)
    print("LowBioPipe - Filter Taxa Module", file=sys.stderr)
    print("=" * 70, file=sys.stderr)

    # Validate inputs
    if not args.taxdump.is_dir():
        sys.exit(f"ERROR: Taxdump directory not found: {args.taxdump}")

    if not args.indir.is_dir():
        sys.exit(f"ERROR: Input directory not found: {args.indir}")

    # Load taxonomy
    start_time = time.time()
    taxonomy = TaxonomyTree(args.taxdump)
    print(f"  Loaded in {time.time() - start_time:.2f}s\n", file=sys.stderr)

    # Load exclusion list
    exclude_taxids = load_exclude_taxids(args.exclude)

    # Build forbidden set (with optional ancestors/descendants)
    forbidden_taxa = taxonomy.expand_forbidden_taxa(
        exclude_taxids,
        include_ancestors=args.include_ancestors,
        include_descendants=args.include_descendants,
    )

    # Process all files
    process_directory(args.indir, args.outdir, forbidden_taxa, args.pattern)

    print("\nAll files processed successfully.", file=sys.stderr)


if __name__ == "__main__":
    main()
