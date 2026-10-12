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

    # To also filter ancestors (up to genus) and descendants:
    ./filter_taxa.py --taxdump /path/to/taxdump \
                     --exclude exclude_taxids.txt \
                     --indir /path/to/per-read-files \
                     --outdir filtered_reads \
                     --include-ancestors \
                     --ancestors-max-rank genus \
                     --include-descendants
"""

import argparse
import re
import sys
import time
from collections import defaultdict, deque
from pathlib import Path

# Major ranks from lowest to highest; the ancestor walk never goes past --ancestors-max-rank
MAJOR_RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
# NCBI renamed 'superkingdom' to 'domain' in 2025; older taxdumps still use it
RANK_ALIASES = {"superkingdom": "domain"}
# root and 'cellular organisms': never removed, whatever the cap
ROOT_TAXIDS = {1, 131567}
# Kraken2 taxid column: '9606', or 'Homo sapiens (taxid 9606)' with --use-names
TAXID_FIELD = re.compile(r"(\d+)|.*\(taxid (\d+)\)")


def rank_level(rank: str) -> int:
    """Position of a rank in MAJOR_RANKS, or -1 for minor ranks and 'no rank'."""
    rank = RANK_ALIASES.get(rank, rank)
    return MAJOR_RANKS.index(rank) if rank in MAJOR_RANKS else -1


class TaxonomyTree:
    """Efficient NCBI taxonomy tree for ancestor/descendant computation."""

    def __init__(self, taxdump_dir: Path):
        """Initialize taxonomy tree from NCBI taxdump directory."""
        self.parent: dict[int, int] = {}  # taxid -> parent_taxid
        self.children: dict[int, list[int]] = defaultdict(
            list
        )  # taxid -> [child_taxids]
        self.rank: dict[int, str] = {}  # taxid -> rank name
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

    def get_ancestors(self, taxid: int, max_rank: str = "genus") -> set[int]:
        """
        Get ancestors of a taxid (including itself), walking up to max_rank.

        The node at max_rank is included; the walk stops before any major rank
        above it and never reaches the root, so e.g. reads assigned to Bacteria
        are never removed.
        """
        cap = rank_level(max_rank)
        ancestors = set()
        current = taxid

        while current in self.parent and current not in ancestors:
            level = rank_level(self.rank.get(current))
            if current != taxid and (current in ROOT_TAXIDS or level > cap):
                break
            ancestors.add(current)
            if level >= cap or self.parent[current] == current:
                break
            current = self.parent[current]

        return ancestors

    def get_descendants(self, taxid: int) -> set[int]:
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
        exclude_taxids: list[int],
        include_ancestors: bool = False,
        include_descendants: bool = False,
        ancestors_max_rank: str = "genus",
    ) -> set[int]:
        """
        Build forbidden taxa set from exclusion list.

        Args:
            exclude_taxids: List of taxids to exclude
            include_ancestors: If True, also exclude ancestors up to ancestors_max_rank
            include_descendants: If True, also exclude all descendants (walk down tree)
            ancestors_max_rank: Highest rank removed when include_ancestors is set

        Returns:
            Set of all forbidden taxids
        """
        forbidden = set()

        mode_desc = "exact taxids"
        if include_ancestors and include_descendants:
            mode_desc = f"exact taxids + ancestors up to {ancestors_max_rank} + descendants"
        elif include_ancestors:
            mode_desc = f"exact taxids + ancestors up to {ancestors_max_rank}"
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

            # Optionally include ancestors (walk up to the cap rank)
            ancestors = (
                self.get_ancestors(taxid, ancestors_max_rank) - {taxid}
                if include_ancestors
                else set()
            )
            to_add.update(ancestors)

            # Optionally include descendants (walk down)
            descendants = (
                self.get_descendants(taxid) - {taxid} if include_descendants else set()
            )
            to_add.update(descendants)

            forbidden.update(to_add)

            rank_name = self.rank.get(taxid, "unknown")

            if include_ancestors or include_descendants:
                parts = []
                if include_ancestors:
                    parts.append(f"+{len(ancestors)} ancestors")
                if include_descendants:
                    parts.append(f"+{len(descendants)} descendants")
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


def load_exclude_taxids(exclude_file: Path) -> list[int]:
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
    input_file: Path, output_file: Path, forbidden_taxa: set[int]
) -> tuple[int, int, int, int]:
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

            # Plain TaxID, or 'Homo sapiens (taxid 9606)' when Kraken2 ran with --use-names
            match = TAXID_FIELD.fullmatch(parts[2].strip())
            if not match:
                malformed_count += 1
                fout.write(line)
                retained_reads += 1
                continue
            taxid = int(match.group(1) or match.group(2))

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
    forbidden_taxa: set[int],
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
        total, removed, retained, _malformed = filter_classified_reads(
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

  # Also filter ancestors (up to genus) and descendants:
  %(prog)s --taxdump /path/to/taxdump \\
           --exclude exclude_taxids.txt \\
           --indir kraken2_reads \\
           --outdir filtered_reads \\
           --include-ancestors \\
           --ancestors-max-rank genus \\
           --include-descendants

The exclusion file should contain one taxid per line:
  9606     # Homo sapiens
  1751056  # Flavobacterium ammonificans
  297      # Hydrogenophilus thermoluteolus

Input files should be Kraken2 per-read output (tab-separated):
  C/U<TAB>readID<TAB>taxid<TAB>length<TAB>LCA k-mer mapping
  (the taxid column may also be "Name (taxid N)", as written by Kraken2 --use-names)

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
        help="Also exclude ancestor taxa, up to --ancestors-max-rank. Default: False",
    )

    parser.add_argument(
        "--ancestors-max-rank",
        choices=MAJOR_RANKS + ["superkingdom"],
        default="genus",
        help="Highest rank removed by --include-ancestors; the root and "
        "'cellular organisms' are never removed. Default: genus",
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
        ancestors_max_rank=args.ancestors_max_rank,
    )

    # Process all files
    process_directory(args.indir, args.outdir, forbidden_taxa, args.pattern)

    print("\nAll files processed successfully.", file=sys.stderr)


if __name__ == "__main__":
    main()
