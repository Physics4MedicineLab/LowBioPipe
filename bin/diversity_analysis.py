#!/usr/bin/env python3
"""
diversity_analysis.py - Comprehensive diversity analysis for microbiome data

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Features:
    - Automatic detection of all taxonomic levels
    - Alpha diversity metrics (observed OTUs, Shannon, Simpson, Chao1)
    - Beta diversity (Bray-Curtis, Jaccard, Aitchison)
    - PCoA ordination plots
    - Heatmaps with hierarchical clustering
    - PERMANOVA statistical tests (if groups provided)

Usage:
    python diversity_analysis.py --indir abundance_dir --outdir results
    python diversity_analysis.py --indir . --outdir results --groups groups.tsv
"""

import argparse
import logging
import sys
from typing import Optional, Dict, Tuple, List
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="[%(levelname)s] %(message)s", stream=sys.stderr
)
logger = logging.getLogger(__name__)

# Valid taxonomic ranks for validation
VALID_RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

# ======================= helpers =======================


def ensure_samples_by_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Transpose DataFrame to samples × features format.

    Args:
        df: DataFrame in features × samples format

    Returns:
        Transposed DataFrame with samples as rows and features as columns
    """
    return df.T.copy()


def drop_empty_samples(df_sxf: pd.DataFrame) -> pd.DataFrame:
    """
    Remove samples with zero total counts.

    Args:
        df_sxf: DataFrame in samples × features format

    Returns:
        DataFrame with empty samples removed
    """
    empty = df_sxf.sum(axis=1) == 0
    if empty.any():
        dropped = list(df_sxf.index[empty])
        logger.warning("Dropping empty samples: %s", dropped)
        df_sxf = df_sxf.loc[~empty]
    return df_sxf


def validate_dataframe(
    df: pd.DataFrame, min_samples: int = 2, min_features: int = 1
) -> bool:
    """
    Validate DataFrame has sufficient data for analysis.

    Args:
        df: DataFrame to validate
        min_samples: Minimum number of samples required
        min_features: Minimum number of features required

    Returns:
        True if valid, False otherwise
    """
    if df.empty:
        logger.error("DataFrame is empty")
        return False
    if df.shape[0] < min_samples:
        logger.error(
            "Insufficient samples: %d (minimum %d required)", df.shape[0], min_samples
        )
        return False
    if df.shape[1] < min_features:
        logger.error(
            "Insufficient features: %d (minimum %d required)", df.shape[1], min_features
        )
        return False
    return True


def compute_alpha(df_sxf: pd.DataFrame) -> pd.DataFrame:
    """
    Compute alpha diversity metrics for all samples.

    Args:
        df_sxf: DataFrame in samples × features format (counts)

    Returns:
        DataFrame with alpha diversity metrics (observed_otus, shannon, simpson, chao1)
    """
    arr = df_sxf.values
    ids = df_sxf.index

    # Compute observed OTUs with fallback
    try:
        observed = alpha_diversity("observed_otus", arr, ids=ids)
        observed = pd.Series(observed, index=ids, name="observed_otus")
    except ValueError as e:
        logger.warning("observed_otus metric failed (%s), computing manually", e)
        observed = pd.Series((arr > 0).sum(axis=1), index=ids, name="observed_otus")

    # Compute other metrics
    metrics = {}
    for m in ["shannon", "simpson", "chao1"]:
        try:
            metrics[m] = alpha_diversity(m, arr, ids=ids)
        except ValueError as e:
            logger.warning("Alpha metric '%s' failed: %s", m, e)
            metrics[m] = pd.Series([np.nan] * len(ids), index=ids)

    alpha_df = pd.DataFrame(metrics, index=ids)
    alpha_df.insert(0, "observed_otus", observed)
    return alpha_df


def save_alpha_plots(alpha_df: pd.DataFrame, outdir: Path, level: str) -> None:
    """
    Generate and save alpha diversity visualization plots.

    Args:
        alpha_df: DataFrame with alpha diversity metrics
        outdir: Output directory path
        level: Taxonomic level name for plot titles
    """
    sns.set_theme(style="whitegrid")
    for metric in alpha_df.columns:
        # Boxplot
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=alpha_df[metric], orient="h")
        sns.stripplot(data=alpha_df[metric], orient="h", color="black", size=4)
        plt.title(f"Alpha Diversity — {metric} — {level}")
        plt.xlabel(metric)
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_boxplot_{level}.png", dpi=300)
        plt.close()

        # Violin plot
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=alpha_df[metric], orient="h", inner="box")
        sns.stripplot(data=alpha_df[metric], orient="h", color="black", size=4)
        plt.title(f"Alpha Diversity — {metric} (violin) — {level}")
        plt.xlabel(metric)
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_violin_{level}.png", dpi=300)
        plt.close()


def compute_beta_and_pcoa(
    df_sxf: pd.DataFrame, metric: str, binarize: bool = False
) -> Tuple[DistanceMatrix, pd.DataFrame, pd.Series]:
    """
    Compute beta diversity distance matrix and perform PCoA ordination.

    Args:
        df_sxf: DataFrame in samples × features format
        metric: Distance metric name (e.g., 'braycurtis', 'jaccard')
        binarize: If True, convert counts to presence/absence before calculation

    Returns:
        Tuple of (distance matrix, PCoA coordinates, proportion explained)
    """
    X = df_sxf.values
    if binarize:
        X = (X > 0).astype(int)
    ids = df_sxf.index
    dist = beta_diversity(metric, X, ids=ids)
    ord_res = pcoa(dist)
    coords = ord_res.samples.copy()
    exp = (ord_res.proportion_explained * 100).round(2)
    return dist, coords, exp


def clr_transform(X: np.ndarray) -> np.ndarray:
    """
    Apply Centered Log-Ratio (CLR) transformation for compositional data.

    The CLR transformation is defined as:
        clr(x) = log(x / geometric_mean(x))

    Args:
        X: 2D array of counts (samples × features)

    Returns:
        CLR-transformed array
    """
    # Add pseudocount to handle zeros
    X_pseudo = X + 1e-10
    log_X = np.log(X_pseudo)
    geometric_mean = log_X.mean(axis=1, keepdims=True)
    return log_X - geometric_mean


def compute_aitchison_and_pcoa(
    df_sxf: pd.DataFrame,
) -> Tuple[DistanceMatrix, pd.DataFrame, pd.Series]:
    """
    Compute Aitchison distance (Euclidean distance in CLR space) and PCoA.

    Aitchison distance is appropriate for compositional data analysis
    where samples represent proportions that sum to a constant.

    Args:
        df_sxf: DataFrame in samples × features format

    Returns:
        Tuple of (distance matrix, PCoA coordinates, proportion explained)
    """
    X = df_sxf.values
    ids = list(df_sxf.index)

    # CLR transform and compute Euclidean distance
    clr_data = clr_transform(X)
    dist_array = pdist(clr_data, metric="euclidean")
    dist_matrix = squareform(dist_array)

    # Create skbio DistanceMatrix for compatibility
    dist = DistanceMatrix(dist_matrix, ids=ids)

    # PCoA
    ord_res = pcoa(dist)
    coords = ord_res.samples.copy()
    exp = (ord_res.proportion_explained * 100).round(2)
    return dist, coords, exp


def plot_pcoa(coords: pd.DataFrame, exp: pd.Series, out_png: Path, title: str) -> None:
    """
    Generate and save PCoA ordination plot.

    Args:
        coords: DataFrame with PCoA coordinates (PC1, PC2, ...)
        exp: Series with proportion of variance explained per axis
        out_png: Output file path
        title: Plot title
    """
    plt.figure(figsize=(7.5, 6.5))
    plt.scatter(coords["PC1"], coords["PC2"], s=80, alpha=0.9)
    plt.xlabel(f"PC1 ({exp.iloc[0]:.2f}%)")
    plt.ylabel(f"PC2 ({exp.iloc[1]:.2f}%)")
    plt.title(f"PCoA — {title}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def make_heatmap(df_sxf: pd.DataFrame, out_png: Path, title: str) -> None:
    """
    Generate hierarchical clustering heatmap.

    Args:
        df_sxf: DataFrame in samples × features format
        out_png: Output file path
        title: Plot title
    """
    g = sns.clustermap(df_sxf, cmap="viridis", figsize=(12, 12))
    g.fig.suptitle(title, y=1.02)
    g.savefig(out_png, dpi=300)
    plt.close(g.fig)


def run_permanova(
    dist_matrix: DistanceMatrix, sample_groups: List[str], permutations: int = 999
) -> object:
    """
    Run PERMANOVA statistical test on distance matrix.

    Args:
        dist_matrix: skbio DistanceMatrix object
        sample_groups: List of group labels for each sample
        permutations: Number of permutations for p-value calculation

    Returns:
        PERMANOVA results object
    """
    df = pd.DataFrame({"group": sample_groups}, index=dist_matrix.ids)
    result = permanova(dist_matrix, df, column="group", permutations=permutations)
    return result


def save_beta_results(
    dist: DistanceMatrix,
    coords: pd.DataFrame,
    exp: pd.Series,
    outdir: Path,
    level: str,
    metric_name: str,
    groups: Optional[Dict[str, str]] = None,
) -> None:
    """
    Save beta diversity results (distance matrix, PCoA, optional PERMANOVA).

    Args:
        dist: Distance matrix
        coords: PCoA coordinates
        exp: Proportion explained per axis
        outdir: Output directory
        level: Taxonomic level
        metric_name: Name of the distance metric
        groups: Optional sample-to-group mapping for PERMANOVA
    """
    # Save distance matrix
    pd.DataFrame(dist.data, index=dist.ids, columns=dist.ids).to_csv(
        outdir / f"beta_{metric_name}_{level}.tsv", sep="\t"
    )

    # Save PCoA coordinates
    coords.to_csv(outdir / f"pcoa_{metric_name}_coords_{level}.tsv", sep="\t")

    # Save PCoA plot
    plot_pcoa(
        coords,
        exp,
        outdir / f"pcoa_{metric_name}_{level}.png",
        f"{level} • {metric_name}",
    )

    # Run PERMANOVA if groups provided
    if groups:
        sample_groups = [groups.get(s, "ungrouped") for s in dist.ids]
        result = run_permanova(dist, sample_groups)
        with open(outdir / f"permanova_{metric_name}_{level}.txt", "w") as f:
            f.write(str(result))


def infer_level_from_path(p: Path) -> str:
    """
    Infer taxonomic level from file path.

    Args:
        p: Path to abundance table

    Returns:
        Inferred taxonomic level or filename stem if not recognized
    """
    stem = p.stem.lower()
    for k in VALID_RANKS:
        if k in stem:
            return k
    return stem


# ======================= pipeline =======================


def process_table(
    path: Path, outdir_base: Path, groups: Optional[Dict[str, str]]
) -> None:
    """
    Process a single abundance table through complete diversity analysis.

    Args:
        path: Path to abundance table (TSV format)
        outdir_base: Base output directory
        groups: Optional sample-to-group mapping for PERMANOVA
    """
    logger.info("Processing: %s", path)
    level = infer_level_from_path(path)
    outdir = outdir_base / level
    outdir.mkdir(parents=True, exist_ok=True)

    # Read and validate data
    try:
        df = pd.read_csv(path, sep="\t", index_col=0)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        logger.error("Failed to read %s: %s", path, e)
        return

    df_sxf = ensure_samples_by_features(df)
    df_sxf = drop_empty_samples(df_sxf)

    if not validate_dataframe(df_sxf, min_samples=2):
        logger.error("Skipping %s due to insufficient data", path)
        return

    # ---- Alpha ----
    alpha_df = compute_alpha(df_sxf)
    alpha_df.to_csv(outdir / f"alpha_metrics_{level}.tsv", sep="\t")
    save_alpha_plots(alpha_df, outdir, level)

    # ---- Heatmap + clustering ----
    make_heatmap(df_sxf, outdir / f"heatmap_{level}.png", f"Heatmap {level}")

    # ---- Beta Bray-Curtis ----
    bray_dist, bray_coords, bray_exp = compute_beta_and_pcoa(df_sxf, "braycurtis")
    save_beta_results(
        bray_dist, bray_coords, bray_exp, outdir, level, "braycurtis", groups
    )

    if groups:
        logger.info("PERMANOVA (Bray-Curtis) done")

    # ---- Beta Jaccard ----
    jacc_dist, jacc_coords, jacc_exp = compute_beta_and_pcoa(
        df_sxf, "jaccard", binarize=True
    )
    save_beta_results(jacc_dist, jacc_coords, jacc_exp, outdir, level, "jaccard")

    # ---- Beta Aitchison ----
    aitc_dist, aitc_coords, aitc_exp = compute_aitchison_and_pcoa(df_sxf)
    save_beta_results(
        aitc_dist, aitc_coords, aitc_exp, outdir, level, "aitchison", groups
    )

    logger.info("Done: %s saved to %s", level, outdir)


# ======================= CLI =======================


def main() -> None:
    """Main entry point for diversity analysis."""
    ap = argparse.ArgumentParser(description="Full diversity analysis pipeline")
    ap.add_argument(
        "--indir", default=".", help="Directory containing abundance tables"
    )
    ap.add_argument("--outdir", default="diversity_outputs", help="Output directory")
    ap.add_argument(
        "--groups", help="Optional TSV mapping samples to groups for PERMANOVA"
    )
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir_base = Path(args.outdir)
    outdir_base.mkdir(parents=True, exist_ok=True)

    # Validate input directory
    if not indir.exists():
        logger.error("Input directory does not exist: %s", indir)
        sys.exit(1)

    # Load groups file if provided
    groups = None
    if args.groups:
        groups_path = Path(args.groups)
        if not groups_path.exists():
            logger.error("Groups file does not exist: %s", groups_path)
            sys.exit(1)
        try:
            gdf = pd.read_csv(
                args.groups, sep="\t", header=None, names=["sample", "group"]
            )
            groups = dict(zip(gdf.sample, gdf.group))
            logger.info("Groups loaded for PERMANOVA (%d samples)", len(groups))
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error("Failed to read groups file: %s", e)
            sys.exit(1)

    # Auto-detect all *_counts_*.tsv files
    tables = sorted(indir.glob("*counts_*.tsv"))
    if not tables:
        logger.error("No abundance tables found in %s", indir)
        sys.exit(1)

    logger.info("Found %d abundance tables", len(tables))
    for table in tables:
        process_table(table, outdir_base, groups)


if __name__ == "__main__":
    sns.set_theme(style="whitegrid")
    main()
