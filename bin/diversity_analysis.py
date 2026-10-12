#!/usr/bin/env python3
"""
diversity_analysis.py - Comprehensive diversity analysis for microbiome data

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Features:
    - Automatic detection of all taxonomic levels
    - Alpha diversity metrics (observed OTUs, Shannon, Simpson, Chao1)
    - Beta diversity (Bray-Curtis, Jaccard, Aitchison)
    - PCoA ordination plots (coloured by group if groups provided)
    - Heatmaps with hierarchical clustering of the most abundant taxa
    - PERMANOVA and PERMDISP statistical tests (if groups provided)

Usage:
    python diversity_analysis.py --indir abundance_dir --outdir results
    python diversity_analysis.py --indir . --outdir results --groups groups.tsv --heatmap-top-n 50
"""

import argparse
import inspect
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import permanova, permdisp
from skbio.stats.ordination import pcoa

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="[%(levelname)s] %(message)s", stream=sys.stderr
)
logger = logging.getLogger(__name__)


def ensure_samples_by_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Transpose DataFrame to samples x features format.

    Args:
        df: DataFrame in features x samples format

    Returns:
        Transposed DataFrame with samples as rows and features as columns
    """
    return df.T.copy()


def drop_empty_samples(df_sxf: pd.DataFrame) -> pd.DataFrame:
    """
    Remove samples with zero total counts.

    Args:
        df_sxf: DataFrame in samples x features format

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
        df_sxf: DataFrame in samples x features format (counts)

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


def group_labels(
    sample_ids: pd.Index, groups: dict[str, str] | None
) -> pd.Series | None:
    """
    Group label for each sample, for plotting.

    Args:
        sample_ids: Sample IDs
        groups: Optional sample-to-group mapping

    Returns:
        Series of group labels ('ungrouped' if missing), or None without groups
    """
    if not groups:
        return None
    return pd.Series(
        [groups.get(s, "ungrouped") for s in sample_ids], index=sample_ids, name="group"
    )


def save_alpha_plots(
    alpha_df: pd.DataFrame,
    outdir: Path,
    level: str,
    groups: dict[str, str] | None = None,
) -> None:
    """
    Generate and save alpha diversity box and violin plots, one box per group.

    Args:
        alpha_df: DataFrame with alpha diversity metrics
        outdir: Output directory path
        level: Taxonomic level name for plot titles
        groups: Optional sample-to-group mapping
    """
    sns.set_theme(style="whitegrid")
    labels = group_labels(alpha_df.index, groups)
    y = "group" if labels is not None else None
    plots = [
        ("boxplot", "", sns.boxplot, {}),
        ("violin", " (violin)", sns.violinplot, {"inner": "box"}),
    ]
    for metric in alpha_df.columns:
        plot_df = alpha_df[[metric]].assign(group=labels) if y else alpha_df[[metric]]
        for kind, suffix, plot_fn, kwargs in plots:
            try:
                plt.figure(figsize=(10, 6))
                plot_fn(data=plot_df, x=metric, y=y, orient="h", **kwargs)
                sns.stripplot(
                    data=plot_df, x=metric, y=y, orient="h", color="black", size=4
                )
                plt.title(f"Alpha Diversity - {metric}{suffix} - {level}")
                plt.xlabel(metric)
                plt.tight_layout()
                plt.savefig(outdir / f"alpha_{metric}_{kind}_{level}.png", dpi=300)
            except Exception:
                logger.exception("Failed to save alpha %s for %s/%s", kind, metric, level)
            finally:
                plt.close()


def compute_beta_and_pcoa(
    df_sxf: pd.DataFrame, metric: str, binarize: bool = False
) -> tuple[DistanceMatrix, pd.DataFrame, pd.Series]:
    """
    Compute beta diversity distance matrix and perform PCoA ordination.

    Args:
        df_sxf: DataFrame in samples x features format
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
        X: 2D array of counts (samples x features)

    Returns:
        CLR-transformed array
    """
    # Add pseudocount to handle zeros (Bayesian uniform prior;
    # Martin-Fernandez et al. 2003, doi:10.1023/A:1023866030544)
    X_pseudo = X + 0.5
    log_X = np.log(X_pseudo)
    geometric_mean = log_X.mean(axis=1, keepdims=True)
    return log_X - geometric_mean


def compute_aitchison_and_pcoa(
    df_sxf: pd.DataFrame,
) -> tuple[DistanceMatrix, pd.DataFrame, pd.Series]:
    """
    Compute Aitchison distance (Euclidean distance in CLR space) and PCoA.

    Aitchison distance is appropriate for compositional data analysis
    where samples represent proportions that sum to a constant.

    Args:
        df_sxf: DataFrame in samples x features format

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


def plot_pcoa(
    coords: pd.DataFrame,
    exp: pd.Series,
    out_png: Path,
    title: str,
    labels: pd.Series | None = None,
) -> None:
    """
    Generate and save PCoA ordination plot with sample IDs next to the points.

    Args:
        coords: DataFrame with PCoA coordinates (PC1, PC2, ...)
        exp: Series with proportion of variance explained per axis
        out_png: Output file path
        title: Plot title
        labels: Optional group label per sample, used to colour the points
    """
    try:
        plt.figure(figsize=(7.5, 6.5))
        hue = labels.reindex(coords.index) if labels is not None else None
        sns.scatterplot(x=coords["PC1"], y=coords["PC2"], hue=hue, s=80, alpha=0.9)
        for sample, (x, y) in coords[["PC1", "PC2"]].iterrows():
            plt.annotate(
                str(sample), (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7
            )
        plt.xlabel(f"PC1 ({exp.iloc[0]:.2f}%)")
        plt.ylabel(f"PC2 ({exp.iloc[1]:.2f}%)")
        plt.title(f"PCoA - {title}")
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
    except Exception:
        logger.exception("Failed to save PCoA plot %s", out_png)
    finally:
        plt.close()


def make_heatmap(
    df_sxf: pd.DataFrame, out_png: Path, title: str, top_n: int = 50
) -> None:
    """
    Generate hierarchical clustering heatmap of log1p counts.

    Args:
        df_sxf: DataFrame in samples x features format (counts)
        out_png: Output file path
        title: Plot title
        top_n: Keep the top_n taxa by mean relative abundance (0 = all taxa)
    """
    if top_n and df_sxf.shape[1] > top_n:
        rel = df_sxf.div(df_sxf.sum(axis=1), axis=0)
        df_sxf = df_sxf[rel.mean().nlargest(top_n).index]
        title = f"{title}, top {top_n} taxa"
    try:
        g = sns.clustermap(np.log1p(df_sxf), cmap="viridis", figsize=(12, 12))
        g.fig.suptitle(title, y=1.02)
        g.savefig(out_png, dpi=300)
    except Exception:
        logger.exception("Failed to save heatmap %s", out_png)
    finally:
        plt.close("all")


def run_group_test(
    test,
    dist_matrix: DistanceMatrix,
    sample_groups: list[str],
    permutations: int = 999,
    seed: int = 42,
    **kwargs,
) -> pd.Series:
    """
    Run a scikit-bio grouping test (permanova or permdisp) on a distance matrix.

    Args:
        test: skbio.stats.distance.permanova or skbio.stats.distance.permdisp
        dist_matrix: skbio DistanceMatrix object
        sample_groups: List of group labels for each sample
        permutations: Number of permutations for p-value calculation
        seed: Random seed for the permutations
        **kwargs: Extra arguments for the test

    Returns:
        Test results
    """
    df = pd.DataFrame({"group": sample_groups}, index=dist_matrix.ids)
    # scikit-bio >= 0.6 ignores the global NumPy seed and takes its own
    if "seed" in inspect.signature(test).parameters:
        kwargs["seed"] = seed
    else:
        np.random.seed(seed)
    return test(dist_matrix, df, column="group", permutations=permutations, **kwargs)


def save_beta_results(
    dist: DistanceMatrix,
    coords: pd.DataFrame,
    exp: pd.Series,
    outdir: Path,
    level: str,
    metric_name: str,
    groups: dict[str, str] | None = None,
    seed: int = 42,
) -> None:
    """
    Save beta diversity results (distance matrix, PCoA, optional PERMANOVA/PERMDISP).

    Args:
        dist: Distance matrix
        coords: PCoA coordinates
        exp: Proportion explained per axis
        outdir: Output directory
        level: Taxonomic level
        metric_name: Name of the distance metric
        groups: Optional sample-to-group mapping for plots, PERMANOVA and PERMDISP
        seed: Random seed for the permutation tests
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
        f"{level} - {metric_name}",
        group_labels(coords.index, groups),
    )

    # Group tests, on the samples listed in the groups file only
    if groups:
        grouped = [s for s in dist.ids if s in groups]
        ungrouped = [s for s in dist.ids if s not in groups]
        if ungrouped:
            logger.warning(
                "PERMANOVA/PERMDISP %s/%s: samples not in groups file are excluded: %s",
                level,
                metric_name,
                ungrouped,
            )
        sample_groups = [groups[s] for s in grouped]
        group_counts = pd.Series(sample_groups).value_counts()
        if len(group_counts) < 2 or len(group_counts) == len(grouped):
            logger.warning(
                "PERMANOVA/PERMDISP %s/%s skipped: needs >= 2 groups and at least one "
                "group with >= 2 samples (%d samples in %d groups)",
                level,
                metric_name,
                len(grouped),
                len(group_counts),
            )
            return
        small_groups = group_counts[group_counts < 3]
        if not small_groups.empty:
            logger.warning(
                "PERMANOVA/PERMDISP %s/%s: groups with <3 samples may yield unreliable "
                "p-values: %s",
                level,
                metric_name,
                dict(small_groups),
            )
        dist_grouped = dist.filter(grouped)
        # PERMDISP tests equal dispersion: a significant PERMANOVA together with a
        # significant PERMDISP may reflect spread rather than location differences.
        # It uses all PCoA axes (scikit-bio >= 0.7.4 otherwise truncates to 10).
        tests = (
            ("permanova", permanova, {}),
            ("permdisp", permdisp, {"dimensions": len(grouped)}),
        )
        for name, test, kwargs in tests:
            try:
                result = run_group_test(
                    test, dist_grouped, sample_groups, seed=seed, **kwargs
                )
            # ArithmeticError: scikit-bio < 0.7.4 permdisp can divide by zero on tied samples
            except (ValueError, ArithmeticError) as e:
                logger.warning("%s %s/%s failed: %s", name, level, metric_name, e)
                continue
            with open(outdir / f"{name}_{metric_name}_{level}.txt", "w") as f:
                f.write(str(result))


def infer_level_from_path(p: Path) -> str:
    """
    Infer taxonomic level from an abundance table name (<prefix>_counts_<rank>.tsv).

    Args:
        p: Path to abundance table

    Returns:
        Taxonomic level after '_counts_', or the filename stem if absent
    """
    return p.stem.lower().rsplit("_counts_", 1)[-1]


def process_table(
    path: Path,
    outdir_base: Path,
    groups: dict[str, str] | None,
    seed: int = 42,
    heatmap_top_n: int = 50,
) -> None:
    """
    Process a single abundance table through complete diversity analysis.

    Args:
        path: Path to abundance table (TSV format)
        outdir_base: Base output directory
        groups: Optional sample-to-group mapping for plots, PERMANOVA and PERMDISP
        seed: Random seed for the permutation tests
        heatmap_top_n: Number of most abundant taxa in the heatmap (0 = all)
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
    save_alpha_plots(alpha_df, outdir, level, groups)

    # ---- Heatmap + clustering (log1p-transformed for visual balance) ----
    make_heatmap(
        df_sxf,
        outdir / f"heatmap_{level}.png",
        f"Heatmap {level} (log1p)",
        heatmap_top_n,
    )

    # ---- Beta Bray-Curtis ----
    bray_dist, bray_coords, bray_exp = compute_beta_and_pcoa(df_sxf, "braycurtis")
    save_beta_results(
        bray_dist, bray_coords, bray_exp, outdir, level, "braycurtis", groups, seed
    )

    if groups:
        logger.info("PERMANOVA (Bray-Curtis) done")

    # ---- Beta Jaccard ----
    jacc_dist, jacc_coords, jacc_exp = compute_beta_and_pcoa(
        df_sxf, "jaccard", binarize=True
    )
    save_beta_results(
        jacc_dist, jacc_coords, jacc_exp, outdir, level, "jaccard", groups, seed
    )

    # ---- Beta Aitchison ----
    aitc_dist, aitc_coords, aitc_exp = compute_aitchison_and_pcoa(df_sxf)
    save_beta_results(
        aitc_dist, aitc_coords, aitc_exp, outdir, level, "aitchison", groups, seed
    )

    if groups:
        logger.warning(
            "Note: PERMANOVA/PERMDISP p-values for %s are not corrected for multiple "
            "comparisons (3 metrics tested). Consider Benjamini-Hochberg FDR "
            "correction when interpreting results.",
            level,
        )

    logger.info("Done: %s saved to %s", level, outdir)


def main() -> None:
    """Main entry point for diversity analysis."""
    ap = argparse.ArgumentParser(description="Full diversity analysis pipeline")
    ap.add_argument(
        "--indir", default=".", help="Directory containing abundance tables"
    )
    ap.add_argument("--outdir", default="diversity_outputs", help="Output directory")
    ap.add_argument(
        "--groups",
        help="Optional TSV mapping samples to groups for PERMANOVA, PERMDISP and plots",
    )
    ap.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducible PERMANOVA/PERMDISP permutations (default: 42)",
    )
    ap.add_argument(
        "--heatmap-top-n", type=int, default=50,
        help="Most abundant taxa (by mean relative abundance) in the heatmap; 0 = all (default: 50)",
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
            # dtype=str keeps numeric IDs ('01') matching the abundance-table headers;
            # gdf["sample"] because gdf.sample is the DataFrame.sample() method
            gdf = pd.read_csv(
                args.groups, sep="\t", header=None, names=["sample", "group"], dtype=str
            )
            groups = dict(zip(gdf["sample"].str.strip(), gdf["group"].str.strip()))
            logger.info("Groups loaded (%d samples)", len(groups))
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
            logger.error("Failed to read groups file: %s", e)
            sys.exit(1)

    # Auto-detect all *_counts_*.tsv files
    tables = sorted(indir.glob("*counts_*.tsv"))
    if not tables:
        logger.error("No abundance tables found in %s", indir)
        sys.exit(1)

    logger.info("Random seed set to %d for reproducible PERMANOVA/PERMDISP", args.seed)

    logger.info("Found %d abundance tables", len(tables))
    for table in tables:
        process_table(table, outdir_base, groups, args.seed, args.heatmap_top_n)


if __name__ == "__main__":
    sns.set_theme(style="whitegrid")
    main()
