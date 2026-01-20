#!/usr/bin/env python3
"""
diversity_analysis.py - Comprehensive diversity analysis for microbiome data

Part of LowBioPipe - A pipeline for low biomass microbiome analysis

Features:
    - Automatic detection of all taxonomic levels
    - Alpha diversity metrics (observed OTUs, Shannon, Simpson, Chao1)
    - Beta diversity (Bray-Curtis, Jaccard)
    - PCoA ordination plots
    - Heatmaps with hierarchical clustering
    - PERMANOVA statistical tests (if groups provided)

Usage:
    python diversity_analysis.py --indir abundance_dir --outdir results
    python diversity_analysis.py --indir . --outdir results --groups groups.tsv
"""

import argparse
from typing import Optional, Dict, Tuple
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova

# ======================= helpers =======================

def ensure_samples_by_features(df: pd.DataFrame) -> pd.DataFrame:
    return df.T.copy()

def drop_empty_samples(df_sxf: pd.DataFrame) -> pd.DataFrame:
    empty = (df_sxf.sum(axis=1) == 0)
    if empty.any():
        print(f"[WARN] Dropping empty samples: {list(df_sxf.index[empty])}")
        df_sxf = df_sxf.loc[~empty]
    return df_sxf

def compute_alpha(df_sxf: pd.DataFrame) -> pd.DataFrame:
    arr = df_sxf.values
    ids = df_sxf.index
    try:
        observed = alpha_diversity("observed_otus", arr, ids=ids)
        observed = pd.Series(observed, index=ids, name="observed_otus")
    except Exception:
        observed = pd.Series((arr > 0).sum(axis=1), index=ids, name="observed_otus")
    metrics = {}
    for m in ["shannon", "simpson", "chao1"]:
        metrics[m] = alpha_diversity(m, arr, ids=ids)
    alpha_df = pd.DataFrame(metrics, index=ids)
    alpha_df.insert(0, "observed_otus", observed)
    return alpha_df

def save_alpha_plots(alpha_df: pd.DataFrame, outdir: Path, level: str):
    sns.set_theme(style="whitegrid")
    for metric in alpha_df.columns:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=alpha_df[metric], orient="h")
        sns.stripplot(data=alpha_df[metric], orient="h", color="black", size=4)
        plt.title(f"Alpha Diversity — {metric} — {level}")
        plt.xlabel(metric)
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_boxplot_{level}.png", dpi=300)
        plt.close()

        plt.figure(figsize=(10, 6))
        sns.violinplot(data=alpha_df[metric], orient="h", inner="box")
        sns.stripplot(data=alpha_df[metric], orient="h", color="black", size=4)
        plt.title(f"Alpha Diversity — {metric} (violin) — {level}")
        plt.xlabel(metric)
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_violin_{level}.png", dpi=300)
        plt.close()


def compute_beta_and_pcoa(df_sxf: pd.DataFrame, metric: str, binarize: bool = False):
    X = df_sxf.values
    if binarize:
        X = (X > 0).astype(int)
    ids = df_sxf.index
    dist = beta_diversity(metric, X, ids=ids)
    ord_res = pcoa(dist)
    coords = ord_res.samples.copy()
    exp = (ord_res.proportion_explained * 100).round(2)
    return dist, coords, exp


def plot_pcoa(coords: pd.DataFrame, exp: pd.Series, out_png: Path, title: str):
    plt.figure(figsize=(7.5, 6.5))
    plt.scatter(coords["PC1"], coords["PC2"], s=80, alpha=0.9)
    plt.xlabel(f"PC1 ({exp.iloc[0]:.2f}%)")
    plt.ylabel(f"PC2 ({exp.iloc[1]:.2f}%)")
    plt.title(f"PCoA — {title}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def make_heatmap(df_sxf: pd.DataFrame, out_png: Path, title: str):
    g = sns.clustermap(df_sxf, cmap="viridis", figsize=(12, 12))
    g.fig.suptitle(title, y=1.02)
    g.savefig(out_png, dpi=300)
    plt.close(g.fig)


def run_permanova(dist_matrix, sample_groups, permutations=999):
    df = pd.DataFrame({"group": sample_groups}, index=dist_matrix.ids)
    result = permanova(dist_matrix, df, column="group", permutations=permutations)
    return result


def infer_level_from_path(p: Path) -> str:
    stem = p.stem.lower()
    keys = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
    for k in keys:
        if k in stem:
            return k
    return stem

# ======================= pipeline =======================

def process_table(path: Path, outdir_base: Path, groups: Optional[Dict[str, str]]):
    print(f"[INFO] Processing: {path}")
    level = infer_level_from_path(path)
    outdir = outdir_base / level
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(path, sep="\t", index_col=0)
    df_sxf = ensure_samples_by_features(df)
    df_sxf = drop_empty_samples(df_sxf)

    # ---- Alpha ----
    alpha_df = compute_alpha(df_sxf)
    alpha_df.to_csv(outdir / f"alpha_metrics_{level}.tsv", sep="\t")
    save_alpha_plots(alpha_df, outdir, level)

    # ---- Heatmap + clustering ----
    make_heatmap(df_sxf, outdir / f"heatmap_{level}.png", f"Heatmap {level}")

    # ---- Beta Bray–Curtis ----
    bray_dist, bray_coords, bray_exp = compute_beta_and_pcoa(df_sxf, "braycurtis")
    pd.DataFrame(bray_dist.data, index=bray_dist.ids, columns=bray_dist.ids)\
        .to_csv(outdir / f"beta_braycurtis_{level}.tsv", sep="\t")
    bray_coords.to_csv(outdir / f"pcoa_braycurtis_coords_{level}.tsv", sep="\t")
    plot_pcoa(bray_coords, bray_exp, outdir / f"pcoa_braycurtis_{level}.png", f"{level} • Bray–Curtis")

    # ---- PERMANOVA (if groups provided) ----
    if groups:
        sample_groups = [groups.get(s, "ungrouped") for s in bray_dist.ids]
        adonis = run_permanova(bray_dist, sample_groups)
        with open(outdir / f"permanova_braycurtis_{level}.txt", "w") as f:
            f.write(str(adonis))
        print("[INFO] PERMANOVA done.")

    # ---- Beta Jaccard ----
    jacc_dist, jacc_coords, jacc_exp = compute_beta_and_pcoa(df_sxf, "jaccard", binarize=True)
    pd.DataFrame(jacc_dist.data, index=jacc_dist.ids, columns=jacc_dist.ids)\
        .to_csv(outdir / f"beta_jaccard_{level}.tsv", sep="\t")
    jacc_coords.to_csv(outdir / f"pcoa_jaccard_coords_{level}.tsv", sep="\t")
    plot_pcoa(jacc_coords, jacc_exp, outdir / f"pcoa_jaccard_{level}.png", f"{level} • Jaccard")

    print(f"[DONE] {level} saved → {outdir}")

# ======================= CLI =======================

def main():
    ap = argparse.ArgumentParser(description="Full diversity analysis pipeline")
    ap.add_argument("--indir", default=".", help="Directory containing abundance tables")
    ap.add_argument("--outdir", default="diversity_outputs", help="Output directory")
    ap.add_argument("--groups", help="Optional TSV mapping samples to groups for PERMANOVA")
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir_base = Path(args.outdir)
    outdir_base.mkdir(parents=True, exist_ok=True)

    # Load groups file if provided
    groups = None
    if args.groups:
        gdf = pd.read_csv(args.groups, sep="\t", header=None, names=["sample", "group"])
        groups = dict(zip(gdf.sample, gdf.group))
        print("[INFO] Groups loaded for PERMANOVA")

    # Auto-detect all *_counts_*.tsv files
    tables = sorted(indir.glob("*counts_*.tsv"))
    if not tables:
        print("[ERROR] No abundance tables found.")
        return

    for table in tables:
        process_table(table, outdir_base, groups)

if __name__ == "__main__":
    sns.set_theme(style="whitegrid")
    main()