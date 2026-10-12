"""Tests for bin/diversity_analysis.py."""

import os
import re
import subprocess
import sys
from pathlib import Path

import diversity_analysis as da
import numpy as np
import pandas as pd
import pytest
import skbio
from conftest import BIN

SKBIO_VERSION = tuple(int(x) for x in skbio.__version__.split(".")[:3])


@pytest.fixture
def counts_dir(tmp_path):
    """Abundance table (taxa x samples) with numeric-looking sample IDs, and a groups file."""
    rng = np.random.default_rng(3)
    counts = rng.poisson(lam=rng.uniform(1, 40, size=30), size=(12, 30))
    samples = [f"{i:02d}" for i in range(12)]
    table = pd.DataFrame(counts.T, index=[f"taxon{i}" for i in range(30)], columns=samples)
    indir = tmp_path / "in"
    indir.mkdir()
    table.to_csv(indir / "abundance_x_counts_species.tsv", sep="\t")
    # 11 grouped samples; '11' is missing from the groups file
    (tmp_path / "groups.tsv").write_text(
        "".join(f"{s}\t{'AB'[i % 2]}\n" for i, s in enumerate(samples[:11]))
    )
    return tmp_path


def run_script(counts_dir, outdir):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return subprocess.run(
        [sys.executable, str(BIN / "diversity_analysis.py"), "--indir", "in",
         "--outdir", outdir, "--groups", "groups.tsv", "--heatmap-top-n", "10"],
        cwd=counts_dir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.skipif(
    SKBIO_VERSION < (0, 7, 4),
    reason="scikit-bio < 0.7.4 permdisp divides by zero on tied samples (pinned version is 0.7.4)",
)
def test_group_tests_reproducible_and_exclude_ungrouped(counts_dir):
    for outdir in ("out1", "out2"):
        result = run_script(counts_dir, outdir)
        assert result.returncode == 0, result.stderr

    species = counts_dir / "out1" / "species"
    for test in ("permanova", "permdisp"):
        for metric in ("braycurtis", "jaccard", "aitchison"):
            text = (species / f"{test}_{metric}_species.txt").read_text()
            assert re.search(r"sample size\s+11\n", text)  # '01'-style IDs matched, '11' excluded
            assert text == (counts_dir / "out2" / "species" / f"{test}_{metric}_species.txt").read_text()
    assert (species / "heatmap_species.png").exists()


def test_single_group_skips_tests(tmp_path):
    dist = da.DistanceMatrix(np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float), ids=["a", "b", "c"])
    coords = pd.DataFrame({"PC1": [0.0, 1.0, 2.0], "PC2": [0.0, 1.0, 0.0]}, index=["a", "b", "c"])
    exp = pd.Series([60.0, 40.0])
    da.save_beta_results(dist, coords, exp, tmp_path, "species", "braycurtis", {"a": "A", "b": "A", "c": "A"})
    assert not list(tmp_path.glob("perm*"))
    assert (tmp_path / "pcoa_braycurtis_species.png").exists()


def test_heatmap_keeps_top_n_by_mean_relative_abundance(monkeypatch):
    captured = {}

    class FakeGrid:
        fig = type("Fig", (), {"suptitle": lambda *a, **k: None})()

        def savefig(self, *args, **kwargs):
            pass

    def fake_clustermap(data, **kwargs):
        captured["columns"] = list(data.columns)
        return FakeGrid()

    monkeypatch.setattr(da.sns, "clustermap", fake_clustermap)
    # 'deep' dominates sample a (huge depth), 'rare' is a large share of the shallow samples
    df = pd.DataFrame(
        {"deep": [1000, 1, 1], "rare": [10, 5, 5], "tiny": [1, 0, 1]}, index=["a", "b", "c"]
    )
    da.make_heatmap(df, Path("unused.png"), "t", top_n=2)
    assert captured["columns"] == ["rare", "deep"]


@pytest.mark.parametrize(
    "name, level",
    [
        ("abundance_lowbiopipe_recentrifuge_counts_species.tsv", "species"),
        ("abundance_x_counts_superkingdom.tsv", "superkingdom"),
        ("abundance_x_counts_any.tsv", "any"),
    ],
)
def test_infer_level_from_path(name, level):
    assert da.infer_level_from_path(Path(name)) == level
