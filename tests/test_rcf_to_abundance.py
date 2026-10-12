"""Tests for bin/rcf_to_abundance.py, on output laid out like Recentrifuge's."""

import logging
import subprocess
import sys

import pandas as pd
import pytest
from conftest import BIN
from rcf_to_abundance import ancestor_at_rank, load_taxonomy, read_recentrifuge

RAW_SAMPLES = ["BLANK1", "S1", "patient_01", "patient_02"]

# UNA counts per sample: raw samples plus Recentrifuge's derived pseudo-samples,
# with distinct values so the tests can tell which column was read
UNA = {
    "BLANK1": {1282: 30},
    "S1": {287: 100, 1280: 10, 46170: 5},
    "patient_01": {287: 50, 1280: 4},
    "patient_02": {562: 70, 287: 7},
    "S1_CTRL_species": {287: 1},
    "S1_EXCLUSIVE_SUMMARY": {1280: 2},
    "S1_CTRL_SUMMARY": {287: 90, 1280: 15},
    "patient_01_CTRL_SUMMARY": {287: 45},
    "patient_02_CTRL_SUMMARY": {562: 60},
    "SHARED_SUMMARY": {287: 999},
}
TAXA = {
    287: ("species", "Pseudomonas aeruginosa"),
    562: ("species", "Escherichia coli"),
    1280: ("species", "Staphylococcus aureus"),
    1282: ("species", "Staphylococcus epidermidis"),
    46170: ("subspecies", "Staphylococcus aureus subsp. aureus"),
}


@pytest.fixture(params=["tsv", "xlsx"])
def rcf_output(request, tmp_path):
    """Write Recentrifuge output the way its save_extra_output() does (-e TSV or FULL)."""
    samples = list(UNA)
    columns = pd.MultiIndex.from_product(
        [samples, ["cnt", "una", "sco"]], names=["Samples", "Stats"]
    ).append(pd.MultiIndex.from_product([["Details"], ["Rank", "Name"]]))
    rows = {
        taxid: [v for s in samples for v in (UNA[s].get(taxid, 0),) * 2 + (30.0,)]
        + [rank, name]
        for taxid, (rank, name) in TAXA.items()
    }
    full = pd.DataFrame.from_dict(rows, orient="index", columns=columns)
    full.index.names = ["Id"]
    stats = pd.DataFrame({s: [1000] for s in RAW_SAMPLES}, index=["seq_read"])
    if request.param == "tsv":
        full.to_csv(tmp_path / "test.rcf.data.tsv", sep="\t")
        stats.to_csv(tmp_path / "test.rcf.stat.tsv", sep="\t")
        return tmp_path / "test.rcf.data.tsv"
    path = tmp_path / "test.rcf.xlsx"
    with pd.ExcelWriter(path) as writer:
        stats.to_excel(writer, sheet_name="_sample_stats")
        full.to_excel(writer, sheet_name="FULL")
    return path


def run_script(args, cwd):
    return subprocess.run(
        [sys.executable, str(BIN / "rcf_to_abundance.py"), *map(str, args)],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def test_raw_source_keeps_every_raw_sample(rcf_output):
    df = read_recentrifuge(str(rcf_output), "raw")
    df.index = df.pop("Taxon").astype(int)
    assert list(df.columns) == RAW_SAMPLES  # patient_01/patient_02 not collapsed
    assert df.loc[287, "S1"] == 100
    assert df.loc[287, "patient_02"] == 7


def test_ctrl_source_uses_ctrl_summary_and_skips_controls(rcf_output):
    df = read_recentrifuge(str(rcf_output), "ctrl", controls=["BLANK1"])
    df.index = df.pop("Taxon").astype(int)
    assert list(df.columns) == ["S1", "patient_01", "patient_02"]
    assert df.loc[287, "S1"] == 90
    assert df.loc[562, "patient_02"] == 60


def test_ctrl_source_warns_about_dropped_samples(rcf_output, caplog):
    # Pretend patient_02 is a sample whose reads were all removed by control filtering
    with caplog.at_level(logging.INFO):
        read_recentrifuge(str(rcf_output), "ctrl", controls=["BLANK1", "patient_02"])
        assert not any(r.levelno == logging.WARNING for r in caplog.records)
        read_recentrifuge(str(rcf_output), "ctrl", controls=[])
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1 and "BLANK1" in warnings[0]


def test_ancestor_at_rank(taxdump):
    parent, rank, _ = load_taxonomy(taxdump)
    assert ancestor_at_rank(46170, "species", parent, rank) == 1280
    assert ancestor_at_rank(287, "phylum", parent, rank) == 1224
    assert ancestor_at_rank(2, "genus", parent, rank) is None


def test_species_aggregation_folds_subspecies(rcf_output, taxdump, tmp_path):
    result = run_script(
        [rcf_output, "--taxdump", taxdump, "--rank", "species", "--aggregate",
         "--min-samples", "1", "--out-prefix", "t"],
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    counts = pd.read_csv(tmp_path / "t_counts_species.tsv", sep="\t", index_col=0)
    assert counts.loc["Staphylococcus aureus [1280]", "S1"] == 15
    assert "SHARED" not in counts.columns


def test_phylum_aggregation(rcf_output, taxdump, tmp_path):
    result = run_script(
        [rcf_output, "--taxdump", taxdump, "--source", "ctrl", "--rank", "phylum",
         "--aggregate", "--min-samples", "1", "--relative", "--out-prefix", "t"],
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    counts = pd.read_csv(tmp_path / "t_counts_phylum.tsv", sep="\t", index_col=0)
    assert counts.loc["Pseudomonadota [1224]", "S1"] == 90
    assert counts.loc["Bacillota [1239]", "S1"] == 15
    rel = pd.read_csv(tmp_path / "t_relative_phylum.tsv", sep="\t", index_col=0)
    assert rel.sum().round(6).eq(1).all()


def test_no_taxa_left_fails_cleanly(rcf_output, taxdump, tmp_path):
    result = run_script(
        [rcf_output, "--taxdump", taxdump, "--rank", "species", "--aggregate",
         "--min-samples", "10", "--out-prefix", "t"],
        tmp_path,
    )
    assert result.returncode == 1
    assert "No taxa at rank 'species'" in result.stderr
