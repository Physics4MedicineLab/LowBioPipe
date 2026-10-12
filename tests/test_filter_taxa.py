"""Tests for bin/filter_taxa.py."""

import pytest
from filter_taxa import TaxonomyTree, filter_classified_reads, load_exclude_taxids


@pytest.fixture
def tree(taxdump):
    return TaxonomyTree(taxdump)


@pytest.mark.parametrize(
    "taxid, max_rank, expected",
    [
        (9606, "genus", {9606, 9605}),
        (9606, "family", {9606, 9605, 9604}),
        (46170, "genus", {46170, 1280, 1279}),  # subspecies -> species -> genus
        (2, "genus", {2}),  # already above the cap: only itself
        (9605, "species", {9605}),
    ],
)
def test_ancestors_stop_at_max_rank(tree, taxid, max_rank, expected):
    assert tree.get_ancestors(taxid, max_rank) == expected


@pytest.mark.parametrize("max_rank", ["domain", "superkingdom"])
def test_ancestors_never_reach_root(tree, max_rank):
    ancestors = tree.get_ancestors(9606, max_rank)
    assert 2759 in ancestors  # Eukaryota ('superkingdom' in this taxdump) is the cap
    assert not ancestors & {1, 131567}


def test_descendants_include_subspecies(tree):
    forbidden = tree.expand_forbidden_taxa([1280], include_descendants=True)
    assert forbidden == {1280, 46170}


def test_exact_mode_only_listed_taxa(tree):
    assert tree.expand_forbidden_taxa([1280, 999999]) == {1280}  # unknown TaxID skipped


def test_exclusion_file_comments_and_blank_lines(test_data):
    assert load_exclude_taxids(test_data / "contaminants.txt") == [9606, 316]


def test_filter_classified_reads(tmp_path):
    reads = tmp_path / "S1.classifiedreads.txt"
    reads.write_text(
        "C\tr1\t9606\t150\t9606:10\n"
        "C\tr2\t287\t150\t287:10\n"
        "U\tr3\t0\t150\t0:10\n"
    )
    out = tmp_path / "S1.classifiedreads.filtered.txt"
    total, removed, retained, malformed = filter_classified_reads(reads, out, {9606})
    assert (total, removed, retained, malformed) == (3, 1, 2, 0)
    assert "r1" not in out.read_text()


def test_filter_classified_reads_with_names(tmp_path):
    """Kraken2 --use-names writes 'Name (taxid N)' in the taxid column."""
    reads = tmp_path / "S1.classifiedreads.txt"
    reads.write_text(
        "C\tr1\tHomo sapiens (taxid 9606)\t150\t9606:10\n"
        "C\tr2\tPseudomonas aeruginosa (taxid 287)\t150\t287:10\n"
        "U\tr3\tunclassified (taxid 0)\t150\t0:10\n"
    )
    out = tmp_path / "S1.classifiedreads.filtered.txt"
    total, removed, retained, malformed = filter_classified_reads(reads, out, {9606})
    assert (total, removed, retained, malformed) == (3, 1, 2, 0)
    assert "r1" not in out.read_text()
