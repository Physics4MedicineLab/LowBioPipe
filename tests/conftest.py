"""Shared fixtures for the bin/ script tests."""

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / "bin"
TEST_DATA = REPO / "test_data"

sys.path.insert(0, str(BIN))


@pytest.fixture
def taxdump() -> Path:
    """Mini NCBI taxdump shipped with the test dataset."""
    return TEST_DATA / "taxdump"


@pytest.fixture
def test_data() -> Path:
    """Root of the synthetic test dataset."""
    return TEST_DATA
