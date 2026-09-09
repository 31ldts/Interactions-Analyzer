"""
Shared fixtures for the Phase 0 baseline test suite.

These tests run against the current monolithic `analyze_interactions.py`
exactly as it is today. They must keep passing, unmodified, throughout the
Phase 2 modularization (see plan-accion-modularizacion-pickit.md) except
where a test itself is found to be wrong, which must be documented.
"""
import copy
import os
import sys

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(TESTS_DIR)
DATA_DIR = os.path.join(TESTS_DIR, "data")

# src-layout: the package is installed (e.g. `pip install -e .`) rather than
# reached via a path hack. Keep a fallback for running tests straight out of
# a checkout without an editable install.
sys.path.insert(0, os.path.join(PROJECT_ROOT, "src"))

from pickit.analyze_interactions import (  # noqa: E402
    AnalyzeInteractions,
    InteractionData,
    TypeMismatchException,
    FileOrDirectoryException,
    InvalidColorException,
    InvalidFileExtensionException,
    InvalidAxisException,
    InvalidFilenameException,
    HeatmapActivityException,
    InvalidModeException,
    MissedActivityException,
)

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(TESTS_DIR)

@pytest.fixture
def analyzer(tmp_path):
    """A fresh AnalyzeInteractions, with input/output dirs pointed at test data / a tmp dir."""
    a = AnalyzeInteractions()
    a.input_directory = DATA_DIR
    a.saving_directory = str(tmp_path)
    return a


@pytest.fixture
def small_matrix():
    """
    A small, deterministic, hand-built interaction matrix: 3 residues x 3 ligands.

    Cell format matches what analyze_files produces:
        "<label_index> |atom-atom(subunit)|" and, for multiple interaction
        types in one cell, joined with '; '.

    Layout (residues in rows, ligands in columns):
                 LIGA (7.5)              LIGB (5.0)          LIGC (0)
        HIS 41   "13 |CA-O1(A)|"         ""                  ""
        ASN 142  "8 |CB-C2(A)|; 13 |CB-N1(A)|"  "13 |CB-N1(A)|"  ""
        GLU 166  ""                      ""                  ""

    Interaction label 8 = 'hydrophobic', 13 = 'hbond' under the default
    arpeggio label ordering used elsewhere in this project.
    """
    matrix = [
        ["", "LIGA (7.5)", "LIGB (5.0)", "LIGC (0)"],
        ["HIS 41", "13 |CA-O1(A)|", "", ""],
        ["ASN 142", "8 |CB-C2(A)|; 13 |CB-N1(A)|", "13 |CB-N1(A)|", ""],
        ["GLU 166", "", "", ""],
    ]
    return matrix


@pytest.fixture
def small_interaction_data(small_matrix):
    """InteractionData wrapping `small_matrix`, using the default arpeggio interaction labels."""
    interactions = (
        ["AMIDEAMIDE", "AMIDERING", "plane-plane", "CARBONPI", "CATIONPI", "METSULPHURPI",
         "covalent", "hydrophobic", "carbonyl", "aromatic", "metal", "ionic", "hbond",
         "DONORPI", "weak_hbond", "polar", "weak_polar", "xbond", "HALOGENPI"]
    )
    colors = ["#e6194B"] * len(interactions)
    return InteractionData(
        colors=colors,
        interactions=interactions,
        ligand=True,
        matrix=copy.deepcopy(small_matrix),
        mode="arpeggio",
        protein=True,
        subunit=False,
        subunits_set={"A"},
    )


@pytest.fixture
def arpeggio_dir():
    """Directory (relative to DATA_DIR) with 2 tiny synthetic Arpeggio JSON files."""
    return "arpeggio_mini"


@pytest.fixture
def activity_csv():
    """Path (relative to DATA_DIR) to a tiny activity CSV matching arpeggio_mini ligands."""
    return "mini_activity.csv"


@pytest.fixture
def template_path():
    return os.path.join(PROJECT_ROOT, "src", "pickit", "template_arpeggio.json")
