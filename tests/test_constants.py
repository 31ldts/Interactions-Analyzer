"""Tests for pickit.constants.

These are new tests introduced during Fase 2, step 1 (extraction of
constants.py). They didn't exist as a dedicated Fase 0 baseline test file
because the constants had no dedicated public API before — they were just
module globals inside analyze_interactions.py, already exercised indirectly
by the baseline suite. These smoke tests lock down the extracted module's
shape so a future edit can't silently drop or rename a constant that other
modules depend on.
"""

from pickit import constants


def test_interaction_labels_and_colors_are_same_length_family():
    # Not required to match 1:1 (COLORS is reused across template modes),
    # but both must be non-empty lists of the expected element type.
    assert constants.INTERACTION_LABELS
    assert all(isinstance(label, str) for label in constants.INTERACTION_LABELS)
    assert constants.COLORS
    assert all(isinstance(color, str) for color in constants.COLORS)


def test_program_modes_contains_arpeggio_and_ichem():
    assert set(constants.PROGRAM_MODES) == {"arpeggio", "ichem"}


def test_delimiters_are_distinct():
    delims = {constants.SAME_DELIM, constants.DIFF_DELIM, constants.GROUP_DELIM}
    assert len(delims) == 3


def test_empty_cell_constants():
    assert constants.EMPTY_CELL == ""
    assert constants.EMPTY_DASH_CELL == "-"


def test_is_not_empty_or_dash_lambda():
    assert constants.is_not_empty_or_dash("13 |CA-O1(A)|") is True
    assert constants.is_not_empty_or_dash(constants.EMPTY_CELL) is False
    assert constants.is_not_empty_or_dash(constants.EMPTY_DASH_CELL) is False


def test_amino_acid_codes_are_three_letter_uppercase():
    assert all(len(code) == 3 and code.isupper() for code in constants.AMINO_ACID_CODES)


def test_default_template_file_name():
    # Documents the defaults used by analyze_files(template_file=...); the
    # caller can always override them with an explicit path.
    #
    # Renamed from a single "template.json" once IChem got its own default
    # template too: DEFAULT_ARPEGGIO_TEMPLATE_FILE is the new canonical
    # name, DEFAULT_TEMPLATE_FILE is kept only as a backward-compatible
    # alias pointing at the same value.
    assert constants.DEFAULT_ARPEGGIO_TEMPLATE_FILE == "template_arpeggio.json"
    assert constants.DEFAULT_TEMPLATE_FILE == constants.DEFAULT_ARPEGGIO_TEMPLATE_FILE
    assert constants.DEFAULT_ICHEM_TEMPLATE_FILE == "template_ichem.json"


def test_reexported_from_analyze_interactions():
    # Backward compatibility: anything importing these names the old way
    # (from pickit.analyze_interactions import X) must keep working.
    from pickit import analyze_interactions as ai

    assert ai.INTERACTION_LABELS is constants.INTERACTION_LABELS
    assert ai.is_not_empty_or_dash is constants.is_not_empty_or_dash
