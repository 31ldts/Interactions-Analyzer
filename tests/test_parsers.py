"""Tests for pickit.parsers (common helpers, arpeggio, ichem).

New tests introduced during Fase 3, when the parsers were split out of
io_mixin.py. Arpeggio behavior is already exercised end-to-end through
test_analyze_files.py (Fase 0 baseline, frozen); these tests target the
parser functions directly and in isolation, including an IChem case using
hand-built line content — the Fase 0 notes explicitly flagged that no
IChem sample data was available among the provided files, so
`parse_ichem_file` had 0% dedicated coverage until now.

Extended for the PDBe-API Arpeggio schema and the IChem template feature
(see plan-accion-modularizacion-pickit.md, "processed schema" sub-task):
`flatten_arpeggio_content` / `is_raw_arpeggio_record` / `parse_ligand_key`
(common.py), `parse_arpeggio_file_pdbe_template` (arpeggio.py), and
`parse_ichem_file`'s new `interaction_list` parameter.
"""

from pickit.parsers.arpeggio import parse_arpeggio_file_pdbe_template
from pickit.parsers.common import (
    flatten_arpeggio_content,
    get_protein_ligand,
    is_raw_arpeggio_record,
    modify_cell,
    parse_ligand_key,
    validate_string,
)
from pickit.parsers.ichem import parse_ichem_file


class TestModifyCell:
    def test_empty_cell_adds_first_interaction(self):
        result = modify_cell(text="", interaction="hbond", atoms="CA-O1(A)", interaction_labels=["hbond", "ionic"])
        assert result == "1 |CA-O1(A)|"

    def test_appends_new_interaction_type_to_existing_cell(self):
        result = modify_cell(
            text="1 |CA-O1(A)|", interaction="ionic", atoms="CB-N1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "1 |CA-O1(A)|; 2 |CB-N1(A)|"

    def test_adds_atoms_to_existing_interaction_type(self):
        result = modify_cell(
            text="1 |CA-O1(A)|", interaction="hbond", atoms="CB-N1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "1 |CA-O1(A), CB-N1(A)|"

    def test_does_not_duplicate_identical_atom_pair(self):
        result = modify_cell(
            text="1 |CA-O1(A)|", interaction="hbond", atoms="CA-O1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "1 |CA-O1(A)|"

    def test_unknown_interaction_falls_back_to_last_label(self):
        result = modify_cell(
            text="", interaction="not_a_label", atoms="CA-O1(A)", interaction_labels=["hbond", "ionic"]
        )
        assert result == "2 |CA-O1(A)|"


class TestValidateString:
    def test_accepts_well_formed_residue_string(self):
        assert validate_string("HIS 41-A") is True

    def test_rejects_missing_dash(self):
        assert validate_string("HIS 41") is False

    def test_rejects_lowercase_residue_code(self):
        assert validate_string("his 41-A") is False


class TestGetProteinLigand:
    AA = ["HIS", "ASN", "GLU"]

    def test_returns_protein_then_ligand_when_only_begin_is_protein(self):
        begin = {"label_comp_type": "P", "label_comp_id": "HIS"}
        end = {"label_comp_type": "N", "label_comp_id": "LIG"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is begin and lig is end

    def test_returns_protein_then_ligand_when_only_end_is_protein(self):
        begin = {"label_comp_type": "N", "label_comp_id": "LIG"}
        end = {"label_comp_type": "P", "label_comp_id": "HIS"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is end and lig is begin

    def test_both_polymer_disambiguates_by_amino_acid_code(self):
        begin = {"label_comp_type": "P", "label_comp_id": "XYZ"}  # not a standard amino acid
        end = {"label_comp_type": "P", "label_comp_id": "HIS"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is end and lig is begin

    def test_both_polymer_and_both_standard_amino_acids_returns_none(self):
        begin = {"label_comp_type": "P", "label_comp_id": "HIS"}
        end = {"label_comp_type": "P", "label_comp_id": "ASN"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is None and lig is None

    def test_neither_is_protein_returns_none(self):
        begin = {"label_comp_type": "N", "label_comp_id": "LIG"}
        end = {"label_comp_type": "N", "label_comp_id": "LG2"}
        prot, lig = get_protein_ligand(begin=begin, end=end, amino_acid_codes=self.AA)
        assert prot is None and lig is None


class TestParseIchemFile:
    def test_single_interaction_line_populates_matrix(self):
        # IChem line format: 10 '|'-separated fields; [0]=interaction,
        # [1]=protein atom, [3]=residue ("RES NUM-SUBUNIT"), [4]=ligand atom.
        # "Hydrophobic" is INTERACTION_LABELS[0] (IChem's own label set, not Arpeggio's).
        line = "Hydrophobic|CA|x|HIS 41-A|O1|x|x|x|x|x"
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=[line],
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert aa == {"HIS 41": 0}
        assert cont == 1
        assert subunits_set == {"A"}
        assert matrix[0][0] == "1 |CA-O1(A)|"

    def test_malformed_line_is_ignored(self):
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=["not|enough|fields"],
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert matrix == []
        assert aa == {}

    def test_subunit_true_keeps_subunit_in_residue_label(self):
        line = "Hydrophobic|CA|x|HIS 41-A|O1|x|x|x|x|x"
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=[line],
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=True,
        )
        assert aa == {"HIS 41-A": 0}
        # subunits_set is only populated in the not-subunit branch
        assert subunits_set == set()


class TestParseIchemFileTemplate:
    """`interaction_list=None` (the default) must reproduce the exact
    pre-template behavior above; passing a list restricts which lines are
    kept and which numeric codes get assigned, mirroring the role
    `interaction_list` plays for the Arpeggio parsers."""

    LINES = [
        "Hydrophobic|CA|x|HIS 41-A|O1|x|x|x|x|x",
        "Ionic_PROT|CB|x|ASN 142-A|N1|x|x|x|x|x",
    ]

    def test_no_template_keeps_every_line_with_default_codes(self):
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=self.LINES,
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
        )
        # "Ionic_PROT" is INTERACTION_LABELS[5] -> code 6, same as before
        # this parameter existed.
        assert aa == {"HIS 41": 0, "ASN 142": 1}
        assert matrix == [["1 |CA-O1(A)|"], ["6 |CB-N1(A)|"]]

    def test_template_drops_lines_outside_the_list(self):
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=self.LINES,
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
            interaction_list=["Hydrophobic"],
        )
        # The Ionic_PROT line is dropped entirely -> no ASN 142 row at all.
        assert aa == {"HIS 41": 0}
        assert matrix == [["1 |CA-O1(A)|"]]

    def test_template_reorders_numeric_codes(self):
        # With a restricted/reordered template, codes follow the template's
        # own order rather than the full INTERACTION_LABELS order.
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=self.LINES,
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
            interaction_list=["Ionic_PROT", "Hydrophobic"],
        )
        assert matrix[aa["HIS 41"]][0] == "2 |CA-O1(A)|"
        assert matrix[aa["ASN 142"]][0] == "1 |CB-N1(A)|"

    def test_empty_template_drops_every_line(self):
        matrix, aa, cont, subunits_set = parse_ichem_file(
            content=self.LINES,
            index=0,
            files=["complex1.txt"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            protein=True,
            ligand=True,
            subunit=False,
            interaction_list=[],
        )
        assert matrix == []
        assert aa == {}


class TestIsRawArpeggioRecord:
    def test_raw_schema_record_is_detected(self):
        assert is_raw_arpeggio_record({"interacting_entities": "INTER", "bgn": {}, "end": {}}) is True

    def test_processed_schema_record_is_not_raw(self):
        record = {"ligand_atoms": ["N19"], "interaction_type": "atom-atom", "end": {}}
        assert is_raw_arpeggio_record(record) is False


class TestParseLigandKey:
    def test_splits_chain_resnum_and_code(self):
        assert parse_ligand_key("A_954_A3M") == {
            "chain_id": "A",
            "author_residue_number": 954,
            "chem_comp_id": "A3M",
        }

    def test_ligand_code_can_itself_contain_underscores(self):
        assert parse_ligand_key("B_12_A_3M") == {
            "chain_id": "B",
            "author_residue_number": 12,
            "chem_comp_id": "A_3M",
        }

    def test_malformed_key_raises_value_error(self):
        import pytest

        with pytest.raises(ValueError):
            parse_ligand_key("not_a_key")


class TestFlattenArpeggioContent:
    def test_flat_list_shape_is_passed_through_with_no_ligand_context(self):
        content = [{"interacting_entities": "INTER", "bgn": {}, "end": {}}]
        flat = flatten_arpeggio_content(content)
        assert flat == [(content[0], None)]

    def test_dict_by_ligand_key_shape_resolves_context_from_the_key(self):
        record = {"ligand_atoms": ["N19"], "interaction_type": "atom-atom", "end": {"chem_comp_id": "GLU"}}
        content = {"A_954_A3M": [record]}
        flat = flatten_arpeggio_content(content)
        assert flat == [(record, {"chain_id": "A", "author_residue_number": 954, "chem_comp_id": "A3M"})]

    def test_dict_by_pdb_id_shape_resolves_context_from_the_nested_ligand(self):
        ligand_ctx = {"chain_id": "A", "author_residue_number": 954, "chem_comp_id": "A3M"}
        record = {"ligand_atoms": ["N19"], "interaction_type": "atom-atom", "end": {"chem_comp_id": "GLU"}}
        content = {"1n1m": [{"ligand": ligand_ctx, "interactions": [record]}]}
        flat = flatten_arpeggio_content(content)
        assert flat == [(record, ligand_ctx)]

    def test_empty_ligand_entry_list_is_skipped_without_error(self):
        content = {"A_954_A3M": []}
        assert flatten_arpeggio_content(content) == []


class TestParseArpeggioFilePdbeTemplate:
    INTERACTION_LIST = ["hbond", "hydrophobic", "vdw_clash"]

    def _record(self, chem_comp_id="GLU", details="hbond", ligand_atoms=None, atom_names=None):
        return {
            "ligand_atoms": ligand_atoms or ["N19"],
            "interaction_type": "atom-atom",
            "interaction_details": [details] if isinstance(details, str) else details,
            "end": {
                "chain_id": "A",
                "author_residue_number": 205,
                "chem_comp_id": chem_comp_id,
                "atom_names": atom_names or ["OE2"],
            },
        }

    def _ligand_context(self, chem_comp_id="A3M"):
        return {"chain_id": "A", "author_residue_number": 954, "chem_comp_id": chem_comp_id}

    def test_single_record_populates_matrix_and_ligand_code(self):
        records = [(self._record(), self._ligand_context())]
        matrix, ligand_code, aa, cont, subunits_set = parse_arpeggio_file_pdbe_template(
            records=records,
            index=0,
            files=["f.json"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            exclude_rules=[{"end": {"chem_comp_id": "HOH"}}],
            interaction_list=self.INTERACTION_LIST,
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert ligand_code == "A3M"
        assert aa == {"GLU 205": 0}
        assert matrix == [["1 |OE2-N19(A)|"]]
        assert subunits_set == {"A"}

    def test_excluded_record_is_dropped(self):
        records = [(self._record(chem_comp_id="HOH"), self._ligand_context())]
        matrix, ligand_code, aa, cont, subunits_set = parse_arpeggio_file_pdbe_template(
            records=records,
            index=0,
            files=["f.json"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            exclude_rules=[{"end": {"chem_comp_id": "HOH"}}],
            interaction_list=self.INTERACTION_LIST,
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert matrix == []
        assert aa == {}

    def test_raw_schema_records_are_skipped_not_double_counted(self):
        raw_record = {"interacting_entities": "INTER", "bgn": {}, "end": {}}
        records = [(raw_record, None), (self._record(), self._ligand_context())]
        matrix, ligand_code, aa, cont, subunits_set = parse_arpeggio_file_pdbe_template(
            records=records,
            index=0,
            files=["f.json"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            exclude_rules=[{"end": {"chem_comp_id": "HOH"}}],
            interaction_list=self.INTERACTION_LIST,
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert aa == {"GLU 205": 0}

    def test_interaction_not_in_vocabulary_is_dropped_but_residue_still_created(self):
        # Mirrors parse_arpeggio_file_template's own behavior: an unlisted
        # interaction string doesn't populate a cell, but the residue row is
        # still created since the record itself wasn't excluded.
        records = [(self._record(details="not_in_vocabulary"), self._ligand_context())]
        matrix, ligand_code, aa, cont, subunits_set = parse_arpeggio_file_pdbe_template(
            records=records,
            index=0,
            files=["f.json"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            exclude_rules=[{"end": {"chem_comp_id": "HOH"}}],
            interaction_list=self.INTERACTION_LIST,
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert aa == {"GLU 205": 0}
        assert matrix == [[""]]

    def test_subunit_true_keeps_chain_in_residue_label(self):
        records = [(self._record(), self._ligand_context())]
        matrix, ligand_code, aa, cont, subunits_set = parse_arpeggio_file_pdbe_template(
            records=records,
            index=0,
            files=["f.json"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            exclude_rules=[{"end": {"chem_comp_id": "HOH"}}],
            interaction_list=self.INTERACTION_LIST,
            protein=True,
            ligand=True,
            subunit=True,
        )
        assert aa == {"GLU 205-A": 0}
        assert matrix == [["1 |OE2-N19|"]]

    def test_multiple_atoms_are_comma_joined(self):
        records = [(self._record(atom_names=["OE1", "OE2"], ligand_atoms=["N19", "C20"]), self._ligand_context())]
        matrix, ligand_code, aa, cont, subunits_set = parse_arpeggio_file_pdbe_template(
            records=records,
            index=0,
            files=["f.json"],
            subunits_set=set(),
            cont=0,
            matrix=[],
            aa={},
            exclude_rules=[{"end": {"chem_comp_id": "HOH"}}],
            interaction_list=self.INTERACTION_LIST,
            protein=True,
            ligand=True,
            subunit=False,
        )
        assert matrix == [["1 |OE1,OE2-N19,C20(A)|"]]
