"""End-to-end tests for the PDBe-API ("processed") Arpeggio schema.

Added alongside `parse_arpeggio_file_pdbe_template` / `flatten_arpeggio_content`
(see plan-accion-modularizacion-pickit.md, "processed schema" sub-task):
`analyze_files(mode=ARPEGGIO, ...)` must transparently accept PDBe-API JSON
exports (dict-by-ligand-key or dict-by-PDB-id, both without
`interacting_entities`) alongside the original raw Arpeggio schema, without
the caller having to say which schema a given file uses — a directory may
even mix both.

Deliberately self-contained (writes its own tiny JSON fixtures under
`tmp_path`) rather than depending on the `arpeggio_dir`/`template_path`
conftest fixtures used by test_analyze_files.py, since those fixtures'
exact on-disk content isn't available in this review.
"""

import json

import pytest


def _raw_record(residue_number=205, chem_comp_id="GLU", contact="hbond"):
    return {
        "bgn": {
            "auth_asym_id": "A",
            "auth_atom_id": "CA",
            "auth_seq_id": residue_number,
            "label_comp_id": chem_comp_id,
            "label_comp_type": "P",
        },
        "contact": [contact],
        "type": "atom-atom",
        "interacting_entities": "INTER",
        "end": {
            "auth_asym_id": "A",
            "auth_atom_id": "N19",
            "auth_seq_id": 954,
            "label_comp_id": "A3M",
            "label_comp_type": "N",
        },
    }


def _pdbe_record(residue_number=205, chem_comp_id="GLU", detail="hbond"):
    return {
        "ligand_atoms": ["N19"],
        "interaction_type": "atom-atom",
        "interaction_details": [detail],
        "end": {
            "chain_id": "A",
            "author_residue_number": residue_number,
            "chem_comp_id": chem_comp_id,
            "atom_names": ["OE2"],
        },
    }


class TestPdbeSchemaDetection:
    def test_dict_by_ligand_key_file_is_parsed_without_extra_config(self, analyzer, tmp_path):
        pdbe_dir = tmp_path / "pdbe_mini"
        pdbe_dir.mkdir()
        content = {"A_954_A3M": [_pdbe_record()]}
        (pdbe_dir / "complex1.json").write_text(json.dumps(content))

        data = analyzer.analyze_files(directory=str(pdbe_dir), mode=analyzer.ARPEGGIO)
        row_labels = [row[0] for row in data.matrix[1:]]
        assert "GLU 205" in row_labels

    def test_dict_by_pdb_id_file_is_parsed_without_extra_config(self, analyzer, tmp_path):
        pdbe_dir = tmp_path / "pdbe_mini"
        pdbe_dir.mkdir()
        content = {
            "1n1m": [
                {
                    "ligand": {"chain_id": "A", "author_residue_number": 954, "chem_comp_id": "A3M"},
                    "interactions": [_pdbe_record()],
                }
            ]
        }
        (pdbe_dir / "complex1.json").write_text(json.dumps(content))

        data = analyzer.analyze_files(directory=str(pdbe_dir), mode=analyzer.ARPEGGIO)
        row_labels = [row[0] for row in data.matrix[1:]]
        assert "GLU 205" in row_labels

    def test_water_is_excluded_by_default(self, analyzer, tmp_path):
        pdbe_dir = tmp_path / "pdbe_mini"
        pdbe_dir.mkdir()
        content = {
            "A_954_A3M": [
                _pdbe_record(residue_number=205, chem_comp_id="GLU"),
                _pdbe_record(residue_number=1133, chem_comp_id="HOH"),
            ]
        }
        (pdbe_dir / "complex1.json").write_text(json.dumps(content))

        data = analyzer.analyze_files(directory=str(pdbe_dir), mode=analyzer.ARPEGGIO)
        row_labels = [row[0] for row in data.matrix[1:]]
        assert "GLU 205" in row_labels
        assert not any("HOH" in label for label in row_labels)

    def test_mixed_raw_and_pdbe_files_in_the_same_directory(self, analyzer, tmp_path):
        mixed_dir = tmp_path / "mixed"
        mixed_dir.mkdir()
        (mixed_dir / "raw_complex.json").write_text(json.dumps([_raw_record()]))
        (mixed_dir / "pdbe_complex.json").write_text(json.dumps({"A_954_A3M": [_pdbe_record()]}))

        data = analyzer.analyze_files(directory=str(mixed_dir), mode=analyzer.ARPEGGIO)
        # Both files' ligand columns show up, and both land on the same
        # GLU 205 row since they describe the same residue/interaction.
        header = data.matrix[0]
        assert len([c for c in header if "A3M" in c]) == 2
        row_labels = [row[0] for row in data.matrix[1:]]
        assert row_labels == ["GLU 205"]


class TestPdbeSchemaCustomTemplate:
    def test_custom_processed_exclude_rule_is_honored(self, analyzer, tmp_path):
        pdbe_dir = tmp_path / "pdbe_mini"
        pdbe_dir.mkdir()
        content = {
            "A_954_A3M": [
                _pdbe_record(residue_number=205, chem_comp_id="GLU"),
                _pdbe_record(residue_number=999, chem_comp_id="SO4"),
            ]
        }
        (pdbe_dir / "complex1.json").write_text(json.dumps(content))

        template = {
            "raw": [{"contact": ["hbond"], "interacting_entities": "INTER", "type": None}],
            "processed": {"exclude": [{"end": {"chem_comp_id": "SO4"}}]},
        }
        template_path = tmp_path / "custom_template.json"
        template_path.write_text(json.dumps(template))

        data = analyzer.analyze_files(
            directory=str(pdbe_dir), mode=analyzer.ARPEGGIO, template_file=str(template_path)
        )
        row_labels = [row[0] for row in data.matrix[1:]]
        assert "GLU 205" in row_labels
        assert not any("SO4" in label for label in row_labels)

    def test_old_style_bare_list_template_still_works_and_excludes_water_by_default(self, analyzer, tmp_path):
        # Backward compatibility: a template_file that's still the old bare
        # list of raw-schema rules (no "raw"/"processed" dict wrapper) must
        # keep working, with the default water-exclusion rule applying to
        # any processed-schema file it's used against.
        pdbe_dir = tmp_path / "pdbe_mini"
        pdbe_dir.mkdir()
        content = {
            "A_954_A3M": [
                _pdbe_record(residue_number=205, chem_comp_id="GLU"),
                _pdbe_record(residue_number=1133, chem_comp_id="HOH"),
            ]
        }
        (pdbe_dir / "complex1.json").write_text(json.dumps(content))

        old_style_template = [{"contact": ["hbond"], "interacting_entities": "INTER", "type": None}]
        template_path = tmp_path / "old_template.json"
        template_path.write_text(json.dumps(old_style_template))

        data = analyzer.analyze_files(
            directory=str(pdbe_dir), mode=analyzer.ARPEGGIO, template_file=str(template_path)
        )
        row_labels = [row[0] for row in data.matrix[1:]]
        assert "GLU 205" in row_labels
        assert not any("HOH" in label for label in row_labels)
