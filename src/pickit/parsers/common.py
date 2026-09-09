"""Parsing helpers shared by the Arpeggio and IChem parsers.

Extracted from nested closures inside the original ``analyze_files`` (see
``parsers/__init__.py`` for context). No logic changes — ``modify_cell``
and ``validate_string`` had no closure dependencies beyond module-level
constants, so they move over verbatim. ``get_protein_ligand`` closed over
``self.aa`` (the amino-acid-code list, confusingly named the same as the
unrelated ``aa`` dict used elsewhere in ``analyze_files`` for
residue-name -> column-index bookkeeping); it's renamed to the explicit
parameter ``amino_acid_codes`` here to keep that distinction unambiguous
now that it's no longer implicit via closure.
"""

import re

from ..constants import DIFF_DELIM, GROUP_DELIM, SAME_DELIM


def modify_cell(
    text: str,
    interaction: str,
    atoms: str,
    interaction_labels: list,
) -> str:
    """
    Updates the cell content by adding interaction type and involved atoms.

    Args:
        text (str): Current content of the matrix cell.
        interaction (str): Type of interaction occurring.
        atoms (str): Atoms involved in the interaction.
        interaction_labels (list): List of predefined interaction labels.

    Returns:
        str: Updated cell content with formatted interaction information.
    """
    # Create an interaction_map based on the global INTERACTION_LABELS list
    interaction_map = {label: str(index + 1) for index, label in enumerate(interaction_labels)}

    # Assign the interaction code based on the interaction_map, or default to the last value
    interaction_code = interaction_map.get(interaction, str(len(interaction_labels)))

    # If the text is empty, add the interaction with the provided atoms
    if text == "":
        return f"{interaction_code} {GROUP_DELIM}{atoms}{GROUP_DELIM}"

    # Split the cell content and remove empty parts, keeping existing interactions
    content = text.replace(DIFF_DELIM, "").split(GROUP_DELIM)[:-1]
    exists = False
    cell = ""

    # Check if the interaction already exists and add atoms to the corresponding interaction
    for index, segment in enumerate(content):
        if index % 2 == 0 and interaction_code == segment.strip():
            # Delete repate entries
            entries = content[index + 1].split(SAME_DELIM)
            detected = False
            for entry in entries:
                if entry == atoms:
                    detected = True
                    break
            if not detected:
                content[index + 1] += f"{SAME_DELIM}{atoms}"
            exists = True
            break

    # Rebuild the cell content, preserving existing interactions
    cell = DIFF_DELIM.join(f"{content[i]}{GROUP_DELIM}{content[i + 1]}{GROUP_DELIM}" for i in range(0, len(content), 2))

    # If the interaction didn't exist, append it at the end
    if not exists:
        cell += f"{DIFF_DELIM}{interaction_code} {GROUP_DELIM}{atoms}{GROUP_DELIM}"

    return cell


def validate_string(input_string: str) -> bool:
    """
    Validates a residue string to ensure it follows the expected format.

    Format:
    - Starts with a three-letter amino acid abbreviation.
    - Followed by spaces and a numeric sequence.
    - Ends with a dash and additional characters.

    Args:
        input_string (str): The residue string to validate.

    Returns:
        bool: True if the format is valid, False otherwise.
    """
    # Regular expression to validate the required format
    pattern = r"^[A-Z]{3} +\d+-.+$"

    # Match the input string with the regular expression pattern
    return bool(re.match(pattern, input_string))


def is_raw_arpeggio_record(record: dict) -> bool:
    """
    Tells apart the two Arpeggio-family JSON schemas ``analyze_files`` now
    accepts, at the level of a single interaction record.

    - "raw" schema: the original Arpeggio output (``bgn``/``end``, ``contact``,
      ``type``, ``interacting_entities``). This field is only ever present here.
    - "processed" schema: PDBe-API exports (``ligand_atoms``, ``interaction_type``,
      ``interaction_details``, ``end``). Never has ``interacting_entities``.

    Args:
        record (dict): A single interaction record, already unwrapped from
            whatever container ``flatten_arpeggio_content`` found it in.

    Returns:
        bool: True if ``record`` follows the raw Arpeggio schema.
    """
    return "interacting_entities" in record


def parse_ligand_key(key: str) -> dict:
    """
    Parses a PDBe-API-style ligand key (e.g. ``"A_954_A3M"``) into a ligand
    context dict, used when a "processed"-schema JSON is keyed by ligand
    instead of carrying the ligand identity on each interaction record.

    Args:
        key (str): Dict key of the form ``"{chain_id}_{residue_number}_{chem_comp_id}"``.

    Returns:
        dict: ``{"chain_id": str, "author_residue_number": int, "chem_comp_id": str}``,
        matching the field names used in a record's ``"end"`` dict so both can be
        handled uniformly downstream.

    Raises:
        ValueError: If ``key`` doesn't have at least the three expected parts.
    """
    parts = key.split("_")
    if len(parts) < 3:
        raise ValueError(f"Unrecognized ligand key format: {key!r} (expected 'CHAIN_RESNUM_CODE').")
    chain_id, residue_number = parts[0], parts[1]
    chem_comp_id = "_".join(parts[2:])  # ligand codes can themselves contain '_'
    return {
        "chain_id": chain_id,
        "author_residue_number": int(residue_number),
        "chem_comp_id": chem_comp_id,
    }


def flatten_arpeggio_content(content) -> list[tuple[dict, dict | None]]:
    """
    Normalizes any of the JSON shapes accepted by Arpeggio-mode ``analyze_files``
    into a flat list of ``(interaction_record, ligand_context)`` pairs, so both
    parsers (raw-schema and processed-schema) can work off a single flat list
    regardless of how the source file wraps its records.

    Three shapes are recognized:

    - ``list[dict]``: the original, unwrapped Arpeggio output. Every record is
      self-contained (``bgn``/``end`` both identify their own residue/atom), so
      ``ligand_context`` is always ``None`` here.
    - ``dict[str, list[dict]]``: a PDBe-API export keyed by ligand, e.g.
      ``{"A_954_A3M": [...]}``. ``ligand_context`` is parsed from the key via
      ``parse_ligand_key``.
    - ``dict[str, list[dict]]`` where each list entry is itself
      ``{"ligand": {...}, "interactions": [...]}`` (a PDBe-API export keyed by
      PDB id instead of by ligand). ``ligand_context`` is the nested ``"ligand"``
      dict.

    Args:
        content: Parsed JSON content of one interaction file, in any of the
            three shapes above.

    Returns:
        list[tuple[dict, dict | None]]: Flattened ``(record, ligand_context)``
        pairs, in file order. Each ``record`` still carries whichever schema
        (raw/processed) it originally had — use ``is_raw_arpeggio_record`` to
        tell them apart before parsing.
    """
    if isinstance(content, list):
        return [(record, None) for record in content]

    flat: list[tuple[dict, dict | None]] = []
    for key, entries in content.items():
        if not entries:
            continue
        if "interactions" in entries[0]:
            # dict[pdb_id] -> [{"ligand": ..., "interactions": [...]}, ...]
            for entry in entries:
                ligand_context = entry["ligand"]
                flat.extend((record, ligand_context) for record in entry["interactions"])
        else:
            # dict[ligand_key] -> [interaction_record, ...]
            ligand_context = parse_ligand_key(key)
            flat.extend((record, ligand_context) for record in entries)
    return flat


def get_protein_ligand(begin: dict, end: dict, amino_acid_codes: list[str]) -> tuple[dict, dict]:
    """
    Determines which of the two interacting atoms (`begin`, `end`) belongs
    to the protein and which to the ligand, based on ``label_comp_type``
    ('P' for polymer/protein) and, when both are polymer atoms, whether
    their residue code is a standard amino acid.

    Args:
        begin (dict): The first interacting atom's Arpeggio record.
        end (dict): The second interacting atom's Arpeggio record.
        amino_acid_codes (list[str]): Standard three-letter amino acid codes
            (``self.aa`` on the analyzer), used to disambiguate when both
            atoms are polymer atoms.

    Returns:
        tuple[dict, dict]: ``(protein_atom, ligand_atom)``, or ``(None, None)``
        if the pair can't be classified as one protein atom + one ligand atom.
    """
    if begin["label_comp_type"] == "P" and end["label_comp_type"] == "P":
        if begin["label_comp_id"] in amino_acid_codes and end["label_comp_id"] in amino_acid_codes:
            return None, None
        elif begin["label_comp_id"] in amino_acid_codes:
            return begin, end
        else:
            return end, begin
    elif begin["label_comp_type"] == "P":
        return begin, end
    elif end["label_comp_type"] == "P":
        return end, begin
    else:
        return None, None
