"""Arpeggio interaction file parsers (plain and template-restricted).

Extracted from the ``arpeggio_analysis`` / ``arpeggio_analysis_template``
closures nested inside the original ``analyze_files`` (see
``parsers/__init__.py`` for context). No logic changes: ``protein``,
``ligand``, ``subunit`` were captured from the enclosing method call via
closure and are now explicit parameters instead, and ``get_protein_ligand``
now takes ``amino_acid_codes`` explicitly rather than reading ``self.aa``
through closure.
"""

from ..constants import ARPEGGIO_CONT, ARPEGGIO_INT_ENT, ARPEGGIO_TYPE
from .common import get_protein_ligand, is_raw_arpeggio_record, modify_cell


def parse_arpeggio_file(
    content: list[dict],
    index: int,
    files: list[str],
    subunits_set: set,
    cont: int,
    matrix: list[list[str]],
    aa: dict,
    protein: bool,
    ligand: bool,
    subunit: bool,
    amino_acid_codes: list[str],
) -> tuple[list[list[str]], str, dict, int, set]:
    """
    Processes an Arpeggio interaction file's content (no template
    restriction), extracting relevant data and updating the (in-progress)
    interaction matrix.

    Args:
        content (list[dict]): Parsed JSON records from an Arpeggio output file.
        index (int): Column index for this file within the overall matrix.
        files (list[str]): All filenames being processed (used only to size new rows).
        subunits_set (set): Accumulator of subunit identifiers seen so far.
        cont (int): Next free row index to assign to a newly-seen residue.
        matrix (list[list[str]]): The in-progress interaction matrix (rows = residues).
        aa (dict): Residue name -> row index, built up across files.
        protein (bool): Whether to include the protein atom in the cell's atom string.
        ligand (bool): Whether to include the ligand atom in the cell's atom string.
        subunit (bool): Whether to keep subunits distinct in the residue label.
        amino_acid_codes (list[str]): Standard three-letter amino acid codes,
            passed through to ``get_protein_ligand``.

    Returns:
        tuple[list[list[str]], str, dict, int, set]:
            Updated ``(matrix, ligand_code, aa, cont, subunits_set)``.
    """
    interaction_list = ARPEGGIO_CONT[:2] + ARPEGGIO_TYPE + ARPEGGIO_CONT[2:]
    # interaction_list.sort(key=lambda s: (s.lower(), s.islower()))

    # Filter to obtain entries with interacting_entities == INTER
    inter_set = [elem for elem in content if elem["interacting_entities"] in ARPEGGIO_INT_ENT]
    # Filter to obtain entries with the desired contact or type
    inter_set = [
        elem
        for elem in inter_set
        for contact in elem["contact"]
        if contact in ARPEGGIO_CONT or elem["type"] in ARPEGGIO_TYPE
    ]
    for inter in inter_set:
        if inter["type"] in ARPEGGIO_TYPE:
            contact = [inter["type"]]
        else:
            contact = [conta for conta in inter["contact"] if conta in ARPEGGIO_CONT]
        prot, lig = get_protein_ligand(begin=inter["bgn"], end=inter["end"], amino_acid_codes=amino_acid_codes)
        if prot is not None:
            residue = prot["label_comp_id"] + " " + str(prot["auth_seq_id"])
            prot_atom = prot["auth_atom_id"]
            prot_subunit = prot["auth_asym_id"]
            ligand_code = lig["label_comp_id"]
            lig_atom = lig["auth_atom_id"]

            subunits_set.add(prot_subunit)
            atoms = (
                f"{prot_atom}-{lig_atom}"
                if protein and ligand
                else prot_atom
                if protein
                else lig_atom
                if ligand
                else ""
            )

            if subunit:
                residue += "-" + prot_subunit
            else:
                atoms += f"({prot_subunit})"

            if residue not in aa:
                aa[residue] = cont
                cont += 1
            column = aa[residue]

            # Ensure matrix size and modify cell
            if len(matrix) <= column:
                matrix.append([""] * len(files))

            for interaction in contact:
                matrix[column][index] = modify_cell(
                    text=matrix[column][index],
                    interaction=interaction,
                    atoms=atoms,
                    interaction_labels=interaction_list,
                )

    return matrix, ligand_code, aa, cont, subunits_set


def _matches_template(entry, template) -> bool:
    """Recursively compares an Arpeggio entry against one template rule."""
    if isinstance(template, dict):
        if not isinstance(entry, dict):
            return False
        for key, tmpl_val in template.items():
            if key not in entry:
                return False
            if not _matches_template(entry[key], tmpl_val):
                return False
        return True
    elif isinstance(template, list):
        if template is None:
            return True
        # Aquí asumimos que entry debe ser un valor único que esté en la lista
        return any(value in template for value in entry)
    elif template is None:
        return True
    else:
        return entry == template


def parse_arpeggio_file_template(
    content: list[dict],
    index: int,
    files: list[str],
    subunits_set: set,
    cont: int,
    matrix: list[list[str]],
    aa: dict,
    template: list[dict],
    interaction_list: list[str],
    protein: bool,
    ligand: bool,
    subunit: bool,
    amino_acid_codes: list[str],
) -> tuple[list[list[str]], str, dict, int, set]:
    """
    Processes an Arpeggio interaction file's content, restricted to entries
    matching a user-supplied template (see ``template.json`` /
    ``template_file`` in ``analyze_files``).

    Args:
        content (list[dict]): Parsed JSON records from an Arpeggio output file.
        index (int): Column index for this file within the overall matrix.
        files (list[str]): All filenames being processed (used only to size new rows).
        subunits_set (set): Accumulator of subunit identifiers seen so far.
        cont (int): Next free row index to assign to a newly-seen residue.
        matrix (list[list[str]]): The in-progress interaction matrix (rows = residues).
        aa (dict): Residue name -> row index, built up across files.
        template (list[dict]): Parsed contents of the template JSON file.
        interaction_list (list[str]): Interaction labels allowed by the template,
            in canonical order.
        protein (bool): Whether to include the protein atom in the cell's atom string.
        ligand (bool): Whether to include the ligand atom in the cell's atom string.
        subunit (bool): Whether to keep subunits distinct in the residue label.
        amino_acid_codes (list[str]): Standard three-letter amino acid codes,
            passed through to ``get_protein_ligand``.

    Returns:
        tuple[list[list[str]], str, dict, int, set]:
            Updated ``(matrix, ligand_code, aa, cont, subunits_set)``.
    """
    # Filter to obtain entries with the desired contact or type
    for inter in content:
        for entry in template:
            if _matches_template(inter, entry):
                contact = set()
                if isinstance(inter["type"], list):
                    for conta in inter["type"]:
                        if conta in interaction_list:
                            contact.add(conta)
                else:
                    if inter["type"] in interaction_list:
                        contact.add(inter["type"])
                if isinstance(inter["contact"], list):
                    for conta in inter["contact"]:
                        if conta in interaction_list:
                            contact.add(conta)
                else:
                    if inter["contact"] in interaction_list:
                        contact.add(inter["contact"])
                prot, lig = get_protein_ligand(begin=inter["bgn"], end=inter["end"], amino_acid_codes=amino_acid_codes)
                if prot is not None:
                    residue = prot["label_comp_id"] + " " + str(prot["auth_seq_id"])
                    prot_atom = prot["auth_atom_id"]
                    prot_subunit = prot["auth_asym_id"]
                    ligand_code = lig["label_comp_id"]
                    lig_atom = lig["auth_atom_id"]

                    subunits_set.add(prot_subunit)
                    atoms = (
                        f"{prot_atom}-{lig_atom}"
                        if protein and ligand
                        else prot_atom
                        if protein
                        else lig_atom
                        if ligand
                        else ""
                    )

                    if subunit:
                        residue += "-" + prot_subunit
                    else:
                        atoms += f"({prot_subunit})"

                    if residue not in aa:
                        aa[residue] = cont
                        cont += 1
                    column = aa[residue]

                    # Ensure matrix size and modify cell
                    if len(matrix) <= column:
                        matrix.append([""] * len(files))

                    for interaction in contact:
                        matrix[column][index] = modify_cell(
                            text=matrix[column][index],
                            interaction=interaction,
                            atoms=atoms,
                            interaction_labels=interaction_list,
                        )
                break

    return matrix, ligand_code, aa, cont, subunits_set


def parse_arpeggio_file_pdbe_template(
    records: list[tuple[dict, dict]],
    index: int,
    files: list[str],
    subunits_set: set,
    cont: int,
    matrix: list[list[str]],
    aa: dict,
    exclude_rules: list[dict],
    interaction_list: list[str],
    protein: bool,
    ligand: bool,
    subunit: bool,
) -> tuple[list[list[str]], str, dict, int, set]:
    """
    Processes the "processed"-schema PDBe-API Arpeggio export, restricted
    by a template's exclude rules and interaction vocabulary.

    Args:
        records (list[tuple[dict, dict]]): Flattened ``(record, ligand_context)``
            pairs. Records still following the raw schema (``is_raw_arpeggio_record``
            is True) are skipped — those are handled by ``parse_arpeggio_file_template``
            instead, since a directory can mix files of both schemas.
        index (int): Column index for this file within the overall matrix.
        files (list[str]): All filenames being processed (used only to size new rows).
        subunits_set (set): Accumulator of subunit identifiers seen so far.
        cont (int): Next free row index to assign to a newly-seen residue.
        matrix (list[list[str]]): The in-progress interaction matrix (rows = residues).
        aa (dict): Residue name -> row index, built up across files.
        exclude_rules (list[dict]): Rules (matched the same way as the raw
            template, via structural comparison) whose match causes a record
            to be dropped. Defaults to filtering out water
            (``end.chem_comp_id == "HOH"``) when the template's ``processed``
            section doesn't define this.
        interaction_list (list[str]): Interaction-type vocabulary and numeric-code
            order. Shared with the raw-schema template so raw and processed files
            land in the same code space when a directory mixes both.
        protein (bool): Whether to include the protein atom(s) in the cell's atom string.
        ligand (bool): Whether to include the ligand atom(s) in the cell's atom string.
        subunit (bool): Whether to keep subunits distinct in the residue label.

    Returns:
        tuple[list[list[str]], str, dict, int, set]:
            Updated ``(matrix, ligand_code, aa, cont, subunits_set)``.
    """
    ligand_code = None
    for record, ligand_context in records:
        if is_raw_arpeggio_record(record):
            continue

        if any(_matches_template(record, rule) for rule in exclude_rules):
            continue

        contact = set()
        details = record.get("interaction_details")
        if isinstance(details, list):
            contact.update(d for d in details if d in interaction_list)
        elif details in interaction_list:
            contact.add(details)

        itype = record.get("interaction_type")
        if isinstance(itype, list):
            contact.update(t for t in itype if t in interaction_list)
        elif itype in interaction_list:
            contact.add(itype)

        prot_end = record["end"]
        residue = f'{prot_end["chem_comp_id"]} {prot_end["author_residue_number"]}'
        prot_atom = ",".join(prot_end.get("atom_names", []))
        prot_subunit = prot_end["chain_id"]
        ligand_code = ligand_context["chem_comp_id"]
        lig_atom = ",".join(record.get("ligand_atoms", []))

        subunits_set.add(prot_subunit)
        atoms = (
            f"{prot_atom}-{lig_atom}"
            if protein and ligand
            else prot_atom
            if protein
            else lig_atom
            if ligand
            else ""
        )

        if subunit:
            residue += "-" + prot_subunit
        else:
            atoms += f"({prot_subunit})"

        if residue not in aa:
            aa[residue] = cont
            cont += 1
        column = aa[residue]

        # Ensure matrix size and modify cell
        if len(matrix) <= column:
            matrix.append([""] * len(files))

        for interaction in contact:
            matrix[column][index] = modify_cell(
                text=matrix[column][index],
                interaction=interaction,
                atoms=atoms,
                interaction_labels=interaction_list,
            )

    return matrix, ligand_code, aa, cont, subunits_set
