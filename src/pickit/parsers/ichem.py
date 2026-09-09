"""IChem interaction file parser.

Extracted from the ``ichem_analysis`` closure nested inside the original
``analyze_files`` (see ``parsers/__init__.py`` for context). No logic
changes: ``protein``, ``ligand``, ``subunit`` were captured from the
enclosing method call via closure and are now explicit parameters instead.
"""

from ..constants import GROUP_DELIM, INTERACTION_LABELS
from .common import modify_cell, validate_string


def parse_ichem_file(
    content: list[str],
    index: int,
    files: list[str],
    subunits_set: set,
    cont: int,
    matrix: list[list[str]],
    aa: dict,
    protein: bool,
    ligand: bool,
    subunit: bool,
    interaction_list: list[str] | None = None,
) -> tuple[list[list[str]], dict, int, set]:
    """
    Processes an IChem interaction file's content, extracting relevant data
    and updating the (in-progress) interaction matrix.

    Args:
        content (list[str]): Lines of the IChem interaction file.
        index (int): Column index for this file within the overall matrix.
        files (list[str]): All filenames being processed (used only to size new rows).
        subunits_set (set): Accumulator of subunit identifiers seen so far.
        cont (int): Next free row index to assign to a newly-seen residue.
        matrix (list[list[str]]): The in-progress interaction matrix (rows = residues).
        aa (dict): Residue name -> row index, built up across files.
        protein (bool): Whether to include the protein atom in the cell's atom string.
        ligand (bool): Whether to include the ligand atom in the cell's atom string.
        subunit (bool): Whether to keep subunits distinct in the residue label.
        interaction_list (list[str] | None): Optional IChem template — an
            ordered subset of ``INTERACTION_LABELS`` to restrict parsing to
            (lines whose interaction isn't in this list are skipped) and to
            use as the numeric-code order passed to ``modify_cell``. When
            ``None`` (the default, and the pre-template behavior), every
            line is kept and codes follow ``INTERACTION_LABELS``.

    Returns:
        tuple[list[list[str]], dict, int, set]: Updated ``(matrix, aa, cont, subunits_set)``.
    """
    labels = interaction_list if interaction_list is not None else INTERACTION_LABELS
    for line in content:
        elements = line.split(GROUP_DELIM)
        if len(elements) == 10:
            interaction = elements[0].strip().replace("\t", "")
            if interaction_list is not None and interaction not in interaction_list:
                continue
            residue = elements[3].strip().replace("\t", "")
            if validate_string(residue):
                if not subunit:
                    sections = residue.split("-")
                    residue = sections[0]
                    subunits_set.add(sections[1])

                atoms = (
                    f"{elements[1].strip()}-{elements[4].strip()}"
                    if protein and ligand
                    else elements[1].strip()
                    if protein
                    else elements[4].strip()
                    if ligand
                    else ""
                )
                if not subunit:
                    atoms += f"({sections[1]})"

                if residue not in aa:
                    aa[residue] = cont
                    cont += 1
                column = aa[residue]

                # Ensure matrix size and modify cell
                if len(matrix) <= column:
                    matrix.append([""] * len(files))

                matrix[column][index] = modify_cell(
                    text=matrix[column][index],
                    interaction=interaction,
                    atoms=atoms,
                    interaction_labels=labels,
                )
    return matrix, aa, cont, subunits_set
