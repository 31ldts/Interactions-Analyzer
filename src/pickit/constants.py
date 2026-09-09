"""Global constants used across the pickit package.

Extracted verbatim from the original monolithic ``analyze_interactions.py``
(Fase 2, paso 1 del plan de modularización). No logic changes were made
during this extraction — only the physical location of the code changed.
"""

###########
# Globals #
###########

# Labels for different types of molecular interactions
INTERACTION_LABELS = [
    "Hydrophobic",
    "Aromatic_Face/Face",
    "Aromatic_Edge/Face",
    "HBond_PROT",
    "HBond_LIG",
    "Ionic_PROT",
    "Ionic_LIG",
    "Metal Acceptor",
    "Pi/Cation",
    "Other_Interactions",
]

# List of colors corresponding to different interaction types
COLORS = ["#ff6384", "#36a2eb", "#ffce56", "#4bc0c0", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

# List of available program modes
PROGRAM_MODES = ["ichem", "arpeggio"]

# Available modes for heatmap visualization
HEATMAP_MODES = ["max", "min", "mean", "count", "percent"]

# Arpeggio-specific interaction entities
ARPEGGIO_INT_ENT = ["INTER"]

# Types of molecular contacts identified by Arpeggio
ARPEGGIO_CONT = [
    "AMIDEAMIDE",
    "AMIDERING",
    "CARBONPI",
    "CATIONPI",
    "METSULPHURPI",
    "covalent",
    "hydrophobic",
    "carbonyl",
    "aromatic",
    "metal",
    "ionic",
    "hbond",
    "DONORPI",
    "weak_hbond",
    "polar",
    "weak_polar",
    "xbond",
    "HALOGENPI",
]

# Types of Arpeggio interactions (currently only plane-plane interactions)
ARPEGGIO_TYPE = ["plane-plane"]

# Colors associated with Arpeggio interaction types
ARPEGGIO_COLORS = [
    "#e6194B",
    "#3cb44b",
    "#ffe119",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#46f0f0",
    "#f032e6",
    "#bcf60c",
    "#fabebe",
    "#008080",
    "#e6beff",
    "#9a6324",
    "#fffac8",
    "#800000",
    "#aaffc3",
    "#808000",
    "#ffd8b1",
    "#000075",
    "#808080",
    "#7fffd4",
    "#ff7f50",
    "#6495ed",
    "#dc143c",
    "#9932cc",
    "#8b008b",
    "#556b2f",
    "#ff8c00",
    "#8fbc8f",
    "#483d8b",
    "#2f4f4f",
    "#00ced1",
    "#9400d3",
    "#ff1493",
    "#1e90ff",
    "#b22222",
    "#228b22",
    "#daa520",
    "#4b0082",
    "#cd5c5c",
]

# Constants for delimiters used in interaction data representation
SAME_DELIM = ", "  # Separates interactions of the same type.
DIFF_DELIM = "; "  # Separates interactions of different types.
GROUP_DELIM = "|"  # Groups interactions of the same type.

# Constants for handling empty cell values in interaction matrices
EMPTY_CELL = ""  # Represents an originally empty cell.
EMPTY_DASH_CELL = "-"  # Represents a cell that is considered empty after filtering.

# Global list of three-letter amino acid codes
AMINO_ACID_CODES = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]



DEFAULT_ARPEGGIO_TEMPLATE_FILE = "template_arpeggio.json"
DEFAULT_ICHEM_TEMPLATE_FILE = "template_ichem.json"

DEFAULT_TEMPLATE_FILE = DEFAULT_ARPEGGIO_TEMPLATE_FILE

###########
# Helpers #
###########


def is_not_empty_or_dash(cell: str) -> bool:
    """Return True unless `cell` is an empty cell or a post-filter dash cell."""
    return not (cell == EMPTY_DASH_CELL or cell == EMPTY_CELL)
