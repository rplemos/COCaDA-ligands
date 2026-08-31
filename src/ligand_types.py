"""
Atom property vectors for non-amino-acid components.

GENERATED FILE - do not edit by hand.
Produced by tools/generate_atom_types.py from LUNA.fdef via RDKit.
Regenerate rather than patching; hand edits are lost and break the
provenance chain the validation gate depends on.

Vector: [Hydrophobic, Aromatic, Positive, Negative, Donor, Acceptor]
"""

ligand_contact_types = {
    # ACT
    'ACT:C': [0, 0, 0, 0, 0, 0],
    'ACT:CH3': [1, 0, 0, 0, 0, 0],
    'ACT:O': [0, 0, 0, 1, 0, 1],
    'ACT:OXT': [0, 0, 0, 1, 0, 1],
    # GOL
    'GOL:C1': [0, 0, 0, 0, 0, 0],
    'GOL:C2': [0, 0, 0, 0, 0, 0],
    'GOL:C3': [0, 0, 0, 0, 0, 0],
    'GOL:O1': [0, 0, 0, 0, 1, 1],
    'GOL:O2': [0, 0, 0, 0, 1, 1],
    'GOL:O3': [0, 0, 0, 0, 1, 1],
    # HEM
    'HEM:C1A': [0, 1, 0, 0, 0, 0],
    'HEM:C1B': [0, 0, 0, 0, 0, 0],
    'HEM:C1C': [0, 1, 0, 0, 0, 0],
    'HEM:C1D': [0, 0, 0, 0, 0, 0],
    'HEM:C2A': [1, 1, 0, 0, 0, 0],
    'HEM:C2B': [1, 0, 0, 0, 0, 0],
    'HEM:C2C': [1, 1, 0, 0, 0, 0],
    'HEM:C2D': [1, 0, 0, 0, 0, 0],
    'HEM:C3A': [1, 1, 0, 0, 0, 0],
    'HEM:C3B': [1, 0, 0, 0, 0, 0],
    'HEM:C3C': [1, 1, 0, 0, 0, 0],
    'HEM:C3D': [1, 0, 0, 0, 0, 0],
    'HEM:C4A': [0, 1, 0, 0, 0, 0],
    'HEM:C4B': [0, 0, 0, 0, 0, 0],
    'HEM:C4C': [0, 1, 0, 0, 0, 0],
    'HEM:C4D': [0, 0, 0, 0, 0, 0],
    'HEM:CAA': [1, 0, 0, 0, 0, 0],
    'HEM:CAB': [1, 0, 0, 0, 0, 0],
    'HEM:CAC': [1, 0, 0, 0, 0, 0],
    'HEM:CAD': [1, 0, 0, 0, 0, 0],
    'HEM:CBA': [1, 0, 0, 0, 0, 0],
    'HEM:CBB': [1, 0, 0, 0, 0, 0],
    'HEM:CBC': [1, 0, 0, 0, 0, 0],
    'HEM:CBD': [1, 0, 0, 0, 0, 0],
    'HEM:CGA': [0, 0, 0, 0, 0, 0],
    'HEM:CGD': [0, 0, 0, 0, 0, 0],
    'HEM:CHA': [1, 0, 0, 0, 0, 0],
    'HEM:CHB': [1, 0, 0, 0, 0, 0],
    'HEM:CHC': [1, 0, 0, 0, 0, 0],
    'HEM:CHD': [1, 0, 0, 0, 0, 0],
    'HEM:CMA': [1, 0, 0, 0, 0, 0],
    'HEM:CMB': [1, 0, 0, 0, 0, 0],
    'HEM:CMC': [1, 0, 0, 0, 0, 0],
    'HEM:CMD': [1, 0, 0, 0, 0, 0],
    'HEM:FE': [0, 0, 0, 0, 0, 0],
    'HEM:NA': [0, 1, 0, 0, 0, 0],
    'HEM:NB': [0, 0, 0, 0, 0, 0],
    'HEM:NC': [0, 1, 0, 0, 0, 0],
    'HEM:ND': [0, 0, 0, 0, 0, 0],
    'HEM:O1A': [0, 0, 0, 1, 0, 1],
    'HEM:O1D': [0, 0, 0, 1, 0, 1],
    'HEM:O2A': [0, 0, 0, 1, 0, 1],
    'HEM:O2D': [0, 0, 0, 1, 0, 1],
    # SER
    'SER:C': [0, 0, 0, 0, 0, 0],
    'SER:CA': [0, 0, 0, 0, 0, 0],
    'SER:CB': [0, 0, 0, 0, 0, 0],
    'SER:N': [0, 0, 0, 0, 1, 0],
    'SER:O': [0, 0, 0, 0, 0, 1],
    'SER:OG': [0, 0, 0, 0, 1, 1],
    # T25
    'T25:C1': [0, 0, 0, 0, 0, 0],
    'T25:C10': [1, 0, 0, 0, 0, 0],
    'T25:C11': [1, 0, 0, 0, 0, 0],
    'T25:C12': [0, 0, 0, 0, 0, 0],
    'T25:C13': [0, 0, 0, 0, 0, 0],
    'T25:C14': [1, 0, 0, 0, 0, 0],
    'T25:C15': [1, 0, 0, 0, 0, 0],
    'T25:C16': [1, 0, 0, 0, 0, 0],
    'T25:C17': [1, 0, 0, 0, 0, 0],
    'T25:C18': [1, 0, 0, 0, 0, 0],
    'T25:C2': [1, 0, 0, 0, 0, 0],
    'T25:C3': [1, 0, 0, 0, 0, 0],
    'T25:C4': [1, 0, 0, 0, 0, 0],
    'T25:C5': [1, 0, 0, 0, 0, 0],
    'T25:C6': [1, 0, 0, 0, 0, 0],
    'T25:C7': [1, 0, 0, 0, 0, 0],
    'T25:C8': [1, 0, 0, 0, 0, 0],
    'T25:C9': [1, 0, 0, 0, 0, 0],
    'T25:O1': [0, 0, 0, 1, 0, 1],
    'T25:O2': [0, 0, 0, 1, 0, 1],
    'T25:O3': [0, 0, 0, 0, 0, 1],
}


# Rings available for aromatic stacking, as lists of atom names.
#
# A ring qualifies when it has five or six atoms and all of them are
# aromatic or sp2 - planar and conjugated. Filtering on aromaticity alone
# would give heme only two of its four pyrroles, because both the CCD and
# RDKit Kekulise rings B and D as pyrrolenines.
#
# Centroids are per ring, not per fused system: the centroid of a large
# fused plane sits too far from its own edges for a 2-5 A window.
ligand_rings = {
    'HEM': [
        ['C1A', 'C2A', 'C3A', 'C4A', 'NA'],
        ['C1B', 'C2B', 'C3B', 'C4B', 'NB'],
        ['C1C', 'C2C', 'C3C', 'C4C', 'NC'],
        ['C1D', 'C2D', 'C3D', 'C4D', 'ND'],
    ],
}
