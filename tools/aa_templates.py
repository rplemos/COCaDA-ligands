"""
Capped templates for the 20 standard amino acids, used to validate the
atom-type generator against COCaDA's curated `contact_types` table.

Why capped: an amino acid inside a chain has its backbone nitrogen and
carbonyl engaged in peptide bonds. Typing the free amino acid (which is what
the PDB Chemical Component Dictionary stores) would give a free amine and a
free carboxylate, and the resulting donor/acceptor/charge assignments would
not match the residue as it actually occurs in a protein. Each residue is
therefore built as its N-acetyl / N'-methylamide analogue (ACE-X-NME).

Why atom maps: RDKit gives no PDB atom names, so each atom we care about
carries a map number, and `MAPS` translates map number -> PDB atom name.

IMPORTANT - bracketed SMILES atoms carry NO implicit hydrogens. Writing
`[N:1]` instead of `[NH:1]` silently strips the backbone amide proton and
destroys every donor assignment downstream. Every bracketed atom below
states its hydrogen count explicitly.

COCaDA-ligands (v2.0)
"""

# Backbone scaffold. "{R}" is replaced by the residue-specific fragment,
# which must itself supply the mapped alpha carbon.
#   map 1 = N    map 2 = CA    map 3 = C    map 4 = O
BACKBONE = "CC(=O)[NH:1]{R}[C:3](=[O:4])NC"

_BB = {1: "N", 2: "CA", 3: "C", 4: "O"}


def _m(**side):
    """Merge the backbone map with a residue's side-chain map."""
    out = dict(_BB)
    out.update({int(k[1:]): v for k, v in side.items()})
    return out


# fragment, {map number: PDB atom name}
TEMPLATES = {
    "A": ("[C@@H:2]([CH3:5])",
          _m(m5="CB")),

    "R": ("[C@@H:2]([CH2:5][CH2:6][CH2:7][NH:8][C:9](=[NH2+:10])[NH2:11])",
          _m(m5="CB", m6="CG", m7="CD", m8="NE", m9="CZ", m10="NH1", m11="NH2")),

    "N": ("[C@@H:2]([CH2:5][C:6](=[O:7])[NH2:8])",
          _m(m5="CB", m6="CG", m7="OD1", m8="ND2")),

    "D": ("[C@@H:2]([CH2:5][C:6](=[O:7])[O-:8])",
          _m(m5="CB", m6="CG", m7="OD1", m8="OD2")),

    "C": ("[C@@H:2]([CH2:5][SH:6])",
          _m(m5="CB", m6="SG")),

    "Q": ("[C@@H:2]([CH2:5][CH2:6][C:7](=[O:8])[NH2:9])",
          _m(m5="CB", m6="CG", m7="CD", m8="OE1", m9="NE2")),

    "E": ("[C@@H:2]([CH2:5][CH2:6][C:7](=[O:8])[O-:9])",
          _m(m5="CB", m6="CG", m7="CD", m8="OE1", m9="OE2")),

    "G": ("[CH2:2]",
          _m()),

    # Neutral N-delta-1 tautomer (HID), the dominant form at pH 7.4.
    # COCaDA marks BOTH ring nitrogens positive; that is handled by a charge
    # placement rule, not by protonating the template.
    "H": ("[C@@H:2]([CH2:5][c:6]1[nH:7][cH:9][n:10][cH:8]1)",
          _m(m5="CB", m6="CG", m7="ND1", m8="CD2", m9="CE1", m10="NE2")),

    "I": ("[C@@H:2]([C@H:5]([CH2:6][CH3:8])[CH3:7])",
          _m(m5="CB", m6="CG1", m7="CG2", m8="CD1")),

    "L": ("[C@@H:2]([CH2:5][CH:6]([CH3:7])[CH3:8])",
          _m(m5="CB", m6="CG", m7="CD1", m8="CD2")),

    "K": ("[C@@H:2]([CH2:5][CH2:6][CH2:7][CH2:8][NH3+:9])",
          _m(m5="CB", m6="CG", m7="CD", m8="CE", m9="NZ")),

    "M": ("[C@@H:2]([CH2:5][CH2:6][S:7][CH3:8])",
          _m(m5="CB", m6="CG", m7="SD", m8="CE")),

    "F": ("[C@@H:2]([CH2:5][c:6]1[cH:7][cH:9][cH:11][cH:10][cH:8]1)",
          _m(m5="CB", m6="CG", m7="CD1", m8="CD2", m9="CE1", m10="CE2", m11="CZ")),

    # Proline: the backbone nitrogen is part of the ring and is a tertiary
    # amide, so it carries no hydrogen and is correctly not a donor.
    "P": ("RING",
          _m(m5="CB", m6="CG", m7="CD")),

    "S": ("[C@@H:2]([CH2:5][OH:6])",
          _m(m5="CB", m6="OG")),

    "T": ("[C@@H:2]([C@H:5]([OH:6])[CH3:7])",
          _m(m5="CB", m6="OG1", m7="CG2")),

    "W": ("[C@@H:2]([CH2:5][c:6]1[cH:7][nH:9][c:10]2[cH:12][cH:14][cH:13][cH:11][c:8]12)",
          _m(m5="CB", m6="CG", m7="CD1", m8="CD2", m9="NE1",
             m10="CE2", m11="CE3", m12="CZ2", m13="CZ3", m14="CH2")),

    "Y": ("[C@@H:2]([CH2:5][c:6]1[cH:7][cH:9][c:11]([OH:12])[cH:10][cH:8]1)",
          _m(m5="CB", m6="CG", m7="CD1", m8="CD2", m9="CE1",
             m10="CE2", m11="CZ", m12="OH")),

    "V": ("[C@@H:2]([CH:5]([CH3:6])[CH3:7])",
          _m(m5="CB", m6="CG1", m7="CG2")),
}

# Proline needs the nitrogen inside the ring, so it cannot use BACKBONE.
PROLINE_SMILES = "CC(=O)[N:1]1[CH2:7][CH2:6][CH2:5][C@H:2]1[C:3](=[O:4])NC"


def smiles_for(code):
    """Return (SMILES, {map number: PDB atom name}) for a one-letter code."""
    frag, mapping = TEMPLATES[code]
    if frag == "RING":
        return PROLINE_SMILES, mapping
    return BACKBONE.replace("{R}", frag), mapping


ALL_CODES = sorted(TEMPLATES)
