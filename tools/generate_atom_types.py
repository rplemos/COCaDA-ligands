#!/usr/bin/env python3
"""
generate_atom_types.py - offline generator for COCaDA atom property tables.

COCaDA classifies every atom with a six-element binary vector

    [ Hydrophobic, Aromatic, Positive, Negative, Donor, Acceptor ]

which for the twenty standard amino acids is hand-curated in
`src/conditions.py` (167 entries). Extending contact detection to prosthetic
groups, cofactors, ions and ligands means producing the same vector for
arbitrary molecules, and hand curation does not scale: heme alone is 43
atoms and the PDB Chemical Component Dictionary holds ~48 000 components.

This script derives those vectors automatically from molecular structure,
using the LUNA pharmacophore feature definitions (`LUNA.fdef`, MIT licence,
vendored alongside this file) evaluated through RDKit. It runs OFFLINE and
emits a static Python data file. COCaDA's runtime never imports RDKit and
keeps its numpy + psutil dependency set.

    Why offline: PDB and mmCIF atom records carry neither bonds nor
    hydrogens, both of which SMARTS-based typing requires, so runtime
    perception would need a dictionary lookup regardless. A committed data
    file is also reproducible and reviewable in a way a perception step
    running against whatever RDKit the user has installed is not.


VALIDATION GATE
---------------
LUNA's definitions reproduce COCaDA's curated table for 161 of 167 atoms
(96.4%) before any adjustment. All six differences are one systematic
convention: LUNA assigns ionizable-group features to every atom of the
matched group, COCaDA to the chemically charged atoms only. Two
reconciliation rules (below) close the gap to 167/167.

That exact round trip is the acceptance test, and it is the whole warrant
for trusting generated types on molecules with no reference data: there is
no independent way to validate the output for heme, so the evidence has to
be that the same generator reproduces chemistry already validated and
published. Run with --validate; a non-zero exit means do not ship the table.


RECONCILIATION RULES
--------------------
R1  Charge placement.
    An ionizable group feature is assigned to the group's PERIPHERAL
    HETEROATOMS - that is, excluding both its carbons and its central hub
    atom (one bonded to three or more other atoms of the same group).
      - carboxylate  -> the two oxygens, not the carboxyl carbon
      - imidazolium  -> the two ring nitrogens, not the ring carbons
      - phosphate,
        sulfate      -> the oxygens, not the central P or S
    Exception - guanidinium, where the central carbon IS included, because
    the positive charge is genuinely delocalised over CZ by resonance. This
    reproduces COCaDA's curation of R:CZ as positive.

    The hub clause matters for ligands rather than for amino acids: no
    standard residue has a tetrahedral oxyanion, but ATP, NAD, FAD, PO4 and
    SO4 all do, and marking the buried central atom as charged would emit
    electrostatic contacts from an atom shielded by its own oxygens.

R2  Thioether sulfur.
    A sulfur bonded to two carbons is typed acceptor, not hydrophobic.
    NOTE - this one follows COCaDA against the evidence. Both LUNA and
    pmapper independently type thioether sulfur as hydrophobic and the
    chemistry favours them, but changing it would move published
    amino-acid results. The divergence is deliberate, frozen, and belongs
    in the methods section. See decision D6 in the design record.

Both rules are general and apply to ligands as well as amino acids, so the
same chemistry is used throughout the output.


USAGE
-----
    python tools/generate_atom_types.py --validate
    python tools/generate_atom_types.py --components HEM ATP ZN
    python tools/generate_atom_types.py --components HEM --out src/ligand_types.py

COCaDA v1.6+ - ligand support extension
"""

import argparse
import os
import sys

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
FDEF = os.path.join(HERE, "LUNA.fdef")

# ---------------------------------------------------------------- flag model

#: Index of each property in COCaDA's six-element vector.
FLAG_ORDER = ["Hydrophobic", "Aromatic", "Positive", "Negative", "Donor", "Acceptor"]
FLAG_SHORT = ["H", "a", "P", "N", "D", "A"]

#: LUNA feature family -> index in the COCaDA vector. Families absent from
#: this map (WeakDonor, ChalcogenDonor, HalogenDonor, Metal, Amide,
#: Electrophile, Nucleophile, Atom) have no COCaDA equivalent today. Metal
#: and HalogenDonor are the natural next contact types; see the design record.
FAMILY_TO_FLAG = {
    "Hydrophobic": 0,
    "Aromatic": 1,
    "PositivelyIonizable": 2,
    "NegativelyIonizable": 3,
    "Donor": 4,
    "Acceptor": 5,
}

IONIZABLE_FAMILIES = ("PositivelyIonizable", "NegativelyIonizable")

# ------------------------------------------------------- reconciliation rules

#: R1 exception - the guanidinium carbon carries delocalised positive charge.
SMARTS_GUANIDINIUM = Chem.MolFromSmarts("[#7][CX3](=[#7])[#7]")

#: R2 - sulfur bonded to two carbons and nothing else.
SMARTS_THIOETHER = Chem.MolFromSmarts("[SX2]([#6])[#6]")


def _load_factory(fdef=FDEF):
    if not os.path.exists(fdef):
        sys.exit(f"error: feature definitions not found at {fdef}")
    return ChemicalFeatures.BuildFeatureFactory(fdef)


#: A group atom bonded to at least this many other atoms of the same group
#: is treated as the group's central hub (the P of a phosphate, the S of a
#: sulfate) and does not itself carry the charge.
HUB_DEGREE = 3


def _charge_bearing_atoms(mol_h, group_ids, guanidinium_atoms):
    """R1 - which atoms of an ionizable group actually carry the charge.

    Peripheral heteroatoms carry it. Carbons do not, unless the carbon is a
    guanidinium centre. Central hub heteroatoms do not either, so a
    phosphate or sulfate charges its oxygens rather than its P or S.
    """
    group = set(group_ids)
    keep = []
    for idx in group_ids:
        atom = mol_h.GetAtomWithIdx(idx)
        if idx in guanidinium_atoms:
            keep.append(idx)                      # resonance exception
            continue
        if atom.GetSymbol() == "C":
            continue                              # carbons never carry it
        in_group_degree = sum(
            1 for nb in atom.GetNeighbors() if nb.GetIdx() in group)
        if in_group_degree >= HUB_DEGREE:
            continue                              # central hub atom
        keep.append(idx)
    # A group that loses every atom would silently drop its charge; fall
    # back to the raw group rather than emitting nothing.
    return keep or list(group_ids)


def type_molecule(mol, factory=None, name_for=None):
    """Assign COCaDA property vectors to the atoms of an RDKit molecule.

    Args:
        mol: RDKit Mol. Bonds and formal charges must already be correct;
            hydrogens are added internally.
        factory: an RDKit feature factory. Built from LUNA.fdef if omitted.
        name_for: callable taking a heavy-atom index and returning that
            atom's name, or None to exclude it from the output. Defaults to
            the atom's index as a string.

    Returns:
        dict mapping atom name -> six-element list of 0/1.
    """
    factory = factory or _load_factory()
    if name_for is None:
        name_for = str

    mol_h = Chem.AddHs(mol)

    guanidinium = set()
    for match in mol_h.GetSubstructMatches(SMARTS_GUANIDINIUM):
        guanidinium.update(match)

    result = {}

    def mark(idx, flag):
        name = name_for(idx)
        if name is None:
            return
        result.setdefault(name, [0] * 6)[flag] = 1

    for feature in factory.GetFeaturesForMol(mol_h):
        family = feature.GetFamily()
        flag = FAMILY_TO_FLAG.get(family)
        if flag is None:
            continue
        ids = feature.GetAtomIds()
        if family in IONIZABLE_FAMILIES:                       # R1
            ids = _charge_bearing_atoms(mol_h, ids, guanidinium)
        for idx in ids:
            mark(idx, flag)

    for match in mol_h.GetSubstructMatches(SMARTS_THIOETHER):  # R2
        sulfur = match[0]
        name = name_for(sulfur)
        if name is None:
            continue
        vec = result.setdefault(name, [0] * 6)
        vec[0] = 0   # not hydrophobic
        vec[5] = 1   # acceptor

    # Atoms matched by no feature still need an all-zero entry, so that a
    # missing key means "unknown component" rather than "apolar atom".
    for atom in mol.GetAtoms():
        name = name_for(atom.GetIdx())
        if name is not None:
            result.setdefault(name, [0] * 6)

    return result


def format_vector(vec):
    """Render a property vector as 'HaPNDA', dots for unset flags."""
    return "".join(FLAG_SHORT[i] if vec[i] else "." for i in range(6))


# ------------------------------------------------------------- validation

def validate(verbose=True):
    """Regenerate the twenty standard amino acids and diff against the
    curated table. Returns the list of disagreements (empty means pass)."""
    sys.path.insert(0, HERE)
    sys.path.insert(0, os.path.join(HERE, os.pardir, "src"))
    import aa_templates as templates

    try:
        from conditions import contact_types as curated
    except ImportError:
        try:
            from conditions_ref import contact_types as curated
        except ImportError:
            sys.exit("error: cannot import the curated contact_types table")

    factory = _load_factory()
    diffs = []
    checked = 0

    for code in templates.ALL_CODES:
        smiles, mapping = templates.smiles_for(code)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            diffs.append((f"{code}:*", None, None, "SMILES failed to parse"))
            continue

        mol_h = Chem.AddHs(mol)
        idx_to_map = {a.GetIdx(): a.GetAtomMapNum() for a in mol_h.GetAtoms()}

        def name_for(idx, _m=idx_to_map, _map=mapping):
            return _map.get(_m.get(idx, 0))

        generated = type_molecule(mol, factory, name_for)

        for atom_name in mapping.values():
            key = f"{code}:{atom_name}"
            if key not in curated:
                continue
            checked += 1
            got = generated.get(atom_name, [0] * 6)
            expected = curated[key]
            if got != expected:
                diffs.append((key, expected, got, ""))

    if verbose:
        passed = checked - len(diffs)
        print(f"round-trip against curated contact_types: "
              f"{passed}/{checked} atoms reproduced exactly")
        if diffs:
            print("\nDISAGREEMENTS (curated vs generated):")
            for key, exp, got, note in diffs:
                if note:
                    print(f"   {key:<10} {note}")
                else:
                    print(f"   {key:<10} curated {format_vector(exp)}   "
                          f"generated {format_vector(got)}")
            print("\nGATE FAILED - do not ship this table.")
        else:
            print("GATE PASSED - generator reproduces the curated table exactly.")

    return diffs


# --------------------------------------------------------------- components

#: Minimal built-in component set so the generator can be exercised before a
#: full Chemical Component Dictionary ingest exists. Real runs must read the
#: CCD, which supplies authoritative atom NAMES, bond orders and aromatic
#: flags. SMILES alone gives no PDB atom names, and RDKit's aromaticity
#: perception of hand-written macrocycle SMILES is unreliable - heme is the
#: known problem case. See the design record.
BUILTIN_COMPONENTS = {
    "HOH": "O",
    "ZN":  "[Zn+2]",
    "MG":  "[Mg+2]",
    "NA":  "[Na+]",
    "CL":  "[Cl-]",
    "CA":  "[Ca+2]",
    "SO4": "[O-]S(=O)(=O)[O-]",
    "PO4": "[O-]P(=O)([O-])[O-]",
    "ACT": "CC(=O)[O-]",
    "EDO": "OCCO",
    "GOL": "OCC(O)CO",
}


def generate_components(names, factory=None):
    """Type a set of built-in components. Returns {comp_id: {atom: vector}}."""
    factory = factory or _load_factory()
    out = {}
    for comp_id in names:
        smiles = BUILTIN_COMPONENTS.get(comp_id)
        if smiles is None:
            print(f"warning: {comp_id} not in the built-in set; "
                  f"CCD ingest not yet implemented", file=sys.stderr)
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"warning: {comp_id} SMILES failed to parse", file=sys.stderr)
            continue
        # Without CCD atom names, fall back to element+ordinal (O1, O2, ...).
        counts = {}
        names_by_idx = {}
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol().upper()
            counts[sym] = counts.get(sym, 0) + 1
            names_by_idx[atom.GetIdx()] = (
                sym if mol.GetNumAtoms() == 1 else f"{sym}{counts[sym]}")
        out[comp_id] = type_molecule(mol, factory, lambda i: names_by_idx.get(i))
    return out


def write_table(table, path):
    """Emit a generated table as an importable Python module."""
    lines = [
        '"""',
        "Atom property vectors for non-amino-acid components.",
        "",
        "GENERATED FILE - do not edit by hand.",
        "Produced by tools/generate_atom_types.py from LUNA.fdef via RDKit.",
        "Regenerate rather than patching; hand edits are lost and break the",
        "provenance chain the validation gate depends on.",
        "",
        "Vector: [Hydrophobic, Aromatic, Positive, Negative, Donor, Acceptor]",
        '"""',
        "",
        "ligand_contact_types = {",
    ]
    for comp_id in sorted(table):
        lines.append(f"    # {comp_id}")
        for atom_name in sorted(table[comp_id]):
            vec = table[comp_id][atom_name]
            lines.append(f"    '{comp_id}:{atom_name}': {vec},")
    lines.append("}")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"wrote {sum(len(v) for v in table.values())} entries "
          f"across {len(table)} components to {path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate COCaDA atom property vectors from molecular structure.")
    parser.add_argument("--validate", action="store_true",
                        help="Reproduce the 20 standard amino acids and diff "
                             "against the curated table. Exits non-zero on any "
                             "disagreement. This is the acceptance gate.")
    parser.add_argument("--components", nargs="+", metavar="ID",
                        help="Component identifiers to type.")
    parser.add_argument("--out", metavar="PATH",
                        help="Write the generated table to this Python file.")
    args = parser.parse_args()

    if not args.validate and not args.components:
        parser.print_help()
        return 0

    if args.validate:
        if validate():
            return 1

    if args.components:
        table = generate_components(args.components)
        if args.out:
            write_table(table, args.out)
        else:
            for comp_id in sorted(table):
                print(f"\n{comp_id}")
                for atom_name in sorted(table[comp_id]):
                    print(f"   {atom_name:<6} {format_vector(table[comp_id][atom_name])}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
