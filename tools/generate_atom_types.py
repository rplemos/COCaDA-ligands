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

    Bridge clause (2026-10-01): a group oxygen or sulfur bonded to two or more
    heavy atoms is an ester or anhydride bridge and is excluded too, so ATP
    charges its seven terminal phosphate oxygens and not O3A, O3B or O5'.
    No standard residue has such a group, so the gate is unaffected.

R2  Thioether sulfur.
    A sulfur bonded to two carbons is typed acceptor, not hydrophobic.
    NOTE - this one follows COCaDA against the evidence. Both LUNA and
    pmapper independently type thioether sulfur as hydrophobic and the
    chemistry favours them, but changing it would move published
    amino-acid results. The divergence is deliberate, frozen, and belongs
    in the methods section. See decision D6 in the design record.

R3  Metal-coordinated donors are not free ions.
    An atom holding a dative bond to a metal has its Positive and Negative
    flags cleared. Its charge is engaged in coordinating the metal and is
    not available to form a salt bridge with a distant residue.
    Found from real data: the CCD gives heme's NA and NC a formal charge of
    -1, correctly, since the porphyrin is a dianion. Without this rule those
    two nitrogens - buried in the macrocycle and bonded to the iron - would
    emit salt bridges and attractive contacts with any nearby lysine or
    arginine, which is chemically wrong.
    No standard amino acid has a metal coordination bond, so this rule
    cannot affect the validation gate; it applies to ligands only.

R4  Ring tautomerism in five-membered rings only. (2026-10-01)
    LUNA marks a ring nitrogen 1,3 to an N-H as both donor and acceptor,
    since the hydrogen may sit on either: right for imidazole, triazole and
    pyrazole, and what makes histidine's ND1 and NE2 reproduce the curated
    table. In six-membered rings the same SMARTS fire on lactams whose minor
    tautomer is negligible - guanine N3 and flavin N1, which carry no
    hydrogen, were donors, and the N-H of guanine N1 and uracil N3 were
    acceptors. Found by comparing against RING's independent donor/acceptor
    assignments, ~3,800 atoms. In a five-membered ring the features hold as
    written; elsewhere the atom's own hydrogens decide - donor only with an
    H, acceptor only without. Not "drop them": LUNA's ordinary donor rule
    excludes acidic N-H, so an imide N-H (uracil N3, flavin N3) is a donor
    only through the tautomer feature, and would otherwise donate nothing.

All four rules are general and apply to ligands as well as amino acids, so
the same chemistry is used throughout the output.


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

#: R4 - LUNA's tautomer features, honoured as such in five-membered rings only.
TAUTOMER_FEATURES = ("TautomerRingDonorN", "TautomerRingAcceptorN")


def _tautomer_applies(atom, family):
    """R4 - whether a tautomer donor/acceptor feature holds for this atom.

    In a five-membered ring (imidazole, triazole, pyrazole) the hydrogen really
    may sit on either nitrogen, and the feature holds as LUNA writes it. In any
    other ring the atom's own structure decides: the donor half only on an N
    that carries a hydrogen, the acceptor half only on one that does not.
    """
    if atom.IsInRingSize(5):
        return True
    has_h = atom.GetTotalNumHs(includeNeighbors=True) > 0
    return has_h if family == "Donor" else not has_h

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

#: A group oxygen or sulfur bonded to two or more heavy atoms is a BRIDGE - the
#: ester O5' of a nucleotide, the anhydride O3A/O3B of ATP, the ester oxygen of
#: a sulfate ester - and does not carry the charge either. LUNA's phosphate and
#: sulfate SMARTS match those atoms as group members (`P(=O)(-O)(-O)-[O-]`), so
#: without this they were typed Negative, formed salt bridges, and took the
#: electrostatic registry's one slot per residue pair from the real anionic
#: oxygens: 1G78's Arg67 salt-bridged FMN at O5' (3.85 A) instead of O1P (2.63 A).
BRIDGE_ELEMENTS = ("O", "S")


def _charge_bearing_atoms(mol_h, group_ids, guanidinium_atoms, bridges=True):
    """R1 - which atoms of an ionizable group actually carry the charge.

    Peripheral heteroatoms carry it. Carbons do not, unless the carbon is a
    guanidinium centre. Central hub heteroatoms do not either, so a
    phosphate or sulfate charges its oxygens rather than its P or S, and nor
    do ester or anhydride bridges between two heavy atoms.
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
        if bridges and atom.GetSymbol() in BRIDGE_ELEMENTS and sum(
                1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1) >= 2:
            continue                              # ester / anhydride bridge
        keep.append(idx)
    # A group that loses every atom would silently drop its charge. The bridge
    # clause can do that on its own - an anionic oxygen that also binds a metal
    # has two heavy neighbours - so first retry without it, which reproduces the
    # rule as it stood before the clause; only then fall back to the raw group.
    if not keep and bridges:
        return _charge_bearing_atoms(mol_h, group_ids, guanidinium_atoms,
                                     bridges=False)
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
        if feature.GetType() in TAUTOMER_FEATURES:             # R4
            ids = [i for i in ids if _tautomer_applies(mol_h.GetAtomWithIdx(i), family)]
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

    for atom in mol_h.GetAtoms():                              # R3
        if not any(b.GetBondType() == Chem.BondType.DATIVE for b in atom.GetBonds()):
            continue
        name = name_for(atom.GetIdx())
        if name is None or name not in result:
            continue
        result[name][2] = 0   # not a free cation
        result[name][3] = 0   # not a free anion

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


def stacking_rings(mol, name_of):
    """Rings of a molecule that can take part in aromatic stacking.

    A ring qualifies when it has five or six atoms and every one of them is
    aromatic or sp2 - that is, the ring is planar and conjugated, which is what
    stacking actually requires.

    Why not "all atoms aromatic". Both the CCD and RDKit write heme with only
    two of its four pyrroles aromatic; rings B and D are Kekulised as
    pyrrolenines. Filtering on aromaticity therefore yields an asymmetric two
    of four, which is an artifact of the Kekule structure rather than
    chemistry. The planar-conjugated test returns all four, symmetrically,
    with no per-component curation. Verified to reject sp3 rings: cyclohexane,
    pyrrolidine and T25's epoxide all yield nothing.

    Why individual rings rather than fused systems. Stacking is scored on
    centroid-to-centroid distance within a 2-5 A window, and the centroid of a
    large fused system sits too far from its own edges. Measured on 2RCL, the
    nearest aromatic residue is 7.46 A from heme's whole-plane centroid but
    5.81 A from the closest pyrrole - a 1.6 A difference that decides whether
    real stacking is found. This also matches LUNA, whose aromatic SMARTS
    match individual five- and six-membered rings.

    Note the resulting asymmetry with amino acids: tryptophan keeps its
    historical treatment as a single fused nine-atom ring, because changing it
    would move published results. Ligands use per-ring centroids.
    """

    from rdkit.Chem import HybridizationType

    try:
        atom_rings = mol.GetRingInfo().AtomRings()
    except Exception:                                         # noqa: BLE001
        from rdkit import Chem as _Chem
        _Chem.FastFindRings(mol)
        atom_rings = mol.GetRingInfo().AtomRings()

    rings = []
    for ring in atom_rings:
        if not 5 <= len(ring) <= 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetIsAromatic() or a.GetHybridization() == HybridizationType.SP2
                   for a in atoms):
            continue
        names = sorted(n for n in (name_of(i) for i in ring) if n)
        if len(names) == len(ring):
            rings.append(names)
    return rings


def generate_from_ccd(path, wanted=None, factory=None, verbose=True,
                      start=0, limit=None):
    """Type components read from the Chemical Component Dictionary.

    Args:
        path: a components.cif, a single-component .cif, or a directory of them.
        wanted: component ids to keep, or None for all.
        factory: an RDKit feature factory.
        verbose: report per-component repairs and failures.

    Returns:
        ({comp_id: {atom_name: vector}}, [notes])

    Amino-acid components are typed as they occur INSIDE A CHAIN: leaving
    atoms are dropped and the backbone is capped. The CCD stores the free
    monomer, whose free amine and free carboxylate would otherwise produce
    donor, acceptor and charge assignments that no residue in a protein
    actually has. Sugars and nucleotides are typed as free molecules; see
    ccd.Component.typed_in_chain for why.
    """
    import ccd as ccd_module

    factory = factory or _load_factory()

    if os.path.isdir(path):
        source = sorted(ccd_module.load_dir(path).items())
    elif os.path.getsize(path) > 50 * 1024 * 1024:
        # The full components.cif is streamed: half a gigabyte, ~51 000 entries.
        source = ((c.id, c) for c in ccd_module.iter_components(path))
    else:
        source = sorted(ccd_module.load(path).items())

    keep = set(wanted) if wanted else None

    table, rings, notes = {}, {}, []
    done = 0
    seen = 0
    for comp_id, component in source:
        if keep is not None and comp_id not in keep:
            continue
        seen += 1
        if seen <= start:
            continue
        if limit is not None and done >= limit:
            break
        try:
            mol, names = ccd_module.for_typing(component)
        except Exception as exc:                              # noqa: BLE001
            notes.append(f"{comp_id}: could not build molecule ({exc})")
            continue

        try:
            typed = type_molecule(mol, factory, lambda i: names.get(i))
        except Exception as exc:                              # noqa: BLE001
            notes.append(f"{comp_id}: typing failed ({exc})")
            continue

        heavy = {}
        for atom in mol.GetAtoms():
            name = names.get(atom.GetIdx())
            if name and atom.GetSymbol() != "H":
                heavy[name] = typed.get(name, [0] * 6)
        table[comp_id] = heavy

        try:
            systems = stacking_rings(mol, lambda i: names.get(i))
        except Exception as exc:                              # noqa: BLE001
            notes.append(f"{comp_id}: ring perception failed ({exc})")
            systems = []
        if systems:
            rings[comp_id] = systems

        done += 1
        if verbose and done % 5000 == 0:
            print(f"  ... {done:,} components typed", file=sys.stderr, flush=True)

        for repair in component.repairs:
            notes.append(f"{comp_id}: {repair}")
        if verbose and component.repairs:
            for repair in component.repairs:
                print(f"  note  {comp_id}: {repair}", file=sys.stderr)

    return table, rings, notes


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


def write_dat(table, rings, path):
    """Emit the compact, gzipped, line-per-component data file.

    A Python dict literal is not an option at CCD scale: 51 000 components and
    1.33 million heavy atoms would be a ~50 MB source file re-parsed by every
    worker process at import. This format is read lazily by src/ligand_data.py,
    which loads only the components present in the structure being processed.

    See that module for the format description.
    """

    import gzip
    sys.path.insert(0, os.path.join(HERE, os.pardir, "src"))
    from ligand_data import pack

    written = 0
    with gzip.open(path, "wt", compresslevel=9) as handle:
        for comp_id in sorted(table):
            atoms = table[comp_id]
            if not atoms:
                continue
            packed = ",".join(f"{name}={pack(atoms[name])}"
                              for name in sorted(atoms))
            systems = ";".join("|".join(ring) for ring in rings.get(comp_id, []))
            handle.write(f"{comp_id}\t{packed}\t{systems}\n")
            written += 1

    size = os.path.getsize(path) / (1024 * 1024)
    print(f"wrote {written:,} components, "
          f"{sum(len(v) for v in table.values()):,} atom entries, "
          f"{sum(len(v) for v in rings.values()):,} ring systems "
          f"-> {path} ({size:.1f} MB gzipped)")


def write_table(table, path, rings=None):
    """Emit a generated table as an importable Python module."""
    rings = rings or {}
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

    lines += [
        "",
        "",
        "# Rings available for aromatic stacking, as lists of atom names.",
        "#",
        "# A ring qualifies when it has five or six atoms and all of them are",
        "# aromatic or sp2 - planar and conjugated. Filtering on aromaticity alone",
        "# would give heme only two of its four pyrroles, because both the CCD and",
        "# RDKit Kekulise rings B and D as pyrrolenines.",
        "#",
        "# Centroids are per ring, not per fused system: the centroid of a large",
        "# fused plane sits too far from its own edges for a 2-5 A window.",
        "ligand_rings = {",
    ]
    for comp_id in sorted(rings):
        lines.append(f"    '{comp_id}': [")
        for system in rings[comp_id]:
            lines.append(f"        {system},")
        lines.append("    ],")
    lines.append("}")

    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print(f"wrote {sum(len(v) for v in table.values())} entries "
          f"across {len(table)} components to {path}"
          + (f"; {sum(len(v) for v in rings.values())} aromatic ring system(s) "
             f"across {len(rings)} component(s)" if rings else ""))


def main():
    parser = argparse.ArgumentParser(
        description="Generate COCaDA atom property vectors from molecular structure.")
    parser.add_argument("--validate", action="store_true",
                        help="Reproduce the 20 standard amino acids and diff "
                             "against the curated table. Exits non-zero on any "
                             "disagreement. This is the acceptance gate.")
    parser.add_argument("--components", nargs="+", metavar="ID",
                        help="Component identifiers to type. With --ccd, selects "
                             "from the dictionary; without it, uses the small "
                             "built-in set.")
    parser.add_argument("--ccd", metavar="PATH",
                        help="Chemical Component Dictionary source: components.cif, "
                             "a single-component .cif, or a directory of them. This "
                             "is the authoritative source of atom names, bond orders "
                             "and aromatic flags.")
    parser.add_argument("--start", type=int, default=0, metavar="N",
                        help="Skip the first N components. With --limit, lets the "
                             "full dictionary be generated in chunks.")
    parser.add_argument("--limit", type=int, default=None, metavar="N",
                        help="Stop after N components.")
    parser.add_argument("--out", metavar="PATH",
                        help="Write the generated table to this Python file.")
    args = parser.parse_args()

    if not args.validate and not args.components and not args.ccd:
        parser.print_help()
        return 0

    if args.validate:
        if validate():
            return 1

    if args.ccd:
        table, rings, notes = generate_from_ccd(args.ccd, args.components,
                                                start=args.start, limit=args.limit)
        if args.out:
            if args.out.endswith(".gz"):
                write_dat(table, rings, args.out)
            else:
                write_table(table, args.out, rings)
        else:
            for comp_id in sorted(table):
                print(f"\n{comp_id}")
                for atom_name in sorted(table[comp_id]):
                    print(f"   {atom_name:<6} {format_vector(table[comp_id][atom_name])}")
        if notes:
            print(f"\n{len(notes)} note(s) recorded during generation", file=sys.stderr)
    elif args.components:
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
