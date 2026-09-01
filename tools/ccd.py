"""
ccd.py - read PDB Chemical Component Dictionary entries and build RDKit
molecules from them, preserving PDB atom names.

The atom-type generator needs three things that a SMILES string cannot supply:

  1. PDB ATOM NAMES. COCaDA keys its property table on `COMP:ATOM`, and those
     names come from the depositor's file. SMILES has no notion of them.
  2. AUTHORITATIVE BONDS. Structure files carry no connectivity, so the CCD is
     the only place the bond orders of a HETATM group can come from.
  3. AUTHORITATIVE AROMATICITY. RDKit's perception of a hand-written macrocycle
     SMILES is unreliable - heme is the standard failure - so ring assignment
     is taken from the CCD's own `pdbx_aromatic_flag` instead of re-perceived.

The CCD ships as one enormous `components.cif` (all ~48 000 components) or as
per-component files from `files.rcsb.org/ligands/`. Both are accepted; the
reader splits on `data_` blocks either way.


THINGS THE REAL DATA FORCES US TO HANDLE
----------------------------------------
These are not hypothetical. Each was found in a real CCD entry.

* METAL COORDINATION. Heme records `FE` bonded to all four pyrrole nitrogens
  with `pdbx_metal_coordination_flag = Y`. Treating those as ordinary single
  bonds overruns nitrogen's valence and the molecule fails to sanitise, so
  they are built as RDKit DATIVE bonds directed from the donor to the metal.
  That flag is also the input for a future metal-coordination contact type.

* INTERNALLY INCONSISTENT CHARGE AND PROTONATION. Heme gives `O2A` a formal
  charge of -1 *and* a bonded hydrogen `H2A`, which is not a valid oxygen.
  `repair_charged_hydrogens` drops the hydrogen and keeps the charge, since
  the deprotonated carboxylate is the chemically intended state. Every repair
  is recorded on the component so it can be reported rather than hidden.

* PARTIAL AROMATICITY. Heme flags rings A and C aromatic but not B and D,
  which is a faithful record of porphyrin's ambiguous Kekule structure rather
  than an error. The reader passes the flags through unchanged; downstream
  code must not assume a symmetric macrocycle.

* LEAVING ATOMS. Polymer components are stored as the FREE amino acid, so SER
  carries `OXT`, `HXT` and a second amide hydrogen `H2`, all flagged
  `pdbx_leaving_atom_flag = Y`. Those atoms are absent from a residue inside a
  chain. Dropping them and capping the open valences reconstructs the residue
  as it actually occurs - which is what the curated amino-acid table describes.

COCaDA v1.6+ - ligand support extension
"""

import os
import re

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

NULL_TOKENS = {"?", "."}

#: CCD `value_order` -> RDKit bond type.
BOND_ORDERS = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "QUAD": Chem.BondType.QUADRUPLE,
    "AROM": Chem.BondType.AROMATIC,
}


# --------------------------------------------------------------- tokenising

_TOKEN = re.compile(r"'([^']*)'|\"([^\"]*)\"|(\S+)")


def _tokenise(line):
    """Split one mmCIF data line, honouring single and double quoting."""
    out = []
    for single, double, bare in _TOKEN.findall(line):
        out.append(single or double or bare)
    return out


def _parse_block(lines):
    """Parse one data_ block into {category: [row dicts]} plus scalar items."""
    categories = {}
    scalars = {}
    index = 0
    total = len(lines)

    while index < total:
        line = lines[index].strip()

        if not line or line.startswith("#"):
            index += 1
            continue

        if line == "loop_":
            index += 1
            headers = []
            while index < total and lines[index].strip().startswith("_"):
                headers.append(lines[index].strip())
                index += 1
            if not headers:
                continue
            category = headers[0].split(".")[0]
            columns = [h.split(".", 1)[1] for h in headers]

            rows = []
            while index < total:
                raw = lines[index]
                stripped = raw.strip()
                if not stripped or stripped.startswith("#") or stripped == "loop_" \
                        or stripped.startswith("_") or stripped.startswith("data_"):
                    break
                values = _tokenise(stripped)
                # A row may wrap; keep consuming until the column count is met.
                while len(values) < len(columns) and index + 1 < total:
                    index += 1
                    values.extend(_tokenise(lines[index].strip()))
                if len(values) == len(columns):
                    rows.append(dict(zip(columns, values)))
                index += 1
            categories.setdefault(category, []).extend(rows)
            continue

        if line.startswith("_"):
            parts = _tokenise(line)
            if len(parts) >= 2:
                scalars[parts[0]] = " ".join(parts[1:])
            index += 1
            continue

        index += 1

    return categories, scalars


# --------------------------------------------------------------- components

class Component:
    """One CCD chemical component."""

    def __init__(self, comp_id, atoms, bonds, scalars=None):
        self.id = comp_id
        self.atoms = atoms          # list of _chem_comp_atom row dicts
        self.bonds = bonds          # list of _chem_comp_bond row dicts
        self.scalars = scalars or {}
        self.repairs = []           # human-readable notes, see repair_* below

    # -- convenience accessors -------------------------------------------

    @property
    def name(self):
        return self.scalars.get("_chem_comp.name", "").strip('"')

    @property
    def type(self):
        return self.scalars.get("_chem_comp.type", "").strip('"')

    @property
    def is_polymer(self):
        """True for components that occur as residues inside a chain."""
        kind = self.type.upper()
        return "PEPTIDE" in kind or "DNA" in kind or "RNA" in kind or "SACCHARIDE" in kind

    def heavy_atom_names(self, drop_leaving=False):
        return [a["atom_id"] for a in self.atoms
                if a.get("type_symbol", "").upper() != "H"
                and not (drop_leaving and a.get("pdbx_leaving_atom_flag") == "Y")]

    def __repr__(self):
        return (f"<Component {self.id} "
                f"{len(self.heavy_atom_names())} heavy atoms, {len(self.bonds)} bonds>")


def parse(text):
    """Parse CCD text (one component or the full components.cif).

    Returns {comp_id: Component}.
    """
    components = {}
    current_id = None
    buffer = []

    def flush():
        if current_id is None:
            return
        categories, scalars = _parse_block(buffer)
        atoms = categories.get("_chem_comp_atom", [])
        bonds = categories.get("_chem_comp_bond", [])
        if atoms:
            components[current_id] = Component(current_id, atoms, bonds, scalars)

    for line in text.splitlines():
        if line.startswith("data_"):
            flush()
            current_id = line[len("data_"):].strip()
            buffer = []
        else:
            buffer.append(line)
    flush()
    return components


def iter_components(path):
    """Yield Components one at a time from a CCD file.

    The full `components.cif` is over half a gigabyte and holds ~51 000
    entries, so it is streamed a data block at a time rather than read into
    memory. Only one component is materialised at once.
    """

    if not os.path.exists(path):
        raise FileNotFoundError(f"CCD file not found: {path}")

    def build(comp_id, buffer):
        categories, scalars = _parse_block(buffer)
        atoms = categories.get("_chem_comp_atom", [])
        if atoms:
            return Component(comp_id, atoms,
                             categories.get("_chem_comp_bond", []), scalars)
        return None

    current_id, buffer = None, []
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("data_"):
                if current_id is not None:
                    component = build(current_id, buffer)
                    if component is not None:
                        yield component
                current_id = line[len("data_"):].strip()
                buffer = []
            else:
                buffer.append(line.rstrip("\n"))

    if current_id is not None:
        component = build(current_id, buffer)
        if component is not None:
            yield component


def load(path):
    """Parse a CCD file from disk into {comp_id: Component}.

    Convenient for single-component files. Use `iter_components` for the full
    dictionary, which is far too large to hold at once.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"CCD file not found: {path}")
    with open(path, errors="replace") as handle:
        return parse(handle.read())


def load_dir(directory):
    """Parse every .cif in a directory of per-component files."""
    components = {}
    for name in sorted(os.listdir(directory)):
        if name.endswith(".cif"):
            components.update(load(os.path.join(directory, name)))
    return components


# ------------------------------------------------------------------ repairs

def repair_charged_hydrogens(mol, component):
    """Drop hydrogens bonded to a formally negative heteroatom.

    Real case: heme records O2A with charge -1 and a bonded hydrogen H2A. The
    deprotonated carboxylate is the intended species, so the charge is kept
    and the hydrogen removed. Without this the molecule will not sanitise.
    """
    doomed = []
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() >= 0 or atom.GetSymbol() == "H":
            continue
        for neighbour in atom.GetNeighbors():
            if neighbour.GetSymbol() == "H":
                doomed.append(neighbour.GetIdx())
                component.repairs.append(
                    f"removed {neighbour.GetProp('atom_name')} bonded to "
                    f"{atom.GetProp('atom_name')} (formal charge "
                    f"{atom.GetFormalCharge()})")
    for idx in sorted(doomed, reverse=True):
        mol.RemoveAtom(idx)
    return mol


# ------------------------------------------------------------ mol building

def to_rdkit(component, drop_leaving=False, cap_peptide=False):
    """Build an RDKit molecule from a component.

    Args:
        component: a Component.
        drop_leaving: remove atoms flagged `pdbx_leaving_atom_flag = Y`. For a
            polymer component this converts the stored free monomer into the
            residue as it occurs inside a chain.
        cap_peptide: for amino-acid components, replace the removed leaving
            atoms with an acetyl on N and an N-methylamide on C, so the
            backbone amide chemistry matches a residue in a chain. Requires
            drop_leaving.

    Returns:
        (mol, {atom index: PDB atom name}). Capping atoms map to None so they
        are typed but never emitted.
    """
    mol = Chem.RWMol()
    index_of = {}
    names = {}

    skip = set()
    if drop_leaving:
        skip = {a["atom_id"] for a in component.atoms
                if a.get("pdbx_leaving_atom_flag") == "Y"}

    for row in component.atoms:
        atom_id = row["atom_id"]
        if atom_id in skip:
            continue
        symbol = row.get("type_symbol", "C").capitalize()
        atom = Chem.Atom(symbol)

        charge = row.get("charge", "0")
        atom.SetFormalCharge(0 if charge in NULL_TOKENS else int(float(charge)))

        if row.get("pdbx_aromatic_flag") == "Y":
            atom.SetIsAromatic(True)

        atom.SetNoImplicit(True)          # the CCD lists hydrogens explicitly
        atom.SetProp("atom_name", atom_id)
        idx = mol.AddAtom(atom)
        index_of[atom_id] = idx
        names[idx] = atom_id

    for row in component.bonds:
        first, second = row.get("atom_id_1"), row.get("atom_id_2")
        if first in skip or second in skip:
            continue
        if first not in index_of or second not in index_of:
            continue

        if row.get("pdbx_metal_coordination_flag") == "Y":
            # Direction matters for a dative bond: donor -> metal.
            begin, end = index_of[first], index_of[second]
            if mol.GetAtomWithIdx(begin).GetSymbol() not in ("N", "O", "S", "P"):
                begin, end = end, begin
            mol.AddBond(begin, end, Chem.BondType.DATIVE)
            continue

        order = BOND_ORDERS.get(row.get("value_order", "SING"), Chem.BondType.SINGLE)
        if row.get("pdbx_aromatic_flag") == "Y":
            order = Chem.BondType.AROMATIC
        mol.AddBond(index_of[first], index_of[second], order)

    repair_charged_hydrogens(mol, component)

    if cap_peptide:
        if not drop_leaving:
            raise ValueError("cap_peptide requires drop_leaving=True")
        _cap_peptide(mol, component)

    names = {a.GetIdx(): (a.GetProp("atom_name") if a.HasProp("atom_name") else None)
             for a in mol.GetAtoms()}

    final = mol.GetMol()
    try:
        Chem.SanitizeMol(final)
    except Exception as exc:                                  # noqa: BLE001
        component.repairs.append(f"full sanitisation failed ({exc}); "
                                 f"retried without kekulisation")
        final = mol.GetMol()
        Chem.SanitizeMol(
            final,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
            catchErrors=True)

    # A partially sanitised molecule can be left with uninitialised ring info,
    # and GetRingInfo() then trips a C++ assertion that takes the whole process
    # down rather than raising. FastFindRings initialises it without requiring
    # a full sanitisation, which matters when walking all ~51 000 components.
    try:
        final.GetRingInfo().NumRings()
    except Exception:                                         # noqa: BLE001
        Chem.FastFindRings(final)

    names = {a.GetIdx(): (a.GetProp("atom_name") if a.HasProp("atom_name") else None)
             for a in final.GetAtoms()}
    return final, names


def _cap_peptide(mol, component):
    """Add ACE on the backbone N and NME on the backbone C.

    The CCD marks backbone atoms explicitly (`pdbx_backbone_atom_flag`), so the
    two attachment points are found from the data rather than by name-guessing,
    which keeps this working for modified residues too.
    """
    backbone = {a["atom_id"] for a in component.atoms
                if a.get("pdbx_backbone_atom_flag") == "Y"}
    if not backbone:
        return

    by_name = {a.GetProp("atom_name"): a.GetIdx()
               for a in mol.GetAtoms() if a.HasProp("atom_name")}

    n_idx = by_name.get("N")
    c_idx = by_name.get("C")

    def add(symbol, charge=0):
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(False)
        atom.SetFormalCharge(charge)
        return mol.AddAtom(atom)

    if n_idx is not None and "N" in backbone:            # acetyl cap
        carbonyl = add("C")
        oxygen = add("O")
        methyl = add("C")
        mol.AddBond(n_idx, carbonyl, Chem.BondType.SINGLE)
        mol.AddBond(carbonyl, oxygen, Chem.BondType.DOUBLE)
        mol.AddBond(carbonyl, methyl, Chem.BondType.SINGLE)

    if c_idx is not None and "C" in backbone:            # N-methylamide cap
        nitrogen = add("N")
        methyl = add("C")
        mol.AddBond(c_idx, nitrogen, Chem.BondType.SINGLE)
        mol.AddBond(nitrogen, methyl, Chem.BondType.SINGLE)
