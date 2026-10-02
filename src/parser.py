"""
COCaDA v1.6
Author: Rafael Lemos - rafaellemos@ufmg.br
Date: 26/03/2026

License: MIT License
"""

from src.classes import Protein, Chain, Residue, Atom

import os
from numpy import mean, array
from numpy.linalg import svd
import re


stacking = {
    'H':[10, 'CG','ND1','CE1','NE2','CD2'],
    'F':[11, 'CG','CD1','CE1','CZ','CE2','CD2'],
    'W':[14, 'CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
    'Y':[12, 'CG','CD1','CE1','CZ','CE2','CD2'],
}

residue_mapping = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# HETATM component ids that are never treated as contact-forming entities.
# Water dominates by count - a large structure carries over a thousand copies,
# so water-water alone would be ~10^6 atom pairs, dwarfing the protein - and
# its positions are the least reliable part of most models. Water-mediated
# bridges are a separate feature, not a side effect of admitting HETATM.
SKIP_HETERO = {'HOH', 'DOD', 'WAT', 'H2O'}

# ... with one exception. A water that COORDINATES A METAL is not incidental
# solvent, it is a ligand of the metal centre: the sixth coordination position
# of a resting-state cytochrome P450 is an aqua ligand, and PDB 1TQN models it
# 2.07 A from the heme iron. Excluding water dropped that bond entirely.
#
# Such waters are admitted, and ONLY they: a water is kept when it lies within
# the metal-coordination window of a metal atom, which is one distance test per
# water and nothing at all in a structure with no metals. They carry an
# all-zero property vector, exactly as metals do, so the only contact they can
# form is the coordination itself - never a hydrogen bond, never a hydrophobic
# contact. Water-mediated bridges remain a separate feature.
#
# 33 of 84 metal atoms in the validation benchmark (39%) have a water within
# the window: 12 of 12 magnesium, 8 of 10 copper, 2 of 16 iron.
WATERS = {'HOH', 'DOD', 'WAT', 'H2O'}

# Modified residues that are deposited as HETATM but are chemically standard
# amino acids, and are mapped onto their parent so they use the Ca path.
# MSE (selenomethionine) is by far the most common residue in this class.
HETERO_AS_STANDARD = {
    'MSE': 'MET',
}

# ... and the atoms that carry a different name in the modified residue. Mapping
# MSE onto MET made its backbone and carbons Met atoms, but its selenium is
# deposited as SE where Met has SD, and `M:SE` has no entry in the curated
# table - so the one atom that distinguishes the residue formed no contacts.
# Renamed to SD, it is typed as the thioether it substitutes for (acceptor) and
# reported as Met, consistently with the residue name, which already prints as
# M. The element stays SE: nothing in contact typing reads it for this atom.
STANDARD_ATOM_NAMES = {
    'MSE': {'SE': 'SD'},
}

# FORCE-FIELD RESIDUE NAMES (D32). Simulation and modelling pipelines name a
# residue's protonation state in its residue name - Amber's HID, HIE, HIP, ASH,
# GLH, LYN, CYX, CYM; CHARMM's HSD, HSE, HSP, LSN, ASPP, GLUP; GROMACS's HISD,
# HISE, HISH, ASPH, GLUH, LYSH, CYS2. Most of them are ALSO Chemical Component
# Dictionary ids, of unrelated molecules (ASH a kinase inhibitor, HIP
# ND1-phosphonohistidine, HSE homoserine, CYX a phosphonic acid), and deposited
# entries use them as such - HIP inside 18 entries' chains, HSE inside 8 and as
# a ligand in 16.
#
# So a residue is read as its native amino acid only when every one of its atoms
# is an atom of that amino acid, CA included. The Amber ASH of a model passes; a
# deposited phosphohistidine HIP (it has a P), homoserine HSE (C3, O3), S-methyl-
# cysteine CYM (CD) or lysine amide LYN (NT) does not, and stays the component it
# is. Chain membership is not required: models write caps, and sometimes
# residues, as separate non-polymer entities. Decided per
# residue, once parsing is done (native_force_field_residues). The protonation
# the name encodes is not kept: COCaDA assigns it from the pH, as for any residue.
FORCE_FIELD_NAMES = {
    'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'ASH': 'ASP', 'GLH': 'GLU',
    'LYN': 'LYS', 'CYX': 'CYS', 'CYM': 'CYS',                          # Amber
    'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS', 'LSN': 'LYS', 'ASPP': 'ASP',
    'GLUP': 'GLU',                                                      # CHARMM
    'HISD': 'HIS', 'HISE': 'HIS', 'HISH': 'HIS', 'ASPH': 'ASP', 'GLUH': 'GLU',
    'LYSH': 'LYS', 'CYS2': 'CYS',                                       # GROMACS
}

# ... and the terminal caps, which stay components but are named, or their atoms
# named, as the Chemical Component Dictionary has them: Amber's NMA is the CCD's
# NME (methylamine; NMA there is a PLP derivative), whose carbon the CCD calls C
# where Amber writes CH3. ACE needs no renaming; it is listed so that a cap
# written as an ATOM record is kept rather than dropped.
FORCE_FIELD_CAPS = {
    'NMA': ('NME', {'N': 'N', 'C': 'C', 'CA': 'C', 'CH3': 'C'}),
    'NME': ('NME', {'N': 'N', 'C': 'C', 'CA': 'C', 'CH3': 'C'}),
    'ACE': ('ACE', {'C': 'C', 'O': 'O', 'CH3': 'CH3'}),
}


def native_force_field_residues(protein):
    """Read force-field-named residues as what they are (D32).

    A residue the parser kept because of its name (see FORCE_FIELD_NAMES) is
    converted when all its atoms belong to the native residue; otherwise it is treated as any residue of that name: a
    HETATM component stays one, an ATOM record with a non-standard name is
    dropped, as it always was. A converted amino acid is kept even when HETATM
    are suppressed, as MSE is, since it is an amino acid either way.
    """
    from src.conditions import contact_types
    native_atoms = {}
    for key in contact_types:
        code, atom = key.split(":")
        native_atoms.setdefault(code, set()).add(atom)

    for chain in protein.chains:
        kept = []
        for residue in chain.residues:
            info = getattr(residue, "force_field", None)
            if info is None:
                kept.append(residue)
                continue
            atom_record = info
            names = {a.atomname for a in residue.atoms}
            source = residue.resname
            if source in FORCE_FIELD_NAMES:
                code = residue_mapping[FORCE_FIELD_NAMES[source]]
                if "CA" in names and names <= native_atoms[code] | {"OXT"}:
                    residue.source_name = source
                    residue.resname = code
                    residue.is_hetero = False
                    residue.atoms[:] = [a for a in residue.atoms if a.atomname != "OXT"]
                    add_standard_ring(residue)
                    kept.append(residue)
                    continue
            elif source in FORCE_FIELD_CAPS:
                cap, renames = FORCE_FIELD_CAPS[source]
                if names <= set(renames):
                    if cap in SKIP_HETERO:
                        continue            # components suppressed (validator)
                    residue.source_name = source
                    residue.resname = cap
                    for atom in residue.atoms:
                        atom.atomname = renames[atom.atomname]
                    kept.append(residue)
                    continue
            if atom_record or source in SKIP_HETERO:
                continue                    # as any residue of that name
            kept.append(residue)
        chain.residues = kept
    protein.chains = [chain for chain in protein.chains if chain.residues]
    return protein


def add_standard_ring(residue):
    """The ring centroid of a standard aromatic residue, as the parser adds it
    while reading - for a residue that became standard after reading."""
    spec = stacking.get(residue.resname)
    if not spec or residue.ring:
        return
    allowed = spec[1:]
    if (all(a.occupancy == 1 for a in residue.atoms if a.atomname in allowed)
            and len({a.atomname for a in residue.atoms}) == spec[0]):
        ring_atoms = array([[a.x, a.y, a.z] for a in residue.atoms if a.atomname in spec])
        if ring_atoms.any():
            residue.atoms.append(centroid(residue, ring_atoms, residue.atoms[0].entity))
            residue.ring = True
            residue.normal_vector = calc_normal_vector(ring_atoms)

# Nucleic acids are deferred: their values are to be imported from Pro-RNA and
# RNApedia, not re-derived here. The standard nucleotides are ATOM records and
# were always skipped - but a MODIFIED nucleotide in the same strand (5CM, BRU,
# 8OG) is a HETATM, and was admitted as a component and typed by the generator:
# half a strand, typed by a method the deferral exists to avoid. Such residues
# are now dropped too, until the import lands.
#
# The test is membership in a nucleic-acid CHAIN, not the component's CCD type:
# the CCD classes free nucleotide ligands as "RNA LINKING" too - GDP, UDP, A3P -
# and those are cofactors, which stay. mmCIF gives each free ligand its own
# entity, so a residue is in a strand when its entity also holds standard
# nucleotides. PDB format shares chain letters between polymer and ligands, so
# there it is a residue of a nucleic chain that SEQRES lists as part of it.
STANDARD_NUCLEOTIDES = {'A', 'C', 'G', 'U', 'I', 'N',
                        'DA', 'DC', 'DG', 'DT', 'DU', 'DI', 'DN'}


def drop_strand_components(protein, in_strand):
    """Remove the hetero residues `in_strand(residue)` places in a nucleic-acid
    chain, and any chain left empty by it. See STANDARD_NUCLEOTIDES."""
    for chain in protein.chains:
        chain.residues = [r for r in chain.residues
                          if not (r.is_hetero and in_strand(r))]
    protein.chains = [chain for chain in protein.chains if chain.residues]
    return protein


#: One mmCIF value: single-quoted, double-quoted or bare. A closing quote counts
#: only when whitespace or the end of the line follows it, which is the CIF rule
#: and the reason `"C1'"` is one value holding a prime rather than two.
_CIF_TOKEN = re.compile(r"""'(.*?)'(?=\s|$)|"(.*?)"(?=\s|$)|(\S+)""")


def cif_tokens(line):
    """Split one mmCIF data line into its values, removing the quotes.

    mmCIF quotes any value that contains a quote character, so every atom name
    carrying a prime is written `"C1'"`, `"O5'"`. A plain whitespace split kept
    the quotes as part of the name; the generated table holds `C1'`, so the
    lookup failed and every such atom formed no contacts at all - the ribose of
    ATP, NAD, FAD, SAH and every nucleotide cofactor or modified nucleotide.
    PDB format has no quoting, so the two formats disagreed on the same entry.

    Most lines carry no quote, and they keep the plain split, which is the
    fast path the parser has always used.
    """
    if '"' not in line and "'" not in line:
        return line.split()
    return [match.group(match.lastindex) for match in _CIF_TOKEN.finditer(line)]


#: ALTLOCS. By default each atom keeps ONE alternate - the highest occupancy,
#: the first on a tie (D6). That loses a contact the other conformer makes: in
#: 4D7B only Thr119 altloc B places OG1 2.52 A from tolcapone's carbonyl, and
#: the tie keeps A, whose CG2 sits there. `--all-altlocs` keeps every alternate
#: and reports each conformer's contacts, labelled; the default is unchanged.

#: Identifiers that stand for "no identifier". RCSB writes `_entry.id XXXX` in
#: every biological-assembly file, and other tools write `?` or `.`; taken as
#: the structure's id, every such file named its output `XXXX_contacts.csv`
#: and each run overwrote the last. When the file offers one of these, the id
#: stays the file name it was given at the start of parsing.
PLACEHOLDER_IDS = {"XXXX", "?", ".", ""}

#: An extended PDB id, `pdb_` and eight characters, of which the first four are
#: zeros for every id issued so far: `pdb_00001abc` is `1abc`.
EXTENDED_PDB_ID = re.compile(r"^pdb_0000([0-9a-z]{4})$", re.IGNORECASE)


def structure_id(entry_id):
    """The name a structure's output is written under, from its `_entry.id`.

    A PDB id is used as it is, and an extended PDB id as the four-character id
    it extends. Anything else is used whole. The last four characters used to be
    taken from every id, which is the PDB id for a deposit but not for a model:
    `AF-P02766-F1` became `6-F1`, so 489 AlphaFold models wrote to 12 output
    names and each batch kept about one CSV in ten (D31).
    """
    match = EXTENDED_PDB_ID.match(entry_id)
    if match:
        return match.group(1)
    return entry_id.replace("/", "_").replace(os.sep, "_")


#: COVALENT LINKS. Structure files carry no bonds inside a residue, but the
#: deposit declares the ones BETWEEN residues - `_struct_conn` in mmCIF, `LINK`
#: in PDB format: an N-glycosylation, a heme c thioether, a covalent inhibitor.
#: Atoms one or two bonds apart across such a link sit 1.3-2.5 A from each
#: other, and COCaDA used to report them as contacts (2GUY: Asn197 ND2 - NAG O5,
#: 2.39 A, "hydrogen bond"). The parser records the covalent links; contacts.py
#: suppresses the pairs around each one (D30). Only `covale*` records: `disulf`
#: is the DS contact type itself and `metalc` is metal coordination, both
#: reported on purpose. A link to a symmetry mate is not inside the file.
IDENTITY_OPERATORS = {"1_555", "1555", "?", ".", ""}


def link_partner(chain, resnum, icode, resname, atomname):
    """One end of a declared link, keyed as the parser keys residues - chain,
    number, insertion code, mapped residue name - with the atom name the parser
    gives it. Returns None when the number is not an integer."""
    try:
        resnum = int(resnum)
    except (TypeError, ValueError):
        return None
    icode = "" if icode in (None, "?", ".") else icode.strip()
    atomname = STANDARD_ATOM_NAMES.get(resname, {}).get(atomname, atomname)
    resname = HETERO_AS_STANDARD.get(resname, resname)
    resname = residue_mapping.get(resname, resname)
    return (chain, resnum, icode, resname), atomname


def struct_conn_links(rows):
    """The covalent links among `_struct_conn` rows (dicts keyed by item)."""
    links = []
    for row in rows:
        if not row.get("conn_type_id", "").startswith("covale"):
            continue
        if (row.get("ptnr1_symmetry", "1_555") not in IDENTITY_OPERATORS
                or row.get("ptnr2_symmetry", "1_555") not in IDENTITY_OPERATORS):
            continue
        ends = []
        for n in (1, 2):
            chain = row.get(f"ptnr{n}_label_asym_id", ".")
            if chain in (".", "?"):
                chain = row.get(f"ptnr{n}_auth_asym_id")
            number = row.get(f"ptnr{n}_auth_seq_id", row.get(f"ptnr{n}_label_seq_id"))
            ends.append(link_partner(chain, number, row.get(f"pdbx_ptnr{n}_PDB_ins_code"),
                                     row.get(f"ptnr{n}_label_comp_id"),
                                     row.get(f"ptnr{n}_label_atom_id")))
        if None not in ends:
            links.append(tuple(ends))
    return links


def pdb_link(line):
    """The covalent link on a PDB `LINK` record, or None.

    PDB format writes metal coordination on LINK records too; those are told
    apart by element when the link is applied (contacts.covalent_exclusions).
    """
    rec = line.ljust(80)
    if rec[59:65].strip() not in IDENTITY_OPERATORS or rec[66:72].strip() not in IDENTITY_OPERATORS:
        return None
    ends = (link_partner(rec[21], rec[22:26], rec[26], rec[17:20].strip(), rec[12:16].strip()),
            link_partner(rec[51], rec[52:56], rec[56], rec[47:50].strip(), rec[42:46].strip()))
    return None if None in ends else ends


def element_of(atomname, explicit=None):
    """Best available element symbol for an atom.

    The explicit element column is preferred wherever the format provides it
    (PDB columns 77-78, mmCIF `type_symbol`). Falling back to the atom name is
    ambiguous: in PDB format an alpha carbon is written ' CA ' and a calcium
    ion 'CA  ', and the parser strips whitespace, so the two become
    indistinguishable once HETATM records are admitted.
    """

    if explicit:
        return explicit.strip().upper()
    name = atomname.strip()
    return name[0].upper() if name else ''


def parse_pdb(pdb_file, all_altlocs=False):
    """
    Parses a PDB file to create a Protein object.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        Protein: The Protein object populated with chain, residue, and atom objects.

    This function reads a PDB file line by line, extracts information about chains, residues,
    and atoms, and constructs a Protein object. It handles special cases such as alternative
    residue names, aromatic residues, and low-quality atoms. The function also assigns a unique
    identifier and title to the protein based on the header information.
    """
    
    current_protein = Protein()
    current_chain = None
    current_residue = None
    current_entity = None
    entity_chains = {}
    entity = None
    waters = []              # held aside; see WATERS
    seqres = {}              # chain -> residue names, for strand membership
    links = []               # declared covalent links; see IDENTITY_OPERATORS
    ph = 7.4
    ph_pattern = re.compile(r'\bPH\b\s*[:\s]\s*([-+]?\d*\.\d+|\d+)')

    with open(pdb_file) as f:
        
        current_protein.id = os.path.basename(pdb_file).split(".")[0]
        
        for line in f:
            line = line.strip()
            
            if line == "ENDMDL":
                break
            
            # for interface checking
            elif line.startswith("COMPND"):
                if "MOL_ID" in line:
                    current_entity = line[-2]
                elif "CHAIN:" in line:
                    chains = line.split(":")[1].strip().replace(";","").replace(" ","")
                    entity_chains[current_entity] = chains.split(",")
        
            elif line.startswith("SEQRES"):
                seqres.setdefault(line[11:12], set()).update(line[19:].split())

            elif line.startswith("LINK "):
                link = pdb_link(line)
                if link:
                    links.append(link)

            elif line.startswith("HEADER"):
                if line[62:].strip().upper() not in PLACEHOLDER_IDS:
                    current_protein.id = line[62:]
                
            elif line.startswith("TITLE"):
                current_protein.set_title(line[10:])
                
            # remark 200 = x-ray; remark 210,215,217 = NMR    
            elif line.startswith("REMARK 200") or line.startswith("REMARK 21"):
                match = ph_pattern.search(line)
                if match:
                    ph_str = match.group(1)
                    if '-' in ph_str or '/' in ph_str or 'NULL' in line.upper():
                        continue
                    ph = float(ph_str)                
                
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                is_hetero_line = line.startswith("HETATM")
                chain_id = line[21]
                
                if entity_chains:
                    for key, chains in entity_chains.items():
                        if chain_id in chains:
                            entity = key
                            break
                        
                resnum = int(line[22:26])
                # Column 27 is the insertion code. Residues 100, 100A and 100B
                # are three residues on one number - antibody CDRs are numbered
                # exactly that way - so it is part of the residue's identity.
                icode = line[26:27].strip() if len(line) > 26 else ""
                # if resnum <= 0:
                #     continue
                resname = line[17:20].strip()
                
                # force-field names: kept, decided per residue at the end (D32)
                force_field = resname in FORCE_FIELD_NAMES or resname in FORCE_FIELD_CAPS

                # modified residues that are chemically standard amino acids
                renames = STANDARD_ATOM_NAMES.get(resname)
                resname = HETERO_AS_STANDARD.get(resname, resname)

                if resname in residue_mapping:
                    resname = residue_mapping[resname]   # one-letter code
                    is_hetero = False
                elif force_field or (is_hetero_line and resname not in SKIP_HETERO):
                    is_hetero = True                     # keep the component id
                else:
                    if is_hetero_line and resname in WATERS:
                        atomname = line[12:16].replace(" ", "")
                        element = element_of(atomname,
                                             line[76:78] if len(line) >= 78 else None)
                        if element not in ('H', 'D'):
                            # Held aside, not admitted: only the ones that turn
                            # out to coordinate a metal are kept. See WATERS.
                            waters.append((chain_id, resnum, icode, resname, atomname,
                                           float(line[30:38]), float(line[38:46]),
                                           float(line[46:54]), float(line[55:60]),
                                           entity, element, line[16:17].strip()))
                    continue

                if current_chain is None or current_chain.id != chain_id:  # new chain
                    if current_residue and len(current_residue.atoms) >= 1: # last residue of previous chain
                        current_chain.residues.append(current_residue) 
                    residues = []
                    current_chain = Chain(chain_id, residues)
                    current_protein.chains.append(current_chain)
                    current_residue = None

                if current_residue is None:  # first residue of the chain
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain,
                                              False, None, is_hetero, icode)
                    #current_chain.residues.append(current_residue)
                
                # A new residue is detected on resnum, insertion code OR
                # resname: a ligand and a protein residue can share a number
                # and only the name separates them, and 100 / 100A share both
                # number and name with only the code between them.
                if (current_residue.resnum != resnum
                        or current_residue.icode != icode
                        or current_residue.resname != resname):
                    if len(current_residue.atoms) >= 1:
                        current_chain.residues.append(current_residue)
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain,
                                              False, None, is_hetero, icode)

                if force_field:   # whether every record of it is ATOM (D32)
                    current_residue.force_field = (getattr(current_residue, "force_field", True)
                                                   and not is_hetero_line)
                                                                
                atomname = line[12:16].replace(" ", "")
                element = element_of(atomname, line[76:78] if len(line) >= 78 else None)
                if renames:
                    atomname = renames.get(atomname, atomname)   # MSE SE -> SD

                # OXT terminates a polypeptide, but it is also a legitimate atom
                # name inside some components (acetate, for one), so the rule is
                # applied to standard residues only.
                if atomname == "OXT" and not is_hetero:
                    continue # the residue is appended at the chain transition
                             # or at END, as in parse_cif
                elif element == 'H' or element == 'D':
                    continue
                
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                occupancy = float(line[55:60])
                altloc = line[16:17].strip()   # column 17
                
                atom = Atom(atomname, x, y, z, occupancy, current_residue, entity,
                            element, altloc) # creates atom

                # Alternate conformations: keep exactly ONE copy of each atom
                # name - the highest occupancy, and the first one on a tie.
                #
                # The previous test compared only against atoms[-1] and skipped
                # only when the new occupancy was STRICTLY lower, so a residue
                # modelled at 0.50/0.50 - very common in X-ray structures - kept
                # both conformers and every one of its contacts was reported
                # twice. Comparing against the whole residue also catches altloc
                # records that are not adjacent in the file.
                #
                # The winner is written back in place rather than appended, so
                # atom ordering is preserved; contacts.py relies on atoms[1]
                # being the alpha carbon.
                #
                # With all_altlocs (--all-altlocs) every alternate is kept
                # instead, one atom per (name, alternate letter), and the
                # detector pairs only atoms of one conformer. See ALTLOCS.
                for index, existing in enumerate(current_residue.atoms):
                    if existing.atomname == atomname and (
                            not all_altlocs or existing.altloc == atom.altloc):
                        if occupancy > existing.occupancy:
                            current_residue.atoms[index] = atom
                        break
                else:
                    current_residue.atoms.append(atom)
                # else:
                #     continue

                # CHECKING FOR AROMATICS
                if current_residue.resname in stacking:
                    allowed = stacking[current_residue.resname][1:]
                    all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                    
                    # if ring has only one conformation and the residue is complete (all atoms populated)
                    if all_atoms_have_occupancy_one and len({a.atomname for a in current_residue.atoms}) == stacking[current_residue.resname][0]:
                        ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms if atom.atomname in stacking[current_residue.resname]])
                        
                        if ring_atoms.any():
                            centroid_atom = centroid(current_residue, ring_atoms, entity)
                            current_residue.atoms.append(centroid_atom)
                            current_residue.ring = True # flags the aromatic residue

                            normal_vector = calc_normal_vector(ring_atoms)
                            current_residue.normal_vector = normal_vector
                            
            elif line.startswith("END"):
                if current_residue is not None and len(current_residue.atoms) >= 1:
                    current_chain.residues.append(current_residue) # appends the last residue  
                # Handling cases where there is no ID
                if current_protein.id is None:
                    id = str(pdb_file).split("/")[-1]
                    id = id.split(".")[0]
                    current_protein.id = id  
    
    native_force_field_residues(current_protein)
    nucleic = {chain: names for chain, names in seqres.items()
               if names & STANDARD_NUCLEOTIDES}
    if nucleic:
        drop_strand_components(
            current_protein,
            lambda r: r.resname in nucleic.get(r.chain.id, ()))
    admit_coordinating_waters(current_protein, waters)
    current_protein.links = links
    return label_alternates(current_protein), ph


def parse_cif(cif_file, all_altlocs=False):
    """
    Parses a .cif file to create a Protein object.

    Args:
        cif_file (str): Path to the cif file.

    Returns:
        Protein: The Protein object populated with chain, residue, and atom objects.

    This function reads a cif file, extracts information about chains, residues, and atoms,
    and constructs a Protein object. It handles specific fields for atom site data, processes
    only the first model for NMR files, and filters low-quality atoms. The function also assigns
    identifiers and titles based on the information from the file.
    """

    valid_atoms = ['N', 'C', 'O', 'S']
    
    current_protein = Protein()
    current_chain = None
    current_residue = None
    waters = []            # held aside; see WATERS
    nucleic_entities = set()  # entities that are nucleic-acid chains
    poly_columns = None       # _entity_poly loop header, while inside it
    poly_scalar_id = None     # _entity_poly written as scalars (one entity)
    conn_columns = None       # _struct_conn loop header, while inside it
    conn_values, conn_rows, conn_text = [], [], None
    conn_scalar = {}          # _struct_conn written as scalars (one link)
    atomsite_block = False # _atom_site. lines
    atominfo_block = False # ATOM        lines
    atom_lines = []
    models = []
    title = None
    title_block = False
    
    ph = 7.4
    experimental_lines = []
    nmr_expt = False

    with open(cif_file) as f:
        
        current_protein.id = os.path.basename(cif_file).split(".")[0]
        
        previous = ""
        for line in f:
            line = line.strip()
            previous, after_loop = line, previous == "loop_"

            if line.startswith("_entry.id"):
                values = cif_tokens(line)
                value = values[1] if len(values) > 1 else ""
                if value.upper() not in PLACEHOLDER_IDS:
                    current_protein.id = structure_id(value)

            # _entity_poly.type names each polymer entity's kind, which is how a
            # strand made only of modified nucleotides - 1KDH's four BRU - is
            # recognised as nucleic; see STANDARD_NUCLEOTIDES. A loop when there
            # are several polymer entities, scalars when there is one.
            if poly_columns is not None:
                if line.startswith("_entity_poly."):
                    poly_columns.append(line.split(".", 1)[1])
                    continue
                if line.startswith("#") or line.startswith("loop_") or line.startswith("_"):
                    poly_columns = None
                elif not line.startswith(";"):
                    values = cif_tokens(line)
                    if ("type" in poly_columns and len(values) > poly_columns.index("type")
                            and values[0].isdigit()):
                        if "nucleotide" in values[poly_columns.index("type")]:
                            nucleic_entities.add(values[0])
                    continue
            # _struct_conn: the links between residues (see IDENTITY_OPERATORS).
            # A loop, whose rows can wrap and can hold ;-delimited text fields,
            # or scalars when the file declares a single link.
            if conn_columns is not None:
                if line.startswith("_struct_conn."):
                    conn_columns.append(line.split(".", 1)[1])
                    continue
                if conn_text is not None:
                    if line.startswith(";"):
                        conn_values.append(" ".join(conn_text).strip() or "?")
                        conn_text = None
                    else:
                        conn_text.append(line)
                elif line.startswith(";"):
                    conn_text = [line[1:]]
                elif line.startswith("#") or line.startswith("loop_") or line.startswith("_"):
                    conn_columns = None
                else:
                    conn_values.extend(cif_tokens(line))
                if conn_columns is not None:
                    while len(conn_values) >= len(conn_columns):
                        conn_rows.append(dict(zip(conn_columns, conn_values)))
                        conn_values = conn_values[len(conn_columns):]
                    continue
            if line.startswith("_struct_conn."):
                values = cif_tokens(line)
                if after_loop:
                    conn_columns, conn_values = [values[0].split(".", 1)[1]], []
                elif len(values) > 1:
                    conn_scalar[values[0].split(".", 1)[1]] = values[1]
                continue

            if line.startswith("_entity_poly.entity_id"):
                values = cif_tokens(line)
                if len(values) > 1:
                    poly_scalar_id = values[1]
                else:
                    poly_columns = ["entity_id"]
                continue
            if line.startswith("_entity_poly.type") and poly_scalar_id is not None:
                values = cif_tokens(line)
                if len(values) > 1 and "nucleotide" in values[1]:
                    nucleic_entities.add(poly_scalar_id)
                continue

            # this title block can definitely be simplified, but there are a lot of edge cases to handle
            # current_protein.set_title() is set before returning the whole object
            if line.startswith("_struct.title"):
                title = line[len("_struct.title"):].strip()
                if title.startswith(";") and title.endswith(";"):
                    title = title[1:-1].strip()
                elif title.startswith("'") and title.endswith("'"):
                    title = title[1:-1].strip()
                elif title == "":
                    title_block = True
            elif title_block:
                if line.startswith(";"):
                    if title == "":
                        title = line[1:].strip()
                    else:
                        title_block = False
                else:
                    title += line.strip()
                    title_block = False
                    
            # Direct pH tags      
            if line.startswith("_exptl_crystal_grow.pH ") or line.startswith("_pdbx_nmr_exptl_sample_conditions.pH "):
                parts = line.split()
                if len(parts) > 1:
                    try:
                        ph = float(parts[1])
                    except ValueError:
                        continue
                    
            # Handling files where pH value is not directly after the line
            elif line.startswith("_pdbx_nmr_exptl_sample_conditions"):
                field = line.split(".")[1]
                experimental_lines.append(field)
                nmr_expt = True             
            elif nmr_expt and line.startswith("1"):
                parts = line.split()
                try:
                    nmr_ph_index = experimental_lines.index("pH")
                    ph = float(parts[nmr_ph_index])
                except (ValueError, IndexError):
                    pass
                nmr_expt = False
                
            if line.startswith("_atom_site.group_PDB"): # entering ATOM definition block
                atomsite_block = True
                line = line.split(".")[1]
                atom_lines.append(line)
                
            elif atomsite_block and line.startswith("_atom_site"):
                line = line.split(".")[1]
                atom_lines.append(line)
                
            elif atomsite_block and (line.startswith("ATOM") or line.startswith("HETATM")): # maps the order of the columns
                atomname_index = atom_lines.index("label_atom_id")
                resname_index = atom_lines.index("label_comp_id")
                chain_index = atom_lines.index("label_asym_id")
                chain_index2 = atom_lines.index("auth_asym_id")
                
                if "auth_seq_id" in atom_lines:
                    resnum_index = atom_lines.index("auth_seq_id")
                else:
                    resnum_index = atom_lines.index("label_seq_id")

                # The insertion code. '?' or '.' means there is none.
                icode_index = (atom_lines.index("pdbx_PDB_ins_code")
                               if "pdbx_PDB_ins_code" in atom_lines else None)

                # The alternate-location indicator, which is what separates two
                # chemical identities modelled at the same position.
                altloc_index = (atom_lines.index("label_alt_id")
                                if "label_alt_id" in atom_lines else None)
                
                x_index = atom_lines.index("Cartn_x")
                y_index = atom_lines.index("Cartn_y")
                z_index = atom_lines.index("Cartn_z")
                occupancy_index = atom_lines.index("occupancy")
                model_index = atom_lines.index("pdbx_PDB_model_num")
                atom_element_index = atom_lines.index("type_symbol")
                entity_index = atom_lines.index("label_entity_id")
                                                               
                atomsite_block = False
                atominfo_block = True
                
            if (line.startswith("ATOM") or line.startswith("HETATM")) and atominfo_block:
                is_hetero_line = line.startswith("HETATM")
                line = cif_tokens(line)

                # Hydrogens are excluded everywhere; the old N/C/O/S allowlist
                # is not applied to HETATM records, which is what previously
                # discarded every metal, phosphorus and halogen in a structure.
                element = element_of("", line[atom_element_index])
                if element in ('H', 'D'):
                    continue
                if not is_hetero_line and element not in valid_atoms:
                    continue
                
                models.append(int(line[model_index]))
                curr_model = int(line[model_index])
                if curr_model != models[0]: # parses only the first model (NMR files)
                    # Flush the residue in progress before leaving. Without this
                    # the last entity of model 1 is discarded, because nothing
                    # after it ever triggers a residue boundary. In 1t18 that
                    # silently dropped HC4, the chromophore - the only ligand in
                    # the structure.
                    if current_residue is not None and len(current_residue.atoms) >= 1:
                        current_chain.residues.append(current_residue)
                        current_residue = None
                    break
                    #return current_protein
                
                if line[chain_index] != ".":
                    chain_id = line[chain_index]
                else:
                    chain_id = line[chain_index2]
                
                resnum = int(line[resnum_index])
                icode = "" if icode_index is None else line[icode_index]
                if icode in (".", "?"):
                    icode = ""
                # if resnum <= 0:
                #     continue
                resname = line[resname_index]

                # force-field names: kept, decided per residue at the end (D32)
                force_field = resname in FORCE_FIELD_NAMES or resname in FORCE_FIELD_CAPS

                # modified residues that are chemically standard amino acids
                renames = STANDARD_ATOM_NAMES.get(resname)
                resname = HETERO_AS_STANDARD.get(resname, resname)

                if resname in residue_mapping:
                    resname = residue_mapping[resname]   # one-letter code
                    is_hetero = False
                elif force_field or (is_hetero_line and resname not in SKIP_HETERO):
                    is_hetero = True                     # keep the component id
                else:
                    if is_hetero_line and resname in WATERS:
                        # Held aside, not admitted: see WATERS.
                        alt = "" if altloc_index is None else line[altloc_index]
                        waters.append((chain_id, resnum, icode, resname,
                                       line[atomname_index],
                                       float(line[x_index]), float(line[y_index]),
                                       float(line[z_index]),
                                       float(line[occupancy_index]),
                                       line[entity_index] if line[entity_index] != "." else chain_id,
                                       element, alt))
                    elif not is_hetero_line and resname in STANDARD_NUCLEOTIDES:
                        # Skipped, but it marks its entity as a nucleic-acid
                        # chain; see STANDARD_NUCLEOTIDES.
                        nucleic_entities.add(line[entity_index]
                                             if line[entity_index] != "." else chain_id)
                    continue

                if current_chain is None or current_chain.id != chain_id:  # new chain
                    if current_residue and len(current_residue.atoms) >= 1: # last residue of previous chain
                        current_chain.residues.append(current_residue) 
                    residues = []
                    current_chain = Chain(chain_id, residues)
                    current_protein.chains.append(current_chain)
                    current_residue = None

                if current_residue is None:  # first residue of the chain
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain,
                                              False, None, is_hetero, icode)
                    #current_chain.residues.append(current_residue)
                
                # A new residue is detected on resnum, insertion code OR
                # resname: a ligand and a protein residue can share a number
                # and only the name separates them, and 100 / 100A share both
                # number and name with only the code between them.
                if (current_residue.resnum != resnum
                        or current_residue.icode != icode
                        or current_residue.resname != resname):
                    if len(current_residue.atoms) >= 1:
                        current_chain.residues.append(current_residue) 
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain,
                                              False, None, is_hetero, icode)

                if force_field:   # whether every record of it is ATOM (D32)
                    current_residue.force_field = (getattr(current_residue, "force_field", True)
                                                   and not is_hetero_line)
                                                                
                atomname = line[atomname_index]
                if renames:
                    atomname = renames.get(atomname, atomname)   # MSE SE -> SD
                # OXT is a legitimate atom name inside some components, so the
                # terminal-oxygen rule applies to standard residues only.
                if atomname == "OXT" and not is_hetero:
                    continue
                    
                x, y, z = float(line[x_index]), float(line[y_index]), float(line[z_index])
                occupancy = float(line[occupancy_index])
                
                if line[entity_index] == ".":
                    entity = chain_id
                else:
                    entity = line[entity_index]

                altloc = ("" if altloc_index is None else line[altloc_index])
                    
                atom = Atom(atomname, x, y, z, occupancy, current_residue, entity,
                            element, altloc) # creates atom

                # Alternate conformations: keep exactly ONE copy of each atom
                # name - the highest occupancy, and the first one on a tie.
                #
                # The previous test compared only against atoms[-1] and skipped
                # only when the new occupancy was STRICTLY lower, so a residue
                # modelled at 0.50/0.50 - very common in X-ray structures - kept
                # both conformers and every one of its contacts was reported
                # twice. Comparing against the whole residue also catches altloc
                # records that are not adjacent in the file.
                #
                # The winner is written back in place rather than appended, so
                # atom ordering is preserved; contacts.py relies on atoms[1]
                # being the alpha carbon.
                #
                # With all_altlocs (--all-altlocs) every alternate is kept
                # instead, one atom per (name, alternate letter), and the
                # detector pairs only atoms of one conformer. See ALTLOCS.
                for index, existing in enumerate(current_residue.atoms):
                    if existing.atomname == atomname and (
                            not all_altlocs or existing.altloc == atom.altloc):
                        if occupancy > existing.occupancy:
                            current_residue.atoms[index] = atom
                        break
                else:
                    current_residue.atoms.append(atom)
                # else:
                #     continue
                                
                # CHECKING FOR AROMATICS
                if current_residue.resname in stacking:
                    allowed = stacking[current_residue.resname][1:]
                    all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                    
                    # if ring has only one conformation and the residue is complete (all atoms populated)
                    if all_atoms_have_occupancy_one and len({a.atomname for a in current_residue.atoms}) == stacking[current_residue.resname][0]:
                        ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms if atom.atomname in stacking[current_residue.resname]])
                        if ring_atoms.any():
                            centroid_atom = centroid(current_residue, ring_atoms, entity)
                            current_residue.atoms.append(centroid_atom)
                            current_residue.ring = True # flags the aromatic residue

                            normal_vector = calc_normal_vector(ring_atoms)
                            current_residue.normal_vector = normal_vector

            elif atominfo_block and line == "#":
                if current_residue is not None and len(current_residue.atoms) >= 1:
                    current_chain.residues.append(current_residue) # appends the last residue
                atominfo_block = False
    
    if title is not None:
        current_protein.set_title(title.title().replace("'","").replace('"','').replace(",","."))
    else:
        current_protein.set_title(None)
        
    if ph == 0 or ph > 14:
        ph = 7.4

    native_force_field_residues(current_protein)
    if nucleic_entities:
        drop_strand_components(
            current_protein,
            lambda r: r.atoms[0].entity in nucleic_entities)
    admit_coordinating_waters(current_protein, waters)
    current_protein.links = struct_conn_links(conn_rows + ([conn_scalar] if conn_scalar else []))
    return label_alternates(current_protein), ph


def build_hetero_rings(residue, ring_table):
    """Centroid and plane normal for each aromatic ring system of a component.

    Ring membership comes from the generated table, which derives it from the
    Chemical Component Dictionary - structure files carry no connectivity, so
    ring perception is impossible at runtime. Each entry is one five- or
    six-membered ring; fused systems are deliberately NOT merged, so heme gives
    four pyrrole planes rather than one porphyrin (F5 in docs/DECISIONS.md).

    Returns a list of (centroid Atom, normal vector). Incomplete rings are
    skipped: a partially resolved ligand has no meaningful plane.
    """

    systems = ring_table.get(residue.resname)
    if not systems:
        return []

    by_name = {atom.atomname: atom for atom in residue.atoms}
    built = []

    for index, atom_names in enumerate(systems, start=1):
        present = [by_name[n] for n in atom_names if n in by_name]
        if len(present) < len(atom_names):
            continue                      # incomplete ring, no reliable plane
        if any(atom.occupancy != 1 for atom in present):
            continue                      # alternate conformations, as for residues

        coords = array([[a.x, a.y, a.z] for a in present])
        if not coords.any():
            continue

        label = "RNG" if len(systems) == 1 else f"RNG{index}"
        centre = mean(coords, axis=0)
        centre_atom = Atom(label, centre[0], centre[1], centre[2], 1, residue,
                           present[0].entity, "R")
        built.append((centre_atom, calc_normal_vector(coords)))

    return built


def centroid(residue, ring_atoms, entity):
    """
    Calculates the centroid of a set of ring atoms and creates a centroid Atom.

    Args:
        residue (Residue): The residue containing the ring atoms.
        ring_atoms (array): Array of ring atom coordinates.

    Returns:
        Atom: A new Atom object representing the centroid of the ring.

    The function computes the centroid of the given ring atoms and creates a new Atom object
    with this centroid position. The atom is labeled as "RNG" for ring centroid.
    """
    
    centroid = mean(ring_atoms, axis = 0)
    centroid_atom = Atom("RNG", centroid[0], centroid[1], centroid[2], 1, residue, entity)
    
    return centroid_atom


def admit_coordinating_waters(protein, waters):
    """Admit the waters that coordinate a metal, and only those.

    A water within the metal-coordination window of a metal atom is a ligand of
    that metal centre, not incidental solvent: the sixth coordination position
    of a resting-state cytochrome P450 is an aqua ligand, 2.07 A from the heme
    iron in PDB 1TQN, and excluding water dropped that bond entirely.

    `waters` holds the records set aside during parsing, as plain tuples -
    nothing is materialised into a Residue until it earns admission, so a
    structure with a thousand waters and no metal pays for a thousand tuples
    and nothing else. With no metals present the function returns immediately.

    An admitted water is a hetero residue with a single oxygen and an all-zero
    property vector, exactly as a metal is, so the only contact it can form is
    the coordination itself. It is appended to its own chain, which keeps it
    out of the way of the polymer.
    """
    if not waters:
        return protein

    from src.conditions import METALS, categories

    metals = [(atom.x, atom.y, atom.z)
              for residue in protein.get_residues() if residue.is_hetero
              for atom in residue.atoms if atom.element in METALS]
    if not metals:
        return protein

    window = categories['metal_coordination'][1]
    chains = {chain.id: chain for chain in protein.chains}

    # One record per position, highest occupancy, first on a tie - the same
    # rule the atom path applies, and water needs it too: a water modelled in
    # two alternate positions arrives here as two records with one residue
    # number, and both would coordinate the metal, reporting one bond twice.
    # 17 structures in a 245,839-entry validation did exactly that.
    best = {}
    for record in waters:
        chain_id, resnum, icode = record[0], record[1], record[2]
        occupancy = record[8]
        key = (chain_id, resnum, icode)
        if key not in best or occupancy > best[key][8]:
            best[key] = record

    for (chain_id, resnum, icode, resname, atomname,
         x, y, z, occupancy, entity, element, altloc) in best.values():
        if not any((x - mx) ** 2 + (y - my) ** 2 + (z - mz) ** 2 <= window ** 2
                   for mx, my, mz in metals):
            continue

        chain = chains.get(chain_id)
        if chain is None:
            chain = Chain(chain_id, [])
            protein.chains.append(chain)
            chains[chain_id] = chain

        residue = Residue(resnum, resname, [], chain, False, None, True, icode)
        residue.atoms.append(Atom(atomname, x, y, z, occupancy, residue,
                                  entity, element, altloc))
        chain.residues.append(residue)

    return protein


ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def label_alternates(protein):
    """Label apart the residues modelled at one and the same position.

    Microheterogeneity: a position carrying two or three different chemical
    identities, each at partial occupancy, distinguished in the file by the
    alternate-location column. 9C21 models position A/73 as the unmodified
    His73 (altloc B, occupancy 0.50) and the 2-methylhistidine HIC73 (altloc A,
    0.50); 3NED models A/67 as three chromophore variants, NRQ, CH6 and EYG at
    altlocs A, B and C. Roughly 0.14% of the PDB.

    Every alternative is kept - both are deposited data, and which one matters
    depends on the question being asked - but without a label their records are
    indistinguishable in the output beyond the residue name. Each gets the
    file's own alternate-location letter appended to its number: `73A`, `73B`.
    A residue whose alternatives carry no letter falls back to a sequential
    one, and a suffix that would collide with another residue's label in the
    same chain is numbered instead.

    Called once at the end of parsing, because crowding is a property of the
    position rather than of the residue. DISPLAY ONLY: pruning, the contact
    registry and the validation guard key on (chain, resnum, icode, resname),
    which already separates these residues, so a label can never change which
    contacts are found.
    """
    for chain in protein.chains:
        positions = {}
        for residue in chain.residues:
            positions.setdefault((residue.resnum, residue.icode), []).append(residue)

        taken = {f"{r.resnum}{r.icode}" for r in chain.residues}
        for (resnum, icode), group in positions.items():
            if len({r.resname for r in group}) < 2:
                continue          # one identity, nothing to separate

            used = set()
            for index, residue in enumerate(group):
                suffix = residue.altloc
                if not suffix or suffix in used:
                    # No letter in the file, or the same letter used twice at
                    # this position: fall back to the first free one.
                    suffix = next((c for c in ALPHABET if c not in used),
                                  str(index + 1))
                while f"{resnum}{icode}{suffix}" in taken:
                    suffix += "1"        # a real insertion code got there first
                used.add(suffix)
                taken.add(f"{resnum}{icode}{suffix}")
                residue.alt_suffix = suffix

    return protein


def calc_normal_vector(ring_atoms):
    """
    Computes the normal vector to the plane of a set of ring atoms.

    Args:
        ring_atoms (array): Array of ring atom coordinates.

    Returns:
        array: The normal vector to the plane of the ring atoms.

    The function calculates the normal vector to the plane defined by the given ring atoms
    using Singular Value Decomposition (SVD). The normal vector is extracted from the last
    row of the V^T matrix obtained from SVD.
    """
    
    centroid = mean(ring_atoms, axis = 0) # axis=0 -> mean through the columns
    centered_ring_atoms = ring_atoms - centroid # normalizes to origin
    
    # Use SVD to calculate the plane
    _, _, vh = svd(centered_ring_atoms) # vh = V^T
    normal_vector = vh[2]  # normal vector is the last row of the V^T matrix
    
    return normal_vector
