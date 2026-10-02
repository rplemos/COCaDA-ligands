"""
COCaDA v1.6
Author: Rafael Lemos - rafaellemos@ufmg.br
Date: 26/03/2026

License: MIT License
"""


class ProcessingContext:
    """
    Stores processing parameters passed by command-line.

    Attributes:
        core (str or None): The core selection for processing. Defaults to None.
        output (str or None): The output file path or identifier. Defaults to None.
        region (list or None): Residue numbers to restrict detection to. Defaults to None.
        chains (list or None): Chain ids to restrict detection to. Defaults to None,
            meaning every chain. It used to default to False, which
            contact_detection iterated and crashed on.
        interface (str or None): Path to an interface residue list. Defaults to None.
    """

    def __init__(self, distances, core=None, output=None, region=None, chains=None, interface=None, modified_distances=False, epsilon=0, ph=7.4, silent=False, uncertainty_flags=None, interchain=False, web=False, all_altlocs=False, min_plddt=None):
        self.core = core
        self.output = output
        self.region = region
        self.chains = chains
        self.interface = interface
        self.distances = distances
        self.modified_distances = modified_distances
        self.epsilon = epsilon
        self.ph = ph
        self.silent = silent
        self.uncertainty_flags = uncertainty_flags
        self.interchain = interchain
        self.web = web
        self.all_altlocs = all_altlocs
        self.min_plddt = min_plddt

class Protein:
    """
    Represents a protein structure, including its title, ID, and chains.

    Attributes:
        title (str or None): The title of the protein.
        id (str or None): The ID of the protein.
        chains (list of Chain): A list of Chain objects representing the chains in the protein.
    """
    
    
    def __init__(self):
        """
        Initializes a new Protein instance with default values.
        """
        
        self.title = None
        self.id = None
        self.chains = []

        # Covalent links the file declares between residues (_struct_conn,
        # LINK), as ((key, atom name), (key, atom name)); see parser.py.
        self.links = []

    def set_title(self, title):
        """
        Sets the title of the protein. If a title is already set, appends the new title.

        Args:
            title (str): The title to set or append to the protein's title.
        """
        
        if self.title is None:
            self.title = title
        else:
            self.title += " " + title.strip()

    def get_chains(self):
        """
        Generator that yields each chain in the protein.

        Yields:
            Chain: A chain object in the protein.
        """
        
        for chain in self.chains:
            yield chain

    def get_residues(self):
        """
        Generator that yields each residue in all chains of the protein.

        Yields:
            Residue: A residue object in the protein's chains.
        """
        
        for chain in self.get_chains():
            for residue in chain.residues:
                yield residue
                
    def true_count(self):
        """
        Counts the total number of residues in the protein.

        Returns:
            int: The total number of residues in the protein.
        """
        
        return sum(1 for _ in self.get_residues())
            

class Chain:
    """
    Represents a chain in a protein, containing residues.

    Attributes:
        id (str): The ID of the chain.
        residues (list of Residue): A list of Residue objects representing the residues in the chain.
    """

    def __init__(self, id, residues):
        """
        Initializes a new Chain instance.
        """
        
        self.id = id
        self.residues = residues


class Residue:
    """
    Represents a residue in a protein chain.

    Attributes:
        resnum (int): The residue number.
        icode (str): The insertion code, '' when there is none. Part of the
            residue's identity: 100, 100A and 100B are three different
            residues sharing one number, which is how antibody CDRs are
            numbered.
        resname (str): The residue name (e.g., 'ALA' for alanine).
        atoms (list of Atom): A list of Atom objects in the residue.
        chain (Chain): The chain to which the residue belongs.
        ring (bool): Indicates whether the residue has a ring structure.
        normal_vector (tuple): The normal vector associated with the residue.
    """
    
    def __init__(self, resnum, resname, atoms, chain, ring, normal_vector,
                 is_hetero=False, icode=""):
        """
        Initializes a new Residue instance.
        """
        
        self.resnum = resnum

        # The insertion code, kept separate so resnum stays an integer for the
        # sequence-distance filter and the -reg region selection, while the
        # residue's identity and its printed label carry the code.
        self.icode = icode or ""
        self.resname = resname
        self.atoms = atoms
        self.chain = chain
        self.ring = ring
        self.normal_vector = normal_vector

        # Non-amino-acid entities: prosthetic groups, cofactors, ligands, ions.
        # These carry a 3-to-5 character component id as resname (e.g. 'HEM')
        # rather than the one-letter code used for standard residues.
        self.is_hetero = is_hetero

        # Aromatic ring systems available for stacking, as a list of
        # (centroid Atom, normal vector) pairs. A standard aromatic residue has
        # exactly one, which keeps its behaviour identical to before. Hetero
        # components can have several, one per five- or six-membered ring - NAD
        # has three, and heme four pyrroles, because the centroid of a fused
        # system sits too far from its own edges for the stacking window (F5
        # in docs/DECISIONS.md). Populated by contacts.set_anchors().
        self.rings = []

        # Pruning anchor. Standard residues use their alpha carbon and the
        # mined Ca-Ca distance table; hetero entities have no Ca, so they use
        # their centroid together with a bounding radius computed from the
        # coordinates in the file. Both are set by contacts.set_anchors().
        self.anchor = None
        self.radius = 0.0

        # Suffix separating this residue from the others modelled at the same
        # position. Empty for all but the ~0.14% of structures that carry
        # microheterogeneity, and set by parser.label_alternates() once the
        # whole structure is known, because crowding is a property of the
        # position rather than of the residue. Display only: the pruning, the
        # contact registry and the validation guard all key on
        # (chain, resnum, icode, resname), which already separates these.
        self.alt_suffix = ""

    @property
    def altloc(self):
        """The residue's alternate-location letter, '' when it has none."""
        for atom in self.atoms:
            if atom.altloc:
                return atom.altloc
        return ""

    @property
    def label(self):
        """The residue number as it should appear in output.

        `100` normally, `100A` when the file gives an insertion code, and
        `73B` for one of two chemical identities modelled at position 73.
        """
        return f"{self.resnum}{self.icode}{self.alt_suffix}"

    @property
    def key(self):
        """Identity within a structure: chain, number, insertion code, name."""
        return (self.chain.id, self.resnum, self.icode, self.resname)
  
      
class Atom:
    """
    Represents an atom in a residue.

    Attributes:
        atomname (str): The name of the atom (e.g., 'CA' for alpha carbon).
        x (float): The x-coordinate of the atom.
        y (float): The y-coordinate of the atom.
        z (float): The z-coordinate of the atom.
        occupancy (float): The occupancy value of the atom.
        residue (Residue): The residue to which the atom belongs.
        entity (int): The entity in which the atom is located.
        altloc (str): The alternate-location indicator, '' when there is none.
    """
    
    def __init__(self, atomname, x, y, z, occupancy, residue, entity, element=None,
                 altloc="", bfactor=None):
        """
        Initializes a new Atom instance.
        """
        
        self.atomname = atomname
        # Element symbol, taken from the file's own element column where the
        # format provides one (PDB columns 77-78, mmCIF `type_symbol`).
        # Required for metal coordination, which is defined by element rather
        # than by the property vector, and the only reliable way to tell a
        # calcium ion ('CA  ') from an alpha carbon (' CA ') in PDB format.
        self.element = (element or (atomname[:1] if atomname else '')).strip().upper()
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.residue = residue
        self.entity = entity

        # Alternate-location indicator: PDB column 17, mmCIF `label_alt_id`.
        # It is what separates two chemical identities modelled at one position
        # - 9C21's His73 (altloc B) from the 2-methylhistidine HIC73 (altloc A)
        # - and it is the letter used to label them apart in the output.
        self.altloc = "" if altloc in (None, ".", "?") else altloc.strip()

        # B-factor column, read only for --plddt, where a model keeps its pLDDT
        # there. None otherwise.
        self.bfactor = bfactor


#: Single-character CCD component identifiers are indistinguishable from the
#: one-letter codes COCaDA writes for standard residues: `K` is both potassium
#: and lysine, `F` both fluoride and phenylalanine, `W` both tungsten and
#: tryptophan, and `A`, `C`, `G`, `U` are both ribonucleotides and amino acids.
#: Nineteen such components exist in the Chemical Component Dictionary. In the
#: output they are written in LOWER CASE, so `k` is the ion and `K` the lysine.
#: Longer identifiers cannot collide and are left exactly as deposited: HEM is
#: HEM. The transformation is applied at output only - the component's real
#: identifier is what the property table is keyed on and what every internal
#: lookup uses.
def display_resname(resname, is_hetero):
    """The residue name as it should appear in output."""
    return resname.lower() if is_hetero and len(resname) == 1 else resname


def csv_field(value):
    """One CSV field, quoted only when it has to be (RFC 4180).

    Needed because a deposited atom name may contain a comma. 105 atom
    identifiers in the Chemical Component Dictionary do - `"C1,"` and `"O3,"`
    in nine sugar-bearing components, where a comma stands in for the usual
    prime - and 31 contain a double quote. Written raw, such a name adds a
    column to its row and silently shifts everything after it.
    """
    text = str(value)
    if any(ch in text for ch in (',', '"', '\n', '\r')):
        return '"' + text.replace('"', '""') + '"'
    return text


class Contact: 
    """
    Represents a contact between two atoms in different residues.

    Attributes:
        id1 (str): The ID of the first residue.
        chain1 (str): The chain of the first residue.
        residue_num1 (int): The residue number of the first residue.
        residue_name1 (str): The residue name of the first residue.
        atom1 (str): The atom name of the first atom.
        id2 (str): The ID of the second residue.
        chain2 (str): The chain of the second residue.
        residue_num2 (int): The residue number of the second residue.
        residue_name2 (str): The residue name of the second residue.
        atom2 (str): The atom name of the second atom.
        distance (float): The distance between the two atoms.
        type (str): The type of contact (e.g., hydrogen bond, hydrophobic).
        atom_object1 (Atom): The Atom object representing the first atom.
        atom_object2 (Atom): The Atom object representing the second atom.
    """
    
    def __init__(self, id1, chain1, residue_num1, residue_name1, atom1, 
                 id2, chain2, residue_num2, residue_name2, atom2, 
                 distance, type, atom_object1, atom_object2, is_uncertain=False,
                 icode1="", icode2="", alt1="", alt2="", altloc1="", altloc2=""):
        """
        Initializes a new Contact instance.
        """
        
        self.id1 = id1
        self.chain1 = chain1
        self.residue_num1 = residue_num1
        self.residue_name1 = residue_name1
        self.atom1 = atom1
        self.id2 = id2
        self.chain2 = chain2
        self.residue_num2 = residue_num2
        self.residue_name2 = residue_name2
        self.atom2 = atom2
        self.distance = distance
        self.type = type
        self.atom_object1 = atom_object1
        self.atom_object2 = atom_object2
        self.is_uncertain = is_uncertain

        # Insertion codes, '' when there are none. Kept beside the numbers
        # rather than folded into them so residue_num stays an integer.
        self.icode1 = icode1 or ""
        self.icode2 = icode2 or ""

        # Suffix separating one of two chemical identities modelled at the same
        # position. Display only - nothing keys on it - so a validation guard
        # comparing (chain, resnum, icode, resname, atom) is unaffected by a
        # position becoming crowded, or ceasing to be, between two runs.
        self.alt1 = alt1 or ""
        self.alt2 = alt2 or ""

        # Alternate-location letter of each ATOM, '' when it has none. Only
        # meaningful, and only printed, with --all-altlocs, where every
        # conformer is reported; by default one is kept per atom.
        self.altloc1 = altloc1 or ""
        self.altloc2 = altloc2 or ""

    @property
    def label1(self):
        """First residue's number as printed: 100, 100A, or 73B."""
        return f"{self.residue_num1}{self.icode1}{self.alt1}"

    @property
    def label2(self):
        return f"{self.residue_num2}{self.icode2}{self.alt2}"

    @staticmethod
    def _is_hetero(atom):
        """Whether an atom belongs to a non-amino-acid component.

        Read from the parsed structure rather than guessed from the residue
        name, which cannot be trusted: a one-character component identifier
        looks exactly like an amino-acid code.
        """
        residue = getattr(atom, "residue", None)
        return bool(getattr(residue, "is_hetero", False))

    @property
    def hetero1(self):
        return self._is_hetero(self.atom_object1)

    @property
    def hetero2(self):
        return self._is_hetero(self.atom_object2)

    @property
    def is_amino_pair(self):
        """True when both partners are standard amino acids.

        This is the `AA-AA` column of the output. It is the reliable way to
        tell a contact of the protein with itself from one involving a
        prosthetic group, cofactor, ion or ligand - reliable in a way reading
        the residue names is not.
        """
        return not (self.hetero1 or self.hetero2)

    def print_text(self, altlocs=False):
        """
        Generates a formatted string describing the contact between two atoms.

        Returns:
            str: A string describing the contact, including chain, residue, atom information, and distance.
        """
        
        map_type = {
            "hydrogen_bond":"HB",
            "hydrophobic":"HY",
            "attractive":"AT",
            "repulsive":"RE",
            "salt_bridge":"SB",
            "disulfide_bond":"DS",
            "stacking-other":"AS",
            "stacking-parallel":"AS", # on v.1 all aromatic stackings will be considered the same
            "stacking-perpendicular":"AS", # need to reimplement later
            "uncertain_attractive": "uAT",
            "uncertain_repulsive": "uRE",
            "uncertain_salt_bridge": "uSB",
            "metal_coordination": "MC"
        }
        
        # Named attributes, deliberately. This used to read self.__dict__ by
        # position - all_values[11] was the type - so adding any field to
        # Contact silently shifted every column of the CSV.
        return ",".join(csv_field(field) for field in (
            self.chain1, self.label1,
            display_resname(self.residue_name1, self.hetero1), self.atom1,
            self.chain2, self.label2,
            display_resname(self.residue_name2, self.hetero2), self.atom2,
            self.distance, map_type[self.type],
            "y" if self.is_amino_pair else "n",
        ) + ((self.altloc1, self.altloc2) if altlocs else ()))

        # verbose
        # return (f"{self.chain1}-{self.residue_num1}{self.residue_name1}:{self.atom1} and "
        #         f"{self.chain2}-{self.residue_num2}{self.residue_name2}:{self.atom2}: "
        #         f"{self.distance} A. {self.type.capitalize()}")

