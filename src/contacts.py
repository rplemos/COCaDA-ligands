"""
COCaDA-ligands (v2.0)
Author: Rafael Lemos - rafaellemos@ufmg.br
Date: 29/10/2026

License: MIT License
"""


from math import dist
from numpy import dot, arccos, degrees
from numpy.linalg import norm
from copy import deepcopy
from collections import defaultdict

from src.classes import Atom, Contact
from src.distances import distances
import src.conditions as conditions
from src.parser import build_hetero_rings, WATERS
import src.ligand_data as ligand_data
import sys


def set_anchors(residues, ceiling, ring_table=None):
    """Give every residue a pruning anchor and, for hetero entities, a radius.

    Standard residues keep the existing behaviour exactly: the anchor is
    atoms[1], the alpha carbon, and pruning uses the mined Ca-Ca distance
    table. Nothing about amino-acid pairs changes.

    Hetero entities have no alpha carbon and no mined table - and none can be
    mined for a component that appears in three structures - so they are given
    a centroid anchor plus the maximum distance from that centroid to any of
    their atoms. Two entities can then be skipped when

        dist(anchor1, anchor2) > radius1 + radius2 + ceiling

    which is an upper bound rather than an observed maximum: unlike the mined
    table it cannot discard a pair that really does have atoms in contact.
    It costs one pass over the atoms of each entity.

    Known limitation - the bound degenerates for elongated entities. A lipid
    or retinoid has a radius large enough that radius1 + radius2 + ceiling
    exceeds the Arg-Arg maximum, at which point it prunes nothing. Such
    entities need splitting or a cell list; see the design record.
    """

    for residue in residues:
        if not residue.is_hetero:
            residue.anchor = residue.atoms[1] if len(residue.atoms) > 1 else None
            # Reach of the side chain from the alpha carbon. Unused on the
            # amino-acid path, which prunes with the mined Ca-Ca table, but
            # required whenever this residue is paired with a hetero entity:
            # leaving it at zero ignores the side chain entirely, so an
            # arginine whose NH1 reaches 7 A past its Ca would be pruned away
            # while genuinely in contact.
            if residue.anchor is not None:
                a = residue.anchor
                residue.radius = max(dist((a.x, a.y, a.z), (o.x, o.y, o.z))
                                     for o in residue.atoms)
            else:
                residue.radius = 0.0
            # Standard aromatic residues already had their centroid appended as
            # the final RNG atom during parsing; expose it through the same
            # list the hetero path uses, so stacking has one code path. This is
            # exactly the ring that was used before, so results do not move.
            if residue.ring and residue.normal_vector is not None:
                residue.rings = [(residue.atoms[-1], residue.normal_vector)]
            continue

        atoms = residue.atoms
        if not atoms:
            residue.anchor = None
            continue

        count = len(atoms)
        cx = sum(a.x for a in atoms) / count
        cy = sum(a.y for a in atoms) / count
        cz = sum(a.z for a in atoms) / count
        residue.anchor = Atom("CEN", cx, cy, cz, 1, residue, atoms[0].entity)
        residue.radius = max(dist((cx, cy, cz), (a.x, a.y, a.z)) for a in atoms)

        # Ring perception is impossible from coordinates alone, so membership
        # comes from the generated Chemical Component Dictionary table.
        residue.rings = build_hetero_rings(residue, ring_table or {})
        if residue.rings:
            residue.ring = True


#: Closer than this, two atoms cannot both be present: no non-covalent contact
#: is this short, and a covalent bond between two separate components would be
#: declared by the deposit, not inferred.
OVERLAP = 1.5

#: Occupancies of two alternatives sum to at most one; the margin absorbs the
#: rounding of refined occupancies (0.49 + 0.52).
OCCUPANCY_SUM = 1.05


def mutually_exclusive(residues):
    """Pairs of hetero residues that are alternatives of one another (D26).

    A ligand on a crystallographic two-fold axis is deposited as two
    half-occupancy orientations, and a biological assembly places both in the
    same site; a ligand and a partial-occupancy buffer molecule can share a
    site the same way. The parser's alternate-conformation rule works inside a
    residue, so nothing knew such residues exclude each other, and the contacts
    between them - hundreds per site, hydrogen bonds at near-zero distance
    among them - were 39.8% of all ligand contacts in the TTR assemblies.

    Two hetero residues are alternatives when an atom of one lies within
    OVERLAP of an atom of the other and their occupancies sum to no more than
    one. No contact is reported between them; each keeps its contacts with
    everything else. Amino-acid residues are never considered, so the
    amino-acid path is untouched. Returns a set of (id(a), id(b)) in the order
    the residues are listed.
    """
    partial = []
    for residue in residues:
        if residue.is_hetero and residue.anchor is not None and residue.atoms:
            occupancy = max(atom.occupancy for atom in residue.atoms)
            if occupancy < OCCUPANCY_SUM:
                partial.append((residue, occupancy))
    pairs = set()
    for i, (a, occ_a) in enumerate(partial):
        for b, occ_b in partial[i + 1:]:
            if occ_a + occ_b > OCCUPANCY_SUM:
                continue
            if dist((a.anchor.x, a.anchor.y, a.anchor.z),
                    (b.anchor.x, b.anchor.y, b.anchor.z)) > a.radius + b.radius + OVERLAP:
                continue
            if any(dist((p.x, p.y, p.z), (q.x, q.y, q.z)) < OVERLAP
                   for p in a.atoms for q in b.atoms):
                pairs.add((id(a), id(b)))
    return pairs


#: Two atoms of one residue this close are taken as bonded. It covers the bonds
#: between C, N, O, S and P; longer ones - S-S, C-Br, C-I, anything to Se - are
#: missed, which leaves a contact reported rather than hiding a real one.
COVALENT_BOND = 1.9


def bonded_neighbourhood(residue, atomname):
    """An atom of a residue and every atom of that residue bonded to it.

    Bonds inside a residue are not in the structure file; at this distance they
    are unambiguous. Under --all-altlocs every copy of the atom counts.
    """
    centres = [a for a in residue.atoms if a.atomname == atomname]
    return {a.atomname for a in residue.atoms if not a.atomname.startswith("RNG")
            for c in centres
            if dist((a.x, a.y, a.z), (c.x, c.y, c.z)) <= COVALENT_BOND}


def covalent_exclusions(residues, links):
    """Atom pairs that sit around a declared covalent bond (D30).

    For a bond A-B between two residues, the pairs formed by A and the atoms
    bonded to it, against B and the atoms bonded to it, are one to three bonds
    apart - 1.3 to 2.5 A, set by the bond, not by any interaction - and are not
    reported. Everything else between the two residues still is: in 2GUY the
    Asn197 ND2 - NAG O5 "hydrogen bond" at 2.39 A goes, ND2 - NAG O7 at 3.4 A
    stays. Ring stacking is untouched.

    Only links with a component on at least one side. A covalent link between
    two standard residues - an isopeptide, a crosslink, a cyclic peptide's
    closing bond - would move amino-acid results and is left for the author.
    Links whose atom is a metal are coordination (PDB writes them on LINK
    records too), reported as MC on purpose.

    Returns {(id(residue1), id(residue2)): {(atom name 1, atom name 2), ...}},
    in both orders.
    """
    if not links:
        return {}
    from src.conditions import METALS
    by_key = {}
    for residue in residues:
        by_key[residue.key] = residue
        source = getattr(residue, "source_name", None)   # a force-field name (D32)
        if source:
            by_key[(residue.chain.id, residue.resnum, residue.icode, source)] = residue
    pairs = {}
    for (key1, name1), (key2, name2) in links:
        first, second = by_key.get(key1), by_key.get(key2)
        if first is None or second is None or first is second:
            continue
        if not (first.is_hetero or second.is_hetero):
            continue
        if any(a.element in METALS for r, n in ((first, name1), (second, name2))
               for a in r.atoms if a.atomname == n):
            continue
        around1 = bonded_neighbourhood(first, name1)
        around2 = bonded_neighbourhood(second, name2)
        pairs.setdefault((id(first), id(second)), set()).update(
            (x, y) for x in around1 for y in around2)
        pairs.setdefault((id(second), id(first)), set()).update(
            (y, x) for x in around1 for y in around2)
    return pairs


#: What a contact predicate reads, by contact type: the general vector (and,
#: for hydrogen bonds, the helix and D29 skips), the electrostatic vector, the
#: table names, or the atoms themselves.
GENERAL, HBOND, ELECTROSTATIC, DISULFIDE, METAL = range(5)

#: The vector of an atom with no flag: a water, or a metal for the
#: electrostatic types. Read only.
NO_FLAGS = (0, 0, 0, 0, 0, 0)


def _kind(contact_type):
    if contact_type == 'hydrogen_bond':
        return HBOND
    if contact_type in ('salt_bridge', 'attractive', 'repulsive'):
        return ELECTROSTATIC
    if contact_type == 'disulfide_bond':
        return DISULFIDE
    if contact_type == 'metal_coordination':
        return METAL
    return GENERAL


def contact_detection(protein, context, uncertainty_flags, local_contact_types):
    """
    Detects contacts between atoms in a given protein.

    Args:
        protein (Protein): The protein object containing chain, residue and atom objects.

    Returns:
        list: A list of Contact objects representing the detected contacts.
    """
        
    region, chains, interface, interchain, argument_distances, modified_distances, epsilon = (
        context.region, context.chains, context.interface, context.interchain, context.distances, context.modified_distances, context.epsilon
    )
    # --all-altlocs: every alternate conformation was kept by the parser, and
    # atoms are paired only within one conformer (see parser.ALTLOCS).
    all_altlocs = getattr(context, "all_altlocs", False)

    prot_chains = []
    for chain in protein.chains:
        prot_chains.append(chain.id)

    # Validate the requested chains ONCE. This check used to sit inside the
    # residue pair loop, where it cost O(pairs x chains): its answer cannot
    # change from one pair to the next, and with -c naming every chain of a
    # 2,121-residue structure it was 6 of the 7.8 seconds of detection - more
    # than the atom loop itself. Admitting HETATM made that worse rather than
    # better, because mmCIF gives each component its own label_asym_id: 3FZN
    # has 4 polymer chains and 33 once its ligands are admitted.
    if chains is not None:
        for chain in chains:
            if chain not in prot_chains:
                sys.exit(f"Error: Chain {chain} is not present in the protein.")

    residues = list(protein.get_residues())
    contacts = []
    interface_res = []
    max_ca_distance = 20.47 # 0.01 higher than the Arg-Arg pair

    contacts_by_pair = defaultdict(list) # empty list under the key if it doesn't exist
    contact_registry = defaultdict(set) # empty set under the key if it doesn't exist

    categories = argument_distances if modified_distances else conditions.categories

    # Global ceiling on any contact distance: the atom-loop gate, the hetero
    # bound, and the reach the Ca table has to be widened by.
    ceiling = max(hi for _, hi in categories.values())
    disulfide_range = categories['disulfide_bond']

    # The Ca table was mined for contacts up to 6 A. A wider window needs every
    # entry widened by the excess, epsilon - once, as the paper describes. It
    # used to be added twice: here, and again in the pruning test below.
    # Derived from the windows themselves, so a caller passing custom distances
    # without an epsilon still gets a table that matches them.
    epsilon = max(epsilon, ceiling - 6)
    if epsilon > 0:
        max_ca_distance += epsilon
        updated_distances = {key: value + epsilon for key, value in distances.items()}
    else:
        updated_distances = distances

    # Load property vectors and ring definitions for the non-amino-acid
    # components actually present, and nothing else. Typically one to ten
    # entries out of the ~51 000 in the Chemical Component Dictionary.
    hetero_ids = {r.resname for r in residues if r.is_hetero}
    hetero_types, hetero_rings = ligand_data.load_for(hetero_ids)
    if hetero_types:
        local_contact_types.update(hetero_types)

    set_anchors(residues, ceiling, hetero_rings)
    exclusive = mutually_exclusive(residues)
    covalent = covalent_exclusions(residues, getattr(protein, "links", None))

    # Per-atom answers, worked out once per structure instead of for every
    # atom pair in the loop below: the table name, whether the atom takes part
    # at all, and the vectors the contact predicates read. The electrostatic
    # types (salt bridge, attractive, repulsive) read their own vector: the
    # pH-resolved one for an uncertain titratable atom, and all zeros for a
    # metal - its charge is expressed as coordination, not as an ion pair, so
    # a zinc 2.5 A from an aspartate does not report both. A water is admitted
    # only to coordinate a metal, so its vectors are all zeros; an atom the
    # table lacks takes part only if it is a metal, by element.
    for residue in residues:
        typed = []
        for atom in residue.atoms:
            key = f"{atom.residue.resname}:{atom.atomname}"
            vector = local_contact_types.get(key)
            metal = atom.element in conditions.METALS
            water = atom.residue.resname in WATERS
            atom.key = key
            atom.uncertain = key in uncertainty_flags
            if vector is None and not metal:
                continue                    # untyped: the RNG pseudo-atom, a deposit's own name
            atom.props = NO_FLAGS if water or vector is None else vector
            if atom.uncertain:
                atom.props_elec = resolve_uncertainty(key, uncertainty_flags, local_contact_types)
            elif metal or water or vector is None:
                atom.props_elec = NO_FLAGS
            else:
                atom.props_elec = vector
            typed.append(atom)
        residue.typed_atoms = typed

    # The contact types with what each predicate reads, in the order of
    # `categories` - the registry below keeps the first salt bridge, so the
    # order is part of the result. Two amino acids cannot coordinate a metal,
    # so their pairs skip that type.
    types_all = [(contact_type, low, high, _kind(contact_type),
                  conditions.contact_conditions[contact_type])
                 for contact_type, (low, high) in categories.items()]
    types_amino = [entry for entry in types_all if entry[3] != METAL]

    if interface:
        with open(interface,"r") as f:
            for line in f:
                interface_res.append(line.strip())
        
    # D1: the outer loop used to iterate over residues[1:] while the inner
    # loop started at i+1, so residues[0] appeared on neither side and the
    # first residue in parse order was never compared with anything. Iterating
    # over the full list yields exactly itertools.combinations(residues, 2).
    for i, residue1 in enumerate(residues):
        for _, residue2 in enumerate(residues[i+1:], start=i+1):
            chain1, chain2 = residue1.chain.id, residue2.chain.id

            if (
                # Same residue. The insertion code is part of the test:
                # 100 and 100A share a number but are two residues, and
                # without the code every contact between them was skipped
                # here as self-contact.
                (residue1.resnum == residue2.resnum
                 and residue1.icode == residue2.icode
                 and chain1 == chain2)  # same residue
                or (region and (residue1.resnum not in region or residue2.resnum not in region))  # outside region
                or (chains and (chain1 not in chains or chain2 not in chains))  # chains not allowed
                or (interchain and chain1 == chain2)  # interchain required but same chain
            ):
                continue

            hetero_pair = False
            bonded = None
            if residue1.is_hetero or residue2.is_hetero:
                # Computed bound. The mined Ca-Ca table covers amino-acid pairs
                # only, and there is no equivalent for a prosthetic group or
                # ligand, so the cutoff is derived from this file's coordinates.
                # Single-atom entities (metal ions, chloride) reach this branch
                # too, which the len(atoms) > 1 test below would have excluded.
                hetero_pair = True
                if exclusive and (id(residue1), id(residue2)) in exclusive:
                    continue        # two alternatives of one site (D26)
                anchor1, anchor2 = residue1.anchor, residue2.anchor
                if anchor1 is None or anchor2 is None:
                    continue
                anchor_distance = dist((anchor1.x, anchor1.y, anchor1.z),
                                       (anchor2.x, anchor2.y, anchor2.z))
                if anchor_distance > residue1.radius + residue2.radius + ceiling:
                    continue
                # One sphere per entity is a loose bound for anything large: a
                # 40-atom heme has a radius of ~9 A, so the pair bound reaches
                # 20 A and every residue inside that sphere ran a full
                # atom-by-atom comparison. On 3FZN, 29 hetero residues out of
                # 2,121 cost 4.6 million distance calls and made detection 4.5
                # times slower than the amino-acid-only baseline. The per-atom
                # bounds below cut that down; they are exact, not heuristic -
                # an atom farther from the other entity's centre than that
                # entity's radius plus the largest contact distance cannot
                # reach any of its atoms.
                reach1 = residue2.radius + ceiling
                reach2 = residue1.radius + ceiling
                # atom pairs around a declared covalent bond (D30)
                if covalent:
                    bonded = covalent.get((id(residue1), id(residue2)))

            elif len(residue1.atoms) > 1 and len(residue2.atoms) > 1:
                ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons

                distance_ca = dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
                
                # filter distant residues (static value then specific values)
                if distance_ca > max_ca_distance:
                    continue
                else:
                    key = ''.join(sorted((residue1.resname, residue2.resname)))
                    if distance_ca > updated_distances[key]:
                        continue

            else:
                continue

            # The residue NAME is part of the key. Without it, a position
            # modelled as two different residues - microheterogeneity, where
            # one alternative is deposited as HETATM and the other as ATOM -
            # collapsed to one key, and the electrostatic registry below let
            # only the first of them keep a salt bridge, attractive or
            # repulsive contact. PDB 4Z0Y is the case: His252 and the modified
            # HS8 share position 252 at occupancies 0.49 and 0.51, and 24 of
            # the histidine's contacts disappeared as soon as HS8 was admitted.
            pair_key = (chain1, residue1.resnum, residue1.icode, residue1.resname,
                        chain2, residue2.resnum, residue2.icode, residue2.resname)

            # CHECKING FOR AROMATIC STACKINGS
            # Every ring system of one residue against every ring system of the
            # other. A standard aromatic residue has exactly one, so amino-acid
            # pairs behave precisely as before; a component may have several.
            for ring1, normal1 in residue1.rings:
              for ring2, normal2 in residue2.rings:

                if interface and ring1.entity == ring2.entity:
                    continue
                                
                distance = dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                angle = calc_angle(normal1, normal2)
                
                aromatic_range = categories['aromatic']
                if aromatic_range[0] <= distance <= aromatic_range[1]:
                    if (160 <= angle < 180) or (0 <= angle < 20):
                        stack_type = "-parallel"
                    elif (80 <= angle < 100):
                        stack_type = "-perpendicular"
                    else:
                        stack_type = "-other"

                    contacts_by_pair[pair_key].append({
                        'protein_id': protein.id,
                        'chain1': chain1,
                        'resnum1': residue1.resnum,
                        'icode1': residue1.icode,
                        'alt1': residue1.alt_suffix,
                        'resname1': residue1.resname,
                        'atomname1': ring1.atomname,
                        'chain2': chain2,
                        'resnum2': residue2.resnum,
                        'icode2': residue2.icode,
                        'alt2': residue2.alt_suffix,
                        'resname2': residue2.resname,
                        'atomname2': ring2.atomname,
                        'distance': float(f"{distance:.2f}"),
                        'type': "stacking"+stack_type,
                        'atom1': ring1,
                        'atom2': ring2,
                    })
                    
            # Atoms that can take part. For two amino acids, every typed atom.
            # For a pair with a component, only those within reach of the
            # other's anchor - the exact bound above, tested once per atom
            # rather than once per atom pair. Filtering keeps the order, and
            # the order decides which salt bridge the registry keeps.
            if hetero_pair:
                atoms1 = [a for a in residue1.typed_atoms
                          if dist((a.x, a.y, a.z), (anchor2.x, anchor2.y, anchor2.z)) <= reach1]
                atoms2 = [a for a in residue2.typed_atoms
                          if dist((a.x, a.y, a.z), (anchor1.x, anchor1.y, anchor1.z)) <= reach2]
                pair_types = types_all
            else:
                atoms1, atoms2 = residue1.typed_atoms, residue2.typed_atoms
                pair_types = types_amino

            # skips alpha-helix h-bonds, but only within a chain: residue
            # numbers in different chains are unrelated
            helix = chain1 == chain2 and abs(residue2.resnum - residue1.resnum) <= 3
            cystines = residue1.resname == 'C' and residue2.resname == 'C'
            # per atom pair only where a conformer or a covalent link can exclude one
            guarded = all_altlocs or bool(bonded)

            for atom1 in atoms1:
                key1 = atom1.key
                for atom2 in atoms2:
                    if guarded:
                        # One conformer at a time: altloc A pairs with A or with
                        # an atom that has no alternate, never with B - the PDB's
                        # convention, and what keeps two states from meeting.
                        if (all_altlocs and atom1.altloc and atom2.altloc
                                and atom1.altloc != atom2.altloc):
                            continue
                        if bonded and (atom1.atomname, atom2.atomname) in bonded:
                            continue    # one to three bonds apart across a link (D30)

                    if interface:
                        residue_interface_key = f"{chain1},{residue1.resnum},{residue1.resname}"
                        # if (atom1.entity == atom2.entity) or (residue_interface_key not in interface_res):
                        #     continue
                        if (residue_interface_key not in interface_res):
                            continue

                    distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))

                    # The largest window, not a literal 6: with -d widening
                    # any window past 6 A, everything between 6 A and the
                    # new maximum used to be dropped here - 21% of contacts
                    # with every window widened by 2 A.
                    if distance <= ceiling:

                        for contact_type, low, high, kind, condition in pair_types:

                            if kind == HBOND:
                                if helix:
                                    continue
                                # D29: a disulfide is one interaction. Its two SG
                                # atoms are also donor and acceptor, so the same
                                # pair was written as a hydrogen bond besides DS -
                                # in published v1.6 as well, and only where the
                                # deposited pH left SG both (within 2 of pKa 8.3).
                                if (cystines and key1 == 'C:SG' and atom2.key == 'C:SG'
                                        and disulfide_range[0] <= distance <= disulfide_range[1]):
                                    continue

                            if low <= distance <= high: # fits the range

                                if kind == ELECTROSTATIC:
                                    props1, props2 = atom1.props_elec, atom2.props_elec
                                elif kind == DISULFIDE:
                                    props1, props2 = key1, atom2.key
                                elif kind == METAL:
                                    props1, props2 = atom1, atom2   # defined by element
                                else:
                                    props1, props2 = atom1.props, atom2.props

                                if not condition(props1, props2):
                                    continue
                                    
                                # The registry admits one salt bridge, attractive and
                                # repulsive per residue pair - per conformer, when
                                # every alternate is reported, so B keeps its own.
                                reg_key = (pair_key + (atom1.altloc or atom2.altloc,)
                                           if all_altlocs else pair_key)
                                stored_types = contact_registry[reg_key] # empty set if reg_key is new
                                                                                            
                                if contact_type == 'salt_bridge':
                                    if 'attractive' in stored_types:
                                        # filters attractives out in-place (salt_bridge has priority)
                                        contacts_by_pair[reg_key] = [c for c in contacts_by_pair[reg_key] if c['type'] != 'attractive']
                                        stored_types.remove('attractive')
                                    if 'salt_bridge' in stored_types:
                                        continue # skip adding duplicate
                                    stored_types.add('salt_bridge')
                                elif contact_type == 'attractive':
                                    if 'attractive' in stored_types or 'salt_bridge' in stored_types:
                                        continue
                                    stored_types.add('attractive')
                                elif contact_type == 'repulsive':
                                    if 'repulsive' in stored_types:
                                        continue
                                    stored_types.add('repulsive')
                                else:
                                    stored_types.add(contact_type)
                                    
                                contacts_by_pair[reg_key].append({
                                    'protein_id': protein.id,
                                    'chain1': chain1,
                                    'resnum1': residue1.resnum,
                                    'icode1': residue1.icode,
                                    'alt1': residue1.alt_suffix,
                                    'resname1': residue1.resname,
                                    'atomname1': atom1.atomname,
                                    'chain2': chain2,
                                    'resnum2': residue2.resnum,
                                    'icode2': residue2.icode,
                                    'alt2': residue2.alt_suffix,
                                    'resname2': residue2.resname,
                                    'atomname2': atom2.atomname,
                                    'distance': float(f"{distance:.2f}"),
                                    'type': contact_type,
                                    'atom1': atom1,
                                    'atom2': atom2,
                                    'is_uncertain': kind == ELECTROSTATIC and (atom1.uncertain or atom2.uncertain),
                                    'altloc1': atom1.altloc,
                                    'altloc2': atom2.altloc,
                                })
    
    contacts, count_types = create_contacts(contacts_by_pair)
    
    return contacts, interface_res, count_types


def show_contacts(contacts, altlocs=False):
    """
    Formats and summarizes contact information to be outputted to a file. Only works with the -o flag.

    Args:
        contacts (list): A list of Contact objects of a given protein.

    Returns:
        str: A formatted string summarizing the contact information.
    """
    
    output = []
    
    header = "Chain1,Res1,ResName1,Atom1,Chain2,Res2,ResName2,Atom2,Distance,Type,AA-AA"
    # Two columns at the right-hand end, and only with --all-altlocs, so the
    # default output is byte-identical and left-indexing readers keep working.
    output.append(header + (",Alt1,Alt2" if altlocs else ""))
    for contact in contacts:
        output.append(contact.print_text(altlocs))
        
    return "\n".join(output) # returns as a string to be written directly into the file


# COCaDA-web exclusive
def count_contacts(contacts):
    """
    Formats and returns the number of contacts for each type. Only works with the -o flag.

    Args:
        contacts (list): A list of Contact objects of a given protein.

    Returns:
        list: A list of the number of contacts for each type.
    """
    
    category_counts = {}
    for contact in contacts:
        category = contact.type
        if category in ['stacking-other', 'stacking-parallel', 'stacking-perpendicular']:
            category = 'aromatic'
        category_counts[category] = category_counts.get(category, 0) + 1
        
    expected_keys = ['hydrogen_bond', 'hydrophobic', 'attractive', 'repulsive', 'salt_bridge', 'disulfide_bond', 'aromatic', "uncertain_attractive", "uncertain_repulsive", "uncertain_salt_bridge", 'metal_coordination']
    values = [category_counts.get(key, 0) for key in expected_keys]

    return values


def calc_angle(vector1, vector2):
    """
    Calculates the angle between two ring vectors of aromatic residues

    Args:
        vector1 (tuple): The first vector (x, y, z).
        vector2 (tuple): The second vector (x, y, z).

    Returns:
        float: The angle between the vectors in degrees.
    """
    
    dot_product = dot(vector1, vector2)
    magnitude_product = norm(vector1) * norm(vector2) # normalizes the dot product

    if magnitude_product == 0:
        return 0.0            # a degenerate plane has no meaningful angle

    # Clamp before arccos. The ring normals come from an SVD and are unit
    # vectors, so the ratio is mathematically in [-1, 1] - but for two nearly
    # parallel rings floating point puts it at 1.0000000000000002 about one
    # time in a hundred, and arccos then warns 'invalid value encountered in
    # arccos' and returns nan. A nan angle fails every branch of the
    # parallel/perpendicular test and fell through to 'stacking-other',
    # mislabelling precisely the most parallel stackings there are.
    ratio = min(1.0, max(-1.0, dot_product / magnitude_product))
    angle = arccos(ratio) # angle in radians

    return degrees(angle)


def change_protonation(ph, silent):
    from src.process import log
    
    pka_table = {
        'R': 12.5,
        'K': 10.5,
        'H': 6.0,
        'D': 3.9,
        'E': 4.3,
        'C': 8.3,
        'Y': 10.1,
    }
    
    pH_sensitive_atoms = {
        'R': ['NE', 'CZ', 'NH1', 'NH2'],
        'K': ['NZ'],
        'H': ['ND1', 'NE2'],
        'D': ['OD1', 'OD2'],
        'E': ['OE1', 'OE2'],
        'C': ['SG'],
        'Y': ['OH'],
    }
    
    uncertainty_flags = {}
    local_contact_types = deepcopy(conditions.contact_types)
    
    for key, value in local_contact_types.items():
        resname, atomname = key.split(":")
        if resname in pka_table and atomname in pH_sensitive_atoms.get(resname, []):
            pka = pka_table[resname]
            delta = abs(ph - pka)
            
            original_pos = value[2]
            original_neg = value[3]
            original_donor = value[4]
            original_acceptor = value[5]
            
            new_pos, new_neg = original_pos, original_neg  # Default: no change
            new_donor, new_acceptor = original_donor, original_acceptor
            
            if resname in ['D', 'E', 'Y']:  # Acidic
                if delta < 2.0:
                    new_pos = 0
                    new_neg = 0
                    uncertainty_flags[key] = {'neg': True}
                else:
                    is_deprotonated = ph > pka
                    new_pos = 0
                    new_neg = 1 if is_deprotonated else 0

            elif resname in ['R', 'K', 'H']:  # Basic
                if delta < 2.0:
                    new_pos = 0
                    new_neg = 0
                    uncertainty_flags[key] = {'pos': True}
                else:
                    is_protonated = ph < pka
                    new_pos = 1 if is_protonated else 0
                    new_neg = 0
            
            elif resname == 'C': # Cysteine thiol
                if delta < 2.0:
                    # uncertain protonation: both donor and acceptor could be possible
                    new_donor = 1
                    new_acceptor = 1
                    uncertainty_flags[key] = {'donor': True, 'acceptor': True}
                else:
                    is_deprotonated = ph > pka
                    new_donor = 0 if is_deprotonated else 1
                    new_acceptor = 1 if is_deprotonated else 0
                        
            if (original_pos != new_pos) or (original_neg != new_neg):
                #log(f"pH {ph:.2f} - {key}: (+{original_pos}, -{original_neg}) → (+{new_pos}, -{new_neg}) - pka: {pka}", silent)
                value[2] = new_pos
                value[3] = new_neg
            if (original_donor != new_donor) or (original_acceptor != new_acceptor):
                #log(f"pH {ph:.2f} - {key}: (D{original_donor}, A{original_acceptor}) → (D{new_donor}, A{new_acceptor}) - pka: {pka}", silent)
                value[4] = new_donor
                value[5] = new_acceptor

    #log("\n", silent)
    return uncertainty_flags, local_contact_types


def resolve_uncertainty(name, uncertainty_flags, local_contact_types):
    flags = uncertainty_flags[name]
    new_atom_props = list(local_contact_types[name])    
    
    if flags.get('pos'):
        new_atom_props[2] = 1
    elif flags.get('neg'):
        new_atom_props[3] = 1

    return new_atom_props


def create_contacts(contacts_by_pair):
    count_types = {
        "hydrogen_bond":["HB",0],
        "hydrophobic":["HY",0],
        "attractive":["AT",0],
        "repulsive":["RE",0],
        "salt_bridge":["SB",0],
        "disulfide_bond":["DS",0],
        "stacking":["AS",0],
        "metal_coordination":["MC",0],
    }
    contacts = []

    for contact_list in contacts_by_pair.values():
        for entry in contact_list:
            
            contact_type = entry['type']
            is_uncertain = entry.get('is_uncertain', False)
            if is_uncertain:
                count_types[contact_type][1] += 1
                contact_type = f"uncertain_{contact_type}"
            elif contact_type.startswith('stacking'):
                count_types['stacking'][1] += 1
            else:
                count_types[contact_type][1] += 1
                
            contact = Contact(
                entry['protein_id'],
                entry['chain1'], entry['resnum1'], entry['resname1'], entry['atomname1'],
                entry['protein_id'],
                entry['chain2'], entry['resnum2'], entry['resname2'], entry['atomname2'],
                entry['distance'],
                contact_type,
                entry['atom1'],
                entry['atom2'],
                icode1=entry.get('icode1', ''),
                icode2=entry.get('icode2', ''),
                alt1=entry.get('alt1', ''),
                alt2=entry.get('alt2', ''),
                altloc1=entry.get('altloc1', ''),
                altloc2=entry.get('altloc2', ''),
            )
            contacts.append(contact)

    return contacts, count_types
