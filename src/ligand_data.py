"""
Lazy access to the generated non-amino-acid atom property table.

The wwPDB Chemical Component Dictionary holds ~51 000 components and
1.33 million heavy atoms. Emitting that as a Python dict literal would be a
50 MB source file that every worker process re-parses at import, which is
unacceptable for a tool whose whole point is speed.

Instead the table is a compact, gzipped, line-per-component data file, and
only the components actually present in the structure being processed are
loaded - typically between one and ten. Everything else is never touched.

FILE FORMAT
-----------
One component per line, tab separated:

    HEM<TAB>FE=0,NA=2,C1A=2,...<TAB>C1A|C2A|C3A|C4A|NA;C1B|C2B|...

  field 1  component id
  field 2  atom name = property bitmask, comma separated
  field 3  aromatic ring systems for stacking, atoms joined by '|',
           rings separated by ';'. Empty when the component has none.

The property vector [Hydrophobic, Aromatic, Positive, Negative, Donor,
Acceptor] is packed into one integer, bit i holding flag i. So 2 is
"aromatic only" and 33 is "hydrophobic + acceptor".

COCaDA v1.6+ - ligand support extension
"""

import gzip
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(HERE, "ligand_types.dat.gz")

#: Bit i of a packed mask holds flag i of the property vector.
FLAG_COUNT = 6

_cache_types = {}       # 'HEM:FE' -> [0, 0, 0, 0, 0, 0]
_cache_rings = {}       # 'HEM'    -> [['C1A', ...], ...]
_loaded = set()         # component ids already read from disk
_missing = set()        # component ids looked up and not present


def unpack(mask):
    """Expand a packed integer back into the six-element property vector."""
    return [(mask >> i) & 1 for i in range(FLAG_COUNT)]


def pack(vector):
    """Pack a six-element property vector into a single integer."""
    return sum(bit << i for i, bit in enumerate(vector))


def load_for(comp_ids):
    """Load the components named, returning (types, rings) for them.

    Args:
        comp_ids: iterable of component identifiers, e.g. {'HEM', 'T25'}.

    Returns:
        ({'HEM:FE': [...], ...}, {'HEM': [[...], ...]})

    Components already read are served from cache, and the data file is only
    opened when something is genuinely new. Unknown identifiers are recorded
    so a structure full of unlisted components does not reread the file for
    every one of them.
    """

    wanted = {c for c in comp_ids if c not in _loaded and c not in _missing}

    if wanted and os.path.exists(DATA_FILE):
        found = set()
        with gzip.open(DATA_FILE, "rt") as handle:
            for line in handle:
                comp_id, _, rest = line.partition("\t")
                if comp_id not in wanted:
                    continue
                atoms_field, _, rings_field = rest.partition("\t")

                for entry in atoms_field.split(","):
                    if not entry:
                        continue
                    name, _, mask = entry.rpartition("=")
                    _cache_types[f"{comp_id}:{name}"] = unpack(int(mask))

                rings_field = rings_field.strip()
                if rings_field:
                    _cache_rings[comp_id] = [ring.split("|")
                                             for ring in rings_field.split(";")]

                found.add(comp_id)
                if found == wanted:
                    break

        _loaded.update(found)
        _missing.update(wanted - found)

    ids = set(comp_ids)
    types = {k: v for k, v in _cache_types.items() if k.split(":", 1)[0] in ids}
    rings = {k: v for k, v in _cache_rings.items() if k in ids}
    return types, rings


def available():
    """True when a generated data file is present."""
    return os.path.exists(DATA_FILE)


def reset_cache():
    """Drop everything cached. For tests."""
    _cache_types.clear()
    _cache_rings.clear()
    _loaded.clear()
    _missing.clear()
