#!/usr/bin/env python3
"""
build_ligand_table.py - one command to build COCaDA's non-amino-acid atom
property table from the wwPDB Chemical Component Dictionary.

    python3 tools/build_ligand_table.py --ccd path/to/components.cif

Takes roughly 20-30 minutes for the full dictionary (~51 000 components,
1.33 million heavy atoms). Needs RDKit; COCaDA's runtime does not.

RESUMABLE. Output is appended to a working file one component at a time and
flushed as it goes, so an interrupted run - Ctrl-C, a crash, a machine
restart - loses nothing. Rerunning the same command picks up where it stopped.
Single pass over the input: nothing is re-parsed, and memory stays flat
regardless of dictionary size.

Per-component failures are isolated. The CCD contains entries RDKit cannot
sanitise, entries whose formal charges contradict their own hydrogens, and
entries whose ring perception fails. Each is recorded and skipped; one bad
component never aborts the run. A report of every skip is written alongside
the table.

WHAT IT PRODUCES

    src/ligand_types.dat      the table COCaDA reads (plain text, ~12 MB)
    src/ligand_types.log      every repair and skip, for the record

Commit both. The table is generated data - regenerate it rather than editing
it by hand, or the provenance the validation gate depends on is broken.

COCaDA-ligands (v2.0)
"""

import argparse
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

DEFAULT_OUT = os.path.join(ROOT, "src", "ligand_types.dat")


def existing_ids(path):
    """Component ids already written to a partial file, for resuming."""
    if not os.path.exists(path):
        return set(), 0
    done = set()
    with open(path) as handle:
        for line in handle:
            comp_id, _, _ = line.partition("\t")
            if comp_id:
                done.add(comp_id)
    return done, os.path.getsize(path)


def type_component(component, factory):
    """Every table entry for one component: ({atom name: vector}, ring systems).

    Three layers, each only adding names, never changing one already set:

    1. the molecule the generator types (ccd.for_typing): amino acids as they
       sit in a chain, everything else whole;
    2. heavy leaving atoms that in-chain typing dropped but deposits keep, OXT
       excepted (Component.kept_leaving_atoms), typed from the free molecule;
    3. unambiguous CCD alternate names (Component.alias_names), each with the
       vector of the atom it names. Ring planes are not aliased: the renamed
       atoms are sugar and phosphate atoms, not the aromatic rings that stack.

    Layer 1 is exactly what the table held before layers 2 and 3 existed, so a
    rebuild keeps every earlier entry byte for byte and only adds.
    """
    import ccd as ccd_module
    from generate_atom_types import type_molecule, stacking_rings

    mol, names = ccd_module.for_typing(component)
    typed = type_molecule(mol, factory, lambda i: names.get(i))
    systems = stacking_rings(mol, lambda i: names.get(i))

    atoms = {}
    for atom in mol.GetAtoms():
        name = names.get(atom.GetIdx())
        if name and atom.GetSymbol() != "H":
            atoms[name] = typed.get(name, [0] * 6)
    if not atoms:
        return atoms, systems

    kept = [n for n in component.kept_leaving_atoms() if n not in atoms]
    if kept:
        # The free build records its own repairs; only the outcome is logged,
        # so the log of an unchanged component reads exactly as before.
        before = list(component.repairs)
        try:
            free, free_names = ccd_module.to_rdkit(component)
            free_typed = type_molecule(free, factory, lambda i: free_names.get(i))
            for name in kept:
                atoms[name] = free_typed.get(name, [0] * 6)
            note = f"kept leaving atom(s) typed from the free molecule: {', '.join(kept)}"
        except Exception as exc:                              # noqa: BLE001
            note = f"kept leaving atom(s) {', '.join(kept)} left untyped: {exc}"
        component.repairs[:] = before + [note]

    for alt, name in component.alias_names().items():
        if name in atoms and alt not in atoms:
            atoms[alt] = atoms[name]

    return atoms, systems


def main():
    parser = argparse.ArgumentParser(
        description="Build COCaDA's ligand atom property table from the CCD.")
    parser.add_argument("--ccd", required=True, metavar="PATH",
                        help="components.cif from the wwPDB Chemical Component "
                             "Dictionary, or a directory of per-component files.")
    parser.add_argument("--out", default=DEFAULT_OUT, metavar="PATH",
                        help=f"Output table. Default: {DEFAULT_OUT}")
    parser.add_argument("--restart", action="store_true",
                        help="Discard any partial run and start over.")
    parser.add_argument("--extend", action="store_true",
                        help="Add only components missing from an existing table, "
                             "keeping everything already there. Use after a CCD "
                             "update, or after a reader fix widens coverage.")
    parser.add_argument("--repair", action="store_true",
                        help="Re-type only the components whose stored line cannot "
                             "be parsed, keeping the rest of the table untouched. "
                             "Use after a format fix: it takes seconds rather than "
                             "the full rebuild's twenty minutes.")
    parser.add_argument("--every", type=int, default=1000, metavar="N",
                        help="Progress interval. Default 1000.")
    args = parser.parse_args()

    if not os.path.exists(args.ccd):
        sys.exit(f"error: no such file: {args.ccd}")

    import ccd as ccd_module
    from generate_atom_types import _load_factory
    from ligand_data import encode_name, pack

    partial = args.out + ".partial"
    logfile = os.path.splitext(args.out)[0]
    logfile = logfile[:-4] if logfile.endswith(".dat") else logfile
    logfile += ".log"

    if args.restart:
        for path in (partial, logfile):
            if os.path.exists(path):
                os.remove(path)

    if args.repair:
        if not os.path.exists(args.out):
            sys.exit(f"error: nothing to repair; {args.out} does not exist")
        if os.path.exists(partial):
            os.remove(partial)
        broken = []
        with open(args.out, encoding="ascii") as table, open(partial, "w") as raw:
            for line in table:
                comp_id, _, rest = line.partition("\t")
                atoms_field, _, _ = rest.partition("\t")
                entries = [e for e in atoms_field.split(",") if e]
                ok = all("=" in e and e.rpartition("=")[2].isdigit()
                         for e in entries)
                if ok:
                    raw.write(line)
                else:
                    broken.append(comp_id)
        if not broken:
            os.remove(partial)
            sys.exit(f"nothing to repair: every line in {args.out} parses")
        print(f"repairing {len(broken)} component(s): {', '.join(broken)}")

    if args.extend and os.path.exists(args.out) and not os.path.exists(partial):
        print(f"extending {os.path.basename(args.out)}: "
              f"seeding the working file from it")
        shutil.copyfile(args.out, partial)

    done, size = existing_ids(partial)
    if done:
        print(f"resuming: {len(done):,} components already written "
              f"({size / 1e6:.1f} MB in {os.path.basename(partial)})")

    factory = _load_factory()
    source = (ccd_module.load_dir(args.ccd).values()
              if os.path.isdir(args.ccd) else ccd_module.iter_components(args.ccd))

    started = time.time()
    written = skipped = seen = 0

    with open(partial, "a") as out, open(logfile, "a") as log:
        # The environment, so a later rebuild that differs can be told apart
        # from drift: the 2026-10-01 rebuild could not say which RDKit the
        # shipped table came from without archaeology.
        from rdkit import __version__ as rdkit_version
        log.write(f"BUILD\t{time.strftime('%Y-%m-%d')}\tRDKit {rdkit_version}, "
                  f"{os.path.basename(args.ccd)}\n")
        for component in source:
            seen += 1
            if component.id in done:
                continue

            try:
                atoms, systems = type_component(component, factory)
            except Exception as exc:                          # noqa: BLE001
                log.write(f"SKIP\t{component.id}\t{exc}\n")
                skipped += 1
                continue

            for repair in component.repairs:
                log.write(f"REPAIR\t{component.id}\t{repair}\n")

            if not atoms:
                continue

            # Atom names are escaped: a deposited name may contain this
            # format's own delimiters. See ligand_data for the rule.
            packed = ",".join(f"{encode_name(n)}={pack(atoms[n])}"
                              for n in sorted(atoms))
            rings = ";".join("|".join(encode_name(a) for a in r) for r in systems)
            out.write(f"{component.id}\t{packed}\t{rings}\n")
            out.flush()                    # survive an interruption
            written += 1

            if written % args.every == 0:
                rate = written / (time.time() - started)
                print(f"  {written:,} written, {skipped:,} skipped, "
                      f"{rate:.0f}/s", flush=True)

    total = written + len(done)
    print(f"\ntyped {total:,} components ({skipped:,} skipped) "
          f"in {(time.time() - started) / 60:.1f} min")

    # Uncompressed on purpose: --extend only appends lines, so each CCD update
    # reaches git as the added lines, not as a new copy of the whole table.
    os.replace(partial, args.out)

    mb = os.path.getsize(args.out) / (1024 * 1024)
    print(f"done: {args.out} ({mb:.1f} MB)")
    print(f"log:  {logfile}")
    print("\nNow run the validation gate, which must still pass:")
    print("    python3 tools/generate_atom_types.py --validate")


if __name__ == "__main__":
    sys.exit(main())
