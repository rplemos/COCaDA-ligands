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

    src/ligand_types.dat.gz   the table COCaDA reads (a few MB)
    src/ligand_types.log      every repair and skip, for the record

Commit both. The table is generated data - regenerate it rather than editing
it by hand, or the provenance the validation gate depends on is broken.

COCaDA v1.6+ - ligand support extension
"""

import argparse
import gzip
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

DEFAULT_OUT = os.path.join(ROOT, "src", "ligand_types.dat.gz")


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
    from generate_atom_types import (_load_factory, type_molecule,
                                     stacking_rings)
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
        with gzip.open(args.out, "rt") as gz, open(partial, "w") as raw:
            for line in gz:
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
        with gzip.open(args.out, "rb") as gz, open(partial, "wb") as raw:
            shutil.copyfileobj(gz, raw)

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
        for component in source:
            seen += 1
            if component.id in done:
                continue

            try:
                polymer = component.is_polymer
                mol, names = ccd_module.to_rdkit(
                    component, drop_leaving=polymer, cap_peptide=polymer)
                typed = type_molecule(mol, factory, lambda i: names.get(i))
                systems = stacking_rings(mol, lambda i: names.get(i))
            except Exception as exc:                          # noqa: BLE001
                log.write(f"SKIP\t{component.id}\t{exc}\n")
                skipped += 1
                continue

            for repair in component.repairs:
                log.write(f"REPAIR\t{component.id}\t{repair}\n")

            atoms = {}
            for atom in mol.GetAtoms():
                name = names.get(atom.GetIdx())
                if name and atom.GetSymbol() != "H":
                    atoms[name] = typed.get(name, [0] * 6)
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

    print(f"compressing -> {args.out}")
    with open(partial, "rb") as raw, gzip.open(args.out, "wb", compresslevel=9) as gz:
        shutil.copyfileobj(raw, gz)
    os.remove(partial)

    mb = os.path.getsize(args.out) / (1024 * 1024)
    print(f"done: {args.out} ({mb:.1f} MB)")
    print(f"log:  {logfile}")
    print("\nNow run the validation gate, which must still pass:")
    print("    python3 tools/generate_atom_types.py --validate")


if __name__ == "__main__":
    sys.exit(main())
