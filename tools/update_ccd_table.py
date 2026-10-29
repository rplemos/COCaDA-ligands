#!/usr/bin/env python3
"""
update_ccd_table.py - extend the component table with a new CCD release,
without re-typing anything already in it.

    python3 tools/update_ccd_table.py --ccd ~/ccd/components.cif            # check only
    python3 tools/update_ccd_table.py --ccd ~/ccd/components.cif --install  # and install

Download components.cif yourself (wwPDB: files.wwpdb.org/pub/pdb/data/
monomers/components.cif.gz, then gunzip). Needs RDKit, like the builder.

WHAT IT DOES, in order; any failure stops it and the live table is untouched:

  1. copies src/ligand_types.dat to a work directory and runs
     build_ligand_table.py --extend on the copy: only components the table does
     not have are typed, with the current generator;
  2. checks the result against the table it started from: every component
     already there must come back byte for byte, so nothing that worked can
     change; new components are listed;
  3. loads every new component through ligand_data, as COCaDA would at runtime;
  4. runs the generator gate (167/167), and the unit tests if tests/ is present;
  5. writes update_report.txt in the work directory;
  6. with --install only: replaces src/ligand_types.dat and appends the
     build log to src/ligand_types.log. Commit both, naming the CCD release.

Components whose CCD definition changed since the table was built are NOT
re-typed (that would change working entries); the report lists how many there
are, from _chem_comp.pdbx_modified_date, so a deliberate re-type can be planned.
Components dropped from the CCD stay in the table, also listed.

COCaDA-ligands (v2.0)
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TABLE = os.path.join(ROOT, "src", "ligand_types.dat")
LOG = os.path.join(ROOT, "src", "ligand_types.log")


def read_table(path):
    """{component id: its full line} - the unit the comparison works on."""
    with open(path, encoding="ascii") as handle:
        return {line.partition("\t")[0]: line for line in handle if line.strip()}


def last_build_date(log_path):
    """Date of the most recent BUILD line in the table's log, or None."""
    date = None
    if os.path.exists(log_path):
        with open(log_path) as handle:
            for line in handle:
                if line.startswith("BUILD\t"):
                    date = line.split("\t")[1]
    return date


def ccd_inventory(ccd_path, since):
    """(all component ids, ids modified after `since`) from the new CCD."""
    sys.path.insert(0, HERE)
    import ccd as ccd_module
    ids, changed = set(), set()
    for component in ccd_module.iter_components(ccd_path):
        ids.add(component.id)
        modified = component.scalars.get("_chem_comp.pdbx_modified_date", "")
        if since and modified[:10] > since:
            changed.add(component.id)
    return ids, changed


def run(cmd, report, label):
    report.append(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    tail = (result.stdout + result.stderr).strip().splitlines()[-3:]
    report.extend(f"    {line}" for line in tail)
    if result.returncode != 0:
        raise SystemExit(f"{label} failed (exit {result.returncode}); live table untouched")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ccd", required=True, metavar="PATH",
                    help="components.cif of the new CCD release (uncompressed)")
    ap.add_argument("--table", default=TABLE, metavar="PATH",
                    help=f"Table to extend. Default: {TABLE}")
    ap.add_argument("--workdir", default=None, metavar="DIR",
                    help="Where the extended copy and the report go. "
                         "Default: a new temporary directory.")
    ap.add_argument("--install", action="store_true",
                    help="Replace the live table if every check passes.")
    args = ap.parse_args()

    if not os.path.exists(args.ccd):
        sys.exit(f"error: no such file: {args.ccd}")
    work = args.workdir or tempfile.mkdtemp(prefix="ccd_update_")
    os.makedirs(work, exist_ok=True)
    candidate = os.path.join(work, "ligand_types.dat")
    candidate_log = os.path.join(work, "ligand_types.log")
    report = [f"CCD update {time.strftime('%Y-%m-%d %H:%M')}",
              f"ccd:   {os.path.abspath(args.ccd)}",
              f"table: {os.path.abspath(args.table)}", ""]

    # 1. extend a copy
    for stale in (candidate + ".partial", candidate_log):
        if os.path.exists(stale):
            os.remove(stale)
    shutil.copyfile(args.table, candidate)
    print("extending a copy of the table (only new components are typed)...")
    run([sys.executable, os.path.join(HERE, "build_ligand_table.py"),
         "--ccd", args.ccd, "--out", candidate, "--extend", "--every", "500"],
        report, "the extension")

    # 2. nothing already in the table may change
    before, after = read_table(args.table), read_table(candidate)
    changed = sorted(c for c in before if after.get(c) != before[c])
    if changed:
        report.append(f"\nFAILED: {len(changed)} existing component(s) changed: "
                      f"{', '.join(changed[:20])}")
        open(os.path.join(work, "update_report.txt"), "w").write("\n".join(report) + "\n")
        sys.exit(f"existing entries changed ({len(changed)}); live table untouched. "
                 f"Report: {work}/update_report.txt")
    added = sorted(set(after) - set(before))
    skipped = []
    if os.path.exists(candidate_log):
        with open(candidate_log) as handle:
            skipped = [line.split("\t")[1] for line in handle if line.startswith("SKIP\t")]
    report += ["", f"components before: {len(before):,}",
               f"components after:  {len(after):,}",
               f"added:             {len(added):,}",
               f"could not be typed: {len(skipped):,}  {' '.join(skipped[:30])}",
               "existing components unchanged: all"]

    # 3. every new component loads as COCaDA would load it
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import ligand_data
    ligand_data.DATA_FILE = candidate
    ligand_data.reset_cache()
    types, _rings = ligand_data.load_for(set(added))
    loaded = {key.split(":", 1)[0] for key in types}
    missing = [c for c in added if c not in loaded]
    if missing:
        sys.exit(f"{len(missing)} new component(s) do not load: {missing[:10]}; "
                 f"live table untouched")
    report.append(f"new components loaded through ligand_data: {len(added):,}")

    # 4. the gate and the tests
    run([sys.executable, os.path.join(HERE, "generate_atom_types.py"), "--validate"],
        report, "the generator gate")
    tests = os.path.join(ROOT, "tests")
    if os.path.isdir(tests):
        run([sys.executable, "-m", "pytest", tests, "-q"], report, "the unit tests")
    else:
        report.append("unit tests: no tests/ in this checkout, not run")

    # information: definitions changed upstream, components gone from the CCD
    since = last_build_date(LOG)
    ccd_ids, modified = ccd_inventory(args.ccd, since)
    stale = sorted(modified & set(before))
    gone = sorted(set(before) - ccd_ids)
    report += ["", f"table last built: {since or 'unknown'}",
               f"already in the table but modified in the CCD since: {len(stale):,} "
               f"(not re-typed){'  e.g. ' + ' '.join(stale[:15]) if stale else ''}",
               f"in the table but no longer in the CCD: {len(gone):,} (kept)",
               "", "new components:", "  " + " ".join(added)]

    report_path = os.path.join(work, "update_report.txt")
    with open(report_path, "w") as handle:
        handle.write("\n".join(report) + "\n")
    print(f"\nall checks passed: +{len(added):,} components, "
          f"{len(skipped)} could not be typed. Report: {report_path}")

    # 6. install
    if not args.install:
        print("not installed (dry run). Re-run with --install to replace the live table.")
        return 0
    shutil.copyfile(candidate, TABLE)
    if os.path.exists(candidate_log):
        with open(candidate_log) as src, open(LOG, "a") as dst:
            dst.write(src.read())
    print(f"installed: {TABLE}\nlog appended: {LOG}\n"
          "Commit both, naming the CCD release, and rerun a validation sample.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
