<p align="center">
  <img src="https://github.com/user-attachments/assets/57f17d48-baf4-4bed-97ab-6b817e31dc26" alt="COCaDA_logomenor">
</p>

# COCαDA-ligands - Fast Interatomic Contact Detection for Proteins and Every Molecule They Bind
---

## 🔬 Description
COCαDA-ligands extends [COCαDA](https://github.com/rplemos/COCaDA) (COntact search pruning by Cα Distance Analysis) from the 20 standard amino acids to **every component in the wwPDB Chemical Component Dictionary (CCD)**: prosthetic groups, cofactors, metal ions, substrates, inhibitors, nucleotides, sugars and modified residues.

COCαDA accelerates contact detection by discarding residue pairs whose **Cα–Cα distance** exceeds an empirically mined maximum for that pair of residue types, before any atom-level distance is computed. COCαDA-ligands keeps that pruning for amino-acid pairs, and adds a second regime for any pair involving a component: a bound computed from the structure's own coordinates (the component's centroid, its radius and the largest contact distance).

Atoms of components are typed from a table generated from the CCD (`src/ligand_types.dat`, 51,044 components), with the same six properties COCαDA uses for amino acids: hydrophobic, aromatic, positive, negative, hydrogen-bond donor and acceptor. The generator must reproduce COCαDA's curated amino-acid table exactly before its output is accepted (see [Tools](#-tools)). Admitting components does not change amino-acid results: each structure run with and without its `HETATM` records gives the same amino-acid contacts.

### 🔍 Contact types detected:
| Code | Contact type |
|------|--------------|
| `HY` | Hydrophobic |
| `HB` | Hydrogen Bond |
| `AT` | Attractive |
| `RE` | Repulsive |
| `SB` | Salt Bridge |
| `DS` | Disulfide Bond |
| `AS` | Aromatic Stacking |
| `MC` | Metal Coordination *(new)* |

With a user-defined pH (`-ph`), electrostatic contacts of amino acids whose protonation is uncertain at that pH are reported as `uAT`, `uRE` or `uSB`.

## 🚀 Features:
- 🧪 **Every CCD component**: ligands, cofactors, ions and modified residues are read from `HETATM` records and typed from the CCD; nothing is curated by hand.
- 🔩 **Metal coordination**: an eighth contact type, keyed on the element, between a metal and a nitrogen, oxygen or sulfur atom (0–2.8 Å by default), e.g. the cysteine–iron bond of a heme protein.
- 💧 **Coordinating water**: water is excluded, except a water molecule within the coordination window of a metal, which is a ligand of that metal.
- 🌀 **Aromatic stacking for components**: one plane per five- or six-membered ring (a heme gives four pyrrole planes).
- ⚡ **Fast Processing**: the Cα pruning of COCαDA, which averages 2.5x faster processing than fixed-cutoff definitions and 6x faster than Biopython's `NeighborSearch` ([Lemos et al., 2025](#-citation)).
- 📂 **PDB and CIF Parsing**: including insertion codes, quoted mmCIF atom names, selenomethionine, alternate conformations and force-field residue names.
- 🧬 **Predicted models**: AlphaFold DB, ModelArchive and other models are read like deposited structures, with an optional pLDDT cutoff.
- ⚙️ **User-defined Distance Cutoffs**: the distance range of each of the eight contact types can be changed by the user.
- 🔋 **User-defined pH customization**: allows a more nuanced analysis of electrostatic interactions.
- 🧠 **Multi-Core Processing**: parallel batch processing across any combination of CPU cores.
- 📊 **CSV Output**: one row per contact, with a column marking contacts between two amino acids.

## 📦 Installation

### 🔧 Prerequisites

- Python ≥ 3.x
- `numpy` ≥ 2.0.1
- `psutil` ≥ 6.0.0 *(required only for multi-core mode)*
- `rdkit` *(required only to rebuild or update the component table; see [Tools](#-tools))*

### 📥 Setup

1. Clone the repository:
   ```sh
   git clone https://github.com/rplemos/COCaDA-ligands.git
   ```

2. Navigate into the project directory:
   ```sh
   cd COCaDA-ligands
   ```

3. Set up a virtual environment (optional but recommended):
    ```sh
    python -m venv venv
    source venv/bin/activate
    ```

4. Install dependencies:
    ```sh
    pip install -r requirements.txt
    ```

## 🛠️ Usage
### ▶️ To run COCαDA-ligands:

1. Ensure you are in the project directory and the virtual environment is activated (if used).

2. Run COCαDA-ligands with the path to your folder or file:
    ```sh
    python cocada.py -f path_to_files/*.{cif/pdb} [-m [CORES]] [-o [OUTPUT_DIR]] [-r] [-d] [-ph] [-c] [-s] [-inter] [-comp | -aa] [-alt] [-pl] [-h]
    ```

### ⚙️ Parameters

| Flag | Description |
|------|-------------|
| `-f`, `--files` | **(Required)** Path(s) to `.cif` or `.pdb` file(s). Wildcards (e.g., `*.cif`) are accepted. |
| `-m`, `--cores` | **(Optional)** Enables Multi-Core mode. <br>• No value = uses all available cores. <br>• `-m X` = use 'X' cores (automatically assigned by OS). <br>• `-m X-Y` = use cores X to Y. <br>• `-m X,Y,Z` = specific multiple cores. |
| `-o`, `--output` | **(Optional)** Outputs detailed results to CSV-formatted files. <br>• No value = saves to `./outputs`. <br>• `-o custom_folder` = saves to specified folder. |
| `-r`, `--region`| **(Optional)** Define only a region of residues to be analyzed. <br>• `-r X-Y` = range of residues from X to Y. <br>• `-r X,Y,Z...` = specific multiple residues. |
| `-d`, `--distances` | **(Optional)** Processes custom contact distances from a `.json` file. By default, uses [`contact_distances.json`](contact_distances.json), which can be altered by the user. <br> Users can also input 16 comma-separated float values representing minimum and maximum values for each of the eight contact types, in the following order: <br> `salt_bridge_min,salt_bridge_max,hydrophobic_min,hydrophobic_max,`<br>`hydrogen_bond_min,hydrogen_bond_max,repulsive_min,repulsive_max,`<br>`attractive_min,attractive_max,disulfide_bond_min,disulfide_bond_max,`<br>`aromatic_min,aromatic_max,metal_coordination_min,metal_coordination_max` <br> The 14-value form of COCαDA is still accepted, and leaves metal coordination at its default range. |
| `-ph`, `--ph`| **(Optional)** Define pH value to be used. Only electrostatic contacts are affected. By default, uses the pH recorded in the file, or 7.4 if none is recorded. <br>• `-ph (0-14)` = values must be between 0 and 14. |
| `-c`, `--chains`| **(Optional)** Define only specific chains to be analyzed.<br>• `-c X` = only chain X. <br>• `-c X,Y,Z` = specific multiple chains.|
| `-s`, `--silent`| **(Optional)** Suppresses non-essential console output. |
| `-inter`, `--interchain`| **(Optional)** Calculates only interchain contacts. |
| `-comp`, `--components-only`| **(Optional)** Calculates only contacts that involve at least one component (ligand, cofactor, ion, metal-coordinating water), leaving out those between two amino acids: the rows of a full run whose `AA-AA` column is `n`. Faster, since pairs of amino acids are never compared. |
| `-aa`, `--amino-acids-only`| **(Optional)** Calculates only contacts between two amino acids, as COCαDA v1.6 does: the rows of a full run whose `AA-AA` column is `y`. Cannot be combined with `-comp`. |
| `-alt`, `--all-altlocs`| **(Optional)** Reports the contacts of every alternate conformation instead of only the highest-occupancy one. Atoms are paired only within one conformer (altloc A with A or with atoms that have none), and two columns, `Alt1` and `Alt2`, are added at the end of each CSV row. |
| `-pl`, `--plddt`| **(Optional)** For predicted models: considers only amino-acid residues whose pLDDT, read from the B-factor column as AlphaFold and other predictors write it, is at least the given value (0–100). A residue's pLDDT is the mean over its atoms. Components (ligands, ions) are never removed. Meaningless for experimental structures, whose B-factor column holds B-factors. |
| `-h`, `--help` | **(Optional)** Shows help and usage instructions. |

Structures with more than 10,000 residues are skipped and listed in `big.csv` in the output folder.

### 📊 Output

With `-o`, each structure gives `<ID>_contacts.csv`, one row per contact:

```
Chain1,Res1,ResName1,Atom1,Chain2,Res2,ResName2,Atom2,Distance,Type,AA-AA
A,442,C,SG,B,508,HEM,FE,2.35,MC,n
```

- **Residue names**: amino acids in one letter (`C`), components by their CCD id (`HEM`). One-character component ids are written in lower case, so `k` is potassium and `K` lysine.
- **Residue numbers** carry the insertion code (`100A`). Alternatives modelled at one position are told apart by their alternate-location letter (`73A`, `73B`).
- **`AA-AA`**: `y` when both partners are standard amino acids, `n` when a component is involved.
- **Aromatic stacking** rows name the ring centroid `RNG`.

`list.csv` summarises each structure: ID, title, size (residues), number of contacts and pH.

### 🧬 Computational models

AlphaFold DB, ModelArchive and other predicted models are read like deposited
structures; the output is named after the model's own id (`AF-P02766-F1`).

- **Confidence.** Use `-pl` to leave out low-confidence residues, e.g. `-pl 70`.
- **The format follows the extension.** A file is read as PDB when its name ends
  in `.pdb` and as mmCIF otherwise; a PDB-format file saved as `.cif` (as some
  ModelArchive downloads are) reads as empty, so rename it. Compressed files
  (`.cif.gz`) must be decompressed first.
- **Force-field residue names** (Amber `ASH`, `HIP`, `CYX`…, CHARMM `HSD`…,
  GROMACS `HISE`…) are read as the residues they name.
- **AlphaFill models must be filtered before analysis.** An AlphaFill file holds
  the ligands transplanted from every PDB homologue at once, at full occupancy —
  a median 30 components per cytochrome P450 — so their contacts with each other
  dominate the output and are artifacts: these ligands are alternatives, never
  present together. Keep only the transplants you want (AlphaFill reports each
  one's source and fit) before running COCαDA-ligands.
- **Metal coordination keeps its 2.8 Å window.** A transplanted or predicted
  metal site is not refined; in AlphaFill P450s the heme iron sits 1.6–12.5 Å
  from the cysteine, and some real Cys–Fe bonds fall just outside the window.

## 🧰 Tools

The component table is generated data: rebuild or extend it with these tools rather than editing it by hand. They need [RDKit](https://www.rdkit.org/) (`pip install rdkit`); running COCαDA-ligands does not. The CCD itself (`components.cif`) is downloaded from the wwPDB: [files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz](https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz), then decompressed.

| Tool | What it does |
|------|--------------|
| [`tools/update_ccd_table.py`](tools/update_ccd_table.py) | Extends the table with a new CCD release, typing only the components it does not have yet. Every component already in the table must come back unchanged, every new one must load as COCαDA-ligands loads it, and the generator gate must pass; only then, with `--install`, is the table replaced. |
| [`tools/build_ligand_table.py`](tools/build_ligand_table.py) | Builds the whole table from the CCD (resumable), or repairs an existing one (`--repair`). |
| [`tools/generate_atom_types.py`](tools/generate_atom_types.py) | The generator: types each component's atoms from its CCD definition with RDKit and the LUNA feature definitions. `--validate` runs the gate: it must reproduce COCαDA's curated amino-acid table exactly. |
| [`tools/ccd.py`](tools/ccd.py) | Reads the CCD and builds an RDKit molecule for each component. |
| [`tools/aa_templates.py`](tools/aa_templates.py) | The amino-acid templates the gate checks the generator against. |

```sh
# a new CCD release: check first, then install
python tools/update_ccd_table.py --ccd components.cif
python tools/update_ccd_table.py --ccd components.cif --install

# rebuild the whole table, or repair it
python tools/build_ligand_table.py --ccd components.cif
python tools/build_ligand_table.py --ccd components.cif --repair

# the gate
python tools/generate_atom_types.py --validate
```

Every build records in [`src/ligand_types.log`](src/ligand_types.log) the RDKit version, the date, and each component that needed a workaround (`REPAIR`) or could not be typed (`SKIP`). The common workarounds are removing a hydrogen that the CCD lists on a charged atom, and, for definitions that RDKit's valence or aromaticity rules reject (for example, carbons bonded to a metal), typing the component without RDKit's kekulisation and aromaticity perception, keeping the CCD's own aromatic flags. The components that cannot be typed are placeholders and ambiguous codes, such as `UNX` (unknown atom), `DUM` (dummy atom), `ASX` and `GLX`.

The pharmacophore feature definitions in `tools/LUNA.fdef` come from [LUNA](https://github.com/keiserlab/LUNA), MIT License; see [`tools/LUNA.LICENSE`](tools/LUNA.LICENSE).

## 🗂️ Versions

| Branch | Version |
|--------|---------|
| `main` | **v2.0**, the current version |
| `v1.9` | v2.0 without the faster contact detection; same results |
| `v1.6` | COCαDA v1.6, as published in [rplemos/COCaDA](https://github.com/rplemos/COCaDA) |

The history of COCαDA, from its first commit, is kept unchanged.

---

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🧾 Citation

- Main paper:

    **Lemos R.P., Mariano D., Silveira S.A. and de Melo-Minardi R.C.** 
    *COCαDA - A Fast and Scalable Algorithm for Interatomic Contact Detection in Proteins Using Cα Distance Matrices.* Frontiers in Bioinformatics 5:1630078, 2025. DOI: [10.3389/fbinf.2025.1630078](https://doi.org/10.3389/fbinf.2025.1630078).

- Conference paper:

    **Lemos R.P., Mariano D., Silveira S.A. and de Melo-Minardi R.C.** 
    *COCαDA - Large-Scale Protein Interatomic Contact Cutoff Optimization by Cα Distance Matrices.*  Proceedings of the XVII Brazilian Symposium on Bioinformatics (BSB), 17, pp. 59–70, 2024. DOI: [https://doi.org/10.5753/bsb.2024.245545](https://doi.org/10.5753/bsb.2024.245545).

## Contact
For any questions or issues, please contact:

Rafael P. Lemos - PhD Student in Bioinformatics @ Federal University of Minas Gerais, Brazil

Email: rafaellemos@ufmg.br

GitHub: https://github.com/rplemos

## Contributions and Acknowledgements
 - Prof. Raquel Cardoso de Melo Minardi, UFMG;
 - Prof. Sabrina de Azevedo Silveira, UFV;
 - Dr. Diego César Batista Mariano, UFMG;
 - All the 'Laboratory of Bioinformatics and Systems' team.
