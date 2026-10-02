<p align="center">
  <img src="https://github.com/user-attachments/assets/57f17d48-baf4-4bed-97ab-6b817e31dc26" alt="COCaDA_logomenor">
</p>

# COCαDA - A Fast and Scalable Algorithm for Interatomic Contact Detection in Proteins Using Cα Distance Matrices
---

## 🔬 Description
COCαDA (COntact search pruning by Cα Distance Analysis) accelerates the calculation of atomic interactions in proteins, by using a set of fine-tuned **Cα distances** between every pair of aminoacid residues.
The code includes a customized parser for both **PDB** and **CIF** files, with support for large files, residue and atom filtering, and geometric analysis (e.g., centroids and normal vectors for aromatic residues). Users can also define their own contact distance cutoffs via the [`contact_distances.json`](contact_distances.json) configuration file.

### 🔍 Contact types detected:
  - Hydrophobic
  - Hydrogen Bond
  - Attractive
  - Repulsive
  - Disulfide Bond
  - Salt Bridge
  - Aromatic Stacking

## 🚀 Features:
- ⚡ **Fast Processing**: COCαDA averages 2.5x faster processing times against Fixed Cutoffs definitions, and 6x faster against Biopython's `NeighborSearch`.
- 📂 **PDB and CIF Parsing**: Efficient parsing for both PDB and CIF files to extract atomic and residue information.
- 🧼 **Residue Filtering**: Ignores low-quality atoms, water molecules, and problematic/erroneous information.
- 🔬 **Interaction Analysis**: Identifies contacts between specified types of atoms based on predefined conditions.
- ⚙️ **User-defined Distance Cutoffs**: Predefined distance cutoffs for all seven contact types can be easily changed by the user.
- 🔋**User-defined pH customization**: Allows a more nuanced analysis of electrostatic interactions.
- 🌀 **Aromatic Stacking Detection**: Computes centroids and normal vectors for aromatic residues, to determine aromatic stacking contacts.
- 🧠 **Multi-Core Processing**: Parallel batch processing across any combination of CPU cores.
- 📊 **CSV Output**: Clean, structured results ideal for post-analysis and exploratory tasks.

## 📦 Installation

### 🔧 Prerequisites

- Python ≥ 3.x
- `numpy` ≥ 2.0.1
- `psutil` ≥ 6.0.0 *(required only for multi-core mode)*

### 📥 Setup

1. Clone the repository:
   ```sh
   git clone https://github.com/LBS-UFMG/COCaDA.git
   ```

2. Navigate into the project directory:
   ```sh
   cd COCaDA
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
### ▶️ To run COCαDA:

1. Ensure you are in the project directory and the virtual environment is activated (if used).

2. Run COCαDA with the path to your folder or file:
    ```sh
    python cocada.py -f path_to_files/*.{cif/pdb} [-m [CORES]] [-o [OUTPUT_DIR]] [-r] [-d] [-ph] [-c] [-s] [-inter] [-alt] [-pl] [-h]
    ```

### ⚙️ Parameters

| Flag | Description |
|------|-------------|
| `-f`, `--files` | **(Required)** Path(s) to `.cif` or `.pdb` file(s). Wildcards (e.g., `*.cif`) are accepted. |
| `-m`, `--mode` | **(Optional)** Enables Multi-Core mode. <br>• No value = uses all available cores. <br>• `-m X` = use 'X' cores (automatically assigned by OS). <br>• `-m X-Y` = use cores X to Y. <br>• `-m X,Y,Z` = specific multiple cores. |
| `-o`, `--output` | **(Optional)** Outputs detailed results to CSV-formatted files. <br>• No value = saves to `./outputs`. <br>• `-o custom_folder` = saves to specified folder. |
| `-r`, `--region`| **(Optional)** Define only a region of residues to be analyzed. <br>• `-r X-Y` = range of residues from X to Y. <br>• `-r X,Y,Z...` = specific multiple residues. |
| `-d`, `--distances` | **(Optional)** Processes custom contact distances from a `.json` file. By default, uses `contact_distances.json`, which can be altered by the user. <br> Users can also input 14 comma-separated float values representing minimum and maximum values for each of the seven contact types, in the following order: <br> `salt_bridge_min,salt_bridge_max,hydrophobic_min,hydrophobic_max,`<br>`hydrogen_bond_min,hydrogen_bond_max,repulsive_min,repulsive_max,`<br>`attractive_min,attractive_max,disulfide_bond_min,disulfide_bond_max,`<br>`aromatic_min,aromatic_max` |
| `-ph`, `--ph`| **(Optional)** Define pH value to be used. Only electrostatic contacts are affected. By default, uses pH value of 7.4. <br>• `-ph (0-14)` = values must be between 0 and 14. |
| `-c`, `--chains`| **(Optional)** Define only specific chains to be analyzed.<br>• `-c X` = only chain X. <br>• `-c X,Y,Z` = specific multiple chains.|
| `-s`, `--silent`| **(Optional)** Suppresses non-essential console output. |
| `-inter`, `--interchain`| **(Optional)** Calculates only interchain contacts. |
| `-alt`, `--all-altlocs`| **(Optional)** Reports the contacts of every alternate conformation instead of only the highest-occupancy one. Atoms are paired only within one conformer (altloc A with A or with atoms that have none), and two columns, `Alt1` and `Alt2`, are added at the end of each CSV row. |
| `-pl`, `--plddt`| **(Optional)** For predicted models: considers only amino-acid residues whose pLDDT, read from the B-factor column as AlphaFold and other predictors write it, is at least the given value (0–100). A residue's pLDDT is the mean over its atoms. Components (ligands, ions) are never removed. Meaningless for experimental structures, whose B-factor column holds B-factors. |
| `-h`, `--help` | **(Optional)** Shows help and usage instructions. |

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
  one's source and fit) before running COCaDA.
- **Metal coordination keeps its 2.8 Å window.** A transplanted or predicted
  metal site is not refined; in AlphaFill P450s the heme iron sits 1.6–12.5 Å
  from the cysteine, and some real Cys–Fe bonds fall just outside the window.

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
