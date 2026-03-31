# Understanding Protein Files

**A Bioinformatics Tutorial on Protein File Formats, Manipulation, and Visualisation in Python**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Where Do Protein Structures Come From?](#2-where-do-protein-structures-come-from)
3. [Protein File Formats](#3-protein-file-formats)
   - 3.1 [PDBx/mmCIF — The Absolute Standard](#31-pdbxmmcif--the-absolute-standard)
   - 3.2 [Legacy PDB Format — The Retiring Classic](#32-legacy-pdb-format--the-retiring-classic)
   - 3.3 [MMTF — The Speed King for the Web](#33-mmtf--the-speed-king-for-the-web)
   - 3.4 [Quick Format Comparison](#34-quick-format-comparison)
4. [Anatomy of a PDB File](#4-anatomy-of-a-pdb-file)
   - 4.1 [Header Section](#41-header-section)
   - 4.2 [Coordinate Section (ATOM / HETATM Records)](#42-coordinate-section-atom--hetatm-records)
   - 4.3 [Connectivity & Other Records](#43-connectivity--other-records)
5. [Anatomy of a mmCIF File](#5-anatomy-of-a-mmcif-file)
6. [Setting Up Your Python Environment](#6-setting-up-your-python-environment)
7. [Fetching Structures from the PDB](#7-fetching-structures-from-the-pdb)
8. [Parsing Protein Files with Biopython](#8-parsing-protein-files-with-biopython)
   - 8.1 [The SMCRA Hierarchy](#81-the-smcra-hierarchy)
   - 8.2 [Parsing a PDB File](#82-parsing-a-pdb-file)
   - 8.3 [Parsing a mmCIF File](#83-parsing-a-mmcif-file)
   - 8.4 [Navigating the Structure Object](#84-navigating-the-structure-object)
9. [Working with Atomic Coordinates](#9-working-with-atomic-coordinates)
   - 9.1 [Extracting Coordinates as NumPy Arrays](#91-extracting-coordinates-as-numpy-arrays)
   - 9.2 [Calculating Distances Between Atoms](#92-calculating-distances-between-atoms)
   - 9.3 [Selecting Specific Residues or Chains](#93-selecting-specific-residues-or-chains)
10. [Basic Structure Analysis](#10-basic-structure-analysis)
    - 10.1 [Ramachandran Plot](#101-ramachandran-plot)
    - 10.2 [B-Factor Analysis](#102-b-factor-analysis)
    - 10.3 [Contact Maps](#103-contact-maps)
11. [3D Visualisation with py3Dmol](#11-3d-visualisation-with-py3dmol)
    - 11.1 [Basic Rendering](#111-basic-rendering)
    - 11.2 [Highlighting Chains, Residues, and Ligands](#112-highlighting-chains-residues-and-ligands)
    - 11.3 [Colouring by B-Factor](#113-colouring-by-b-factor)
12. [Static Visualisation with Matplotlib](#12-static-visualisation-with-matplotlib)
13. [Writing Modified Structures Back to File](#13-writing-modified-structures-back-to-file)
14. [Summary & Next Steps](#14-summary--next-steps)

---

## 1. Introduction

Before we can run a molecular dynamics simulation, we need to understand the *raw material* — the protein structure file. This tutorial will walk you through:

- **What** protein structure files contain and how they are organised.
- **How** the major file formats (PDB, mmCIF, MMTF) differ.
- **How** to load, inspect, manipulate, and visualise structures in Python.

By the end of this tutorial you will be comfortable opening any protein structure, extracting the information you need, and producing publication-quality figures — all from a Python script or Jupyter notebook.

> **Prerequisites:** Basic familiarity with Python (variables, loops, functions) and a working Python 3.9+ installation.

---

## 2. Where Do Protein Structures Come From?

The **Protein Data Bank (PDB)** ([https://www.rcsb.org](https://www.rcsb.org)) is the single global archive for experimentally determined 3D structures of biological macromolecules. As of 2026, it holds over 220,000 structures solved by:

| Method | What It Measures | Typical Resolution |
|---|---|---|
| **X-ray Crystallography** | Electron density of a crystal | 1.0–3.0 Å |
| **Cryo-Electron Microscopy (cryo-EM)** | Electron density of frozen particles | 2.0–4.0 Å |
| **NMR Spectroscopy** | Internuclear distances in solution | N/A (ensemble) |
| **Neutron Diffraction** | Nuclear density (esp. hydrogen positions) | 1.5–2.5 Å |

Every structure is assigned a **4-character alphanumeric PDB ID** (e.g. `1UBQ` for ubiquitin, `6LU7` for SARS-CoV-2 main protease).

---

## 3. Protein File Formats

### 3.1 PDBx/mmCIF — The Absolute Standard

**Extension:** `.cif`

The **macromolecular Crystallographic Information Framework** (mmCIF / PDBx) is the current master format of the PDB archive. Since 2014, every new deposition must be in mmCIF, and since 2019 it is the *only* format guaranteed to be available for all entries.

**Key features:**
- **Dictionary-driven.** Every data item is defined in the PDBx/mmCIF dictionary, making it self-documenting and machine-readable.
- **No column-width limits.** Unlike the legacy PDB format, there is no 80-character line restriction — so atom names, chain IDs, and residue numbers of any length are supported.
- **Handles large complexes.** Entries with >62 chains or >99,999 atoms (e.g. ribosomes, viral capsids) work natively.

**Example snippet (mmCIF):**

```
data_1UBQ
#
_entry.id   1UBQ
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM   1    N  N   MET A 1   27.340  24.430   2.614  1.00  9.67
ATOM   2    C  CA  MET A 1   26.266  25.413   2.842  1.00 10.38
ATOM   3    C  C   MET A 1   26.913  26.639   3.531  1.00  9.62
...
```

The `loop_` keyword introduces a table. Each `_atom_site.*` line is a column header, and the rows that follow are the data — one row per atom.

### 3.2 Legacy PDB Format — The Retiring Classic

**Extension:** `.pdb`

The legacy PDB format dates back to the 1970s and is a **fixed-column-width** text format (80 characters per line). It is still incredibly widely used in software and tutorials but has significant limitations:

**Key features:**
- Simple, human-readable flat text.
- 80-character line limit creates hard caps:
  - Max 99,999 atoms
  - Max 62 unique chain IDs (A–Z, a–z, 0–9)
  - Max 9,999 residues per chain
- **Being phased out** — large structures are only available in mmCIF.

**Example snippet (PDB):**

```
HEADER    CHROMOSOMAL PROTEIN                     28-JUL-98   1UBQ
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  9.62           C
ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  9.62           O
ATOM      5  CB  MET A   1      25.246  24.880   3.845  1.00 13.77           C
...
HETATM 1232  O   HOH A  77      45.747  30.081  19.708  1.00 12.34           O
END
```

Every column has a specific meaning (we will dissect this in [Section 4](#4-anatomy-of-a-pdb-file)).

### 3.3 MMTF — The Speed King for the Web

**Extension:** `.mmtf`

The **Macromolecular Transmission Format** was developed for fast web delivery. It uses a compact binary encoding (MessagePack) to compress coordinates and is typically **5–20× smaller** than the equivalent PDB/mmCIF file.

**Key features:**
- Binary format — not human-readable.
- Excellent for streaming and web visualisation (used by NGL Viewer, Mol*).
- Supports all data in mmCIF but in a compressed representation.

> **Note:** MMTF development has slowed in favour of BinaryCIF (`.bcif`), another binary encoding of mmCIF. Keep an eye on this space.

### 3.4 Quick Format Comparison

| Feature | Legacy PDB | mmCIF (PDBx) | MMTF / BinaryCIF |
|---|---|---|---|
| **Extension** | `.pdb` | `.cif` | `.mmtf` / `.bcif` |
| **Type** | Fixed-width text | Key-value / tabular text | Binary |
| **Human-Readable** | Yes | Yes | No |
| **Max Atoms** | 99,999 | Unlimited | Unlimited |
| **Max Chains** | 62 | Unlimited | Unlimited |
| **File Size** | Medium | Large | Small |
| **Current Status** | Being retired | Official standard | Niche / web-oriented |

---

## 4. Anatomy of a PDB File

Let's dissect a legacy PDB file section by section. Understanding this format deeply will help you with every other format, because the concepts are shared.

### 4.1 Header Section

The first lines of a PDB file contain metadata:

```
HEADER    CHROMOSOMAL PROTEIN                     28-JUL-98   1UBQ              
TITLE     UBIQUITIN                                                             
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: UBIQUITIN;                                                 
COMPND   3 CHAIN: A;                                                            
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: HOMO SAPIENS;                                   
EXPDTA    X-RAY DIFFRACTION                                                     
REMARK   2 RESOLUTION.    1.80 ANGSTROMS.                                       
```

| Record | Purpose |
|---|---|
| `HEADER` | Classification, deposition date, PDB ID |
| `TITLE` | Short description of the structure |
| `COMPND` | Compound information (molecule name, chains) |
| `SOURCE` | Organism and expression system |
| `EXPDTA` | Experimental method |
| `REMARK 2` | Resolution (for diffraction methods) |

### 4.2 Coordinate Section (ATOM / HETATM Records)

This is the heart of the file. Each `ATOM` line describes one atom in a standard residue (amino acid or nucleotide). `HETATM` lines describe atoms in non-standard residues (ligands, modified residues, water molecules, ions).

```
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N
```

Here is the column-by-column breakdown:

| Columns | Field | Example | Description |
|---|---|---|---|
| 1–6 | Record type | `ATOM` | `ATOM` or `HETATM` |
| 7–11 | Atom serial number | `1` | Unique integer for each atom |
| 13–16 | Atom name | `N` | Chemical name (e.g. `CA` = C-alpha, `CB` = C-beta) |
| 17 | Alternate location | ` ` | For alternate conformations (A, B, etc.) |
| 18–20 | Residue name | `MET` | Three-letter amino acid code |
| 22 | Chain ID | `A` | Single character identifying the chain |
| 23–26 | Residue sequence number | `1` | Position in the chain |
| 27 | Insertion code | ` ` | For inserted residues (e.g. `27A`) |
| 31–38 | X coordinate | `27.340` | In Ångströms |
| 39–46 | Y coordinate | `24.430` | In Ångströms |
| 47–54 | Z coordinate | `2.614` | In Ångströms |
| 55–60 | Occupancy | `1.00` | Fraction of time atom is at this position |
| 61–66 | B-factor | `9.67` | Temperature factor (atomic displacement, Å²) |
| 77–78 | Element symbol | `N` | Chemical element |

**Important concepts:**

- **Occupancy** — A value of `1.00` means the atom is always at this position. Values < 1 indicate alternate conformations that sum to 1.
- **B-factor (Temperature Factor)** — Measures atomic displacement / flexibility. High B-factors = more mobile regions. This is heavily used in MD analysis and validation.
- **Chain ID** — Identifies individual polypeptide chains, nucleic acid strands, or ligand groups within the structure.

### 4.3 Connectivity & Other Records

```
CONECT 1179 1211 1222
TER     641      GLY A  76
END
```

| Record | Purpose |
|---|---|
| `TER` | Marks the end of a polypeptide / nucleic acid chain |
| `CONECT` | Explicit bonds (mainly for ligands / HETATM) |
| `SSBOND` | Disulfide bond records |
| `LINK` | Other covalent links between residues |
| `END` | End of file |

---

## 5. Anatomy of a mmCIF File

A mmCIF file is organised into **data blocks** (starting with `data_XXXX`) containing **categories** (prefixed with `_`). Related items share a category prefix.

```
data_1UBQ
#
_entry.id   1UBQ
#
_cell.length_a           50.840
_cell.length_b           42.770
_cell.length_c           28.950
_cell.angle_alpha        90.00
_cell.angle_beta         90.00
_cell.angle_gamma        90.00
#
loop_
_atom_site.group_PDB         # ATOM or HETATM
_atom_site.id                # Atom serial number
_atom_site.type_symbol       # Element
_atom_site.label_atom_id     # Atom name
_atom_site.label_comp_id     # Residue name
_atom_site.label_asym_id     # Chain ID (internal)
_atom_site.auth_asym_id      # Chain ID (author-assigned)
_atom_site.label_seq_id      # Residue number (internal)
_atom_site.auth_seq_id       # Residue number (author-assigned)
_atom_site.Cartn_x           # X coordinate
_atom_site.Cartn_y           # Y coordinate
_atom_site.Cartn_z           # Z coordinate
_atom_site.occupancy         # Occupancy
_atom_site.B_iso_or_equiv    # B-factor
ATOM 1  N N   MET A A 1 1 27.340 24.430  2.614 1.00  9.67
ATOM 2  C CA  MET A A 1 1 26.266 25.413  2.842 1.00 10.38
...
```

**Key differences from legacy PDB:**
- **Two chain IDs:** `label_asym_id` (assigned by the PDB during processing) and `auth_asym_id` (the author's original label). This causes confusion — always be aware of which one your software uses.
- **Two residue numbering schemes:** Same idea — `label_seq_id` vs `auth_seq_id`.
- **Extensible:** You can add custom categories without breaking parsers.

---

## 6. Setting Up Your Python Environment

Create and activate a dedicated environment:

```bash
# Using conda (recommended)
conda create -n protein_files python=3.11 -y
conda activate protein_files

# Install the key packages
pip install biopython numpy pandas matplotlib seaborn py3Dmol requests
```

**Package overview:**

| Package | Purpose |
|---|---|
| `biopython` | Parsing PDB/mmCIF files, structure manipulation |
| `numpy` | Coordinate maths (distances, angles) |
| `pandas` | Tabular data handling |
| `matplotlib` / `seaborn` | Static 2D plotting |
| `py3Dmol` | Interactive 3D molecular visualisation |
| `requests` | Downloading files from the PDB |

---

## 7. Fetching Structures from the PDB

```python
import requests
import os

def fetch_pdb(pdb_id: str, file_format: str = "pdb", save_dir: str = ".") -> str:
    """
    Download a structure file from the RCSB PDB.
    
    Parameters
    ----------
    pdb_id : str
        4-character PDB ID (e.g. '1UBQ').
    file_format : str
        'pdb' for legacy PDB, 'cif' for mmCIF.
    save_dir : str
        Directory to save the downloaded file.
    
    Returns
    -------
    str
        Path to the saved file.
    """
    pdb_id = pdb_id.upper()
    
    if file_format == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        filename = f"{pdb_id}.pdb"
    elif file_format == "cif":
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        filename = f"{pdb_id}.cif"
    else:
        raise ValueError("file_format must be 'pdb' or 'cif'")
    
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    
    filepath = os.path.join(save_dir, filename)
    with open(filepath, "w") as f:
        f.write(response.text)
    
    print(f"Downloaded {filename} ({len(response.text):,} bytes)")
    return filepath

# Download ubiquitin in both formats
pdb_path = fetch_pdb("1UBQ", file_format="pdb")
cif_path = fetch_pdb("1UBQ", file_format="cif")
```

---

## 8. Parsing Protein Files with Biopython

Biopython's `Bio.PDB` module is the workhorse for structural bioinformatics in Python.

### 8.1 The SMCRA Hierarchy

Biopython represents structures using a **hierarchical data model** called **SMCRA**:

```
Structure
  └── Model          (e.g. NMR ensemble member, or just Model 0 for X-ray)
        └── Chain    (A, B, C, ...)
              └── Residue   (MET1, GLN2, ILE3, ...)
                    └── Atom    (N, CA, C, O, CB, ...)
```

Every level is a Python object you can iterate over, index, and query.

### 8.2 Parsing a PDB File

```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)  # QUIET=True suppresses warnings
structure = parser.get_structure("ubiquitin", "1UBQ.pdb")

# Basic info
print(f"Structure ID: {structure.id}")
print(f"Number of models: {len(structure)}")
print(f"Number of chains: {len(list(structure.get_chains()))}")
print(f"Number of residues: {len(list(structure.get_residues()))}")
print(f"Number of atoms: {len(list(structure.get_atoms()))}")
```

**Output:**
```
Structure ID: ubiquitin
Number of models: 1
Number of chains: 1
Number of residues: 76
Number of atoms: 660
```

### 8.3 Parsing a mmCIF File

```python
from Bio.PDB.MMCIFParser import MMCIFParser

cif_parser = MMCIFParser(QUIET=True)
structure_cif = cif_parser.get_structure("ubiquitin_cif", "1UBQ.cif")

# The resulting Structure object is identical in interface
print(f"Atoms from CIF: {len(list(structure_cif.get_atoms()))}")
```

### 8.4 Navigating the Structure Object

```python
# Access the first model
model = structure[0]

# Access chain A
chain_a = model["A"]

# Access residue 48 (leucine) on chain A
# The residue ID is a tuple: (hetero-flag, sequence_id, insertion_code)
residue_48 = chain_a[(" ", 48, " ")]
print(f"Residue: {residue_48.get_resname()} {residue_48.get_id()[1]}")

# Access the alpha-carbon of that residue
ca_atom = residue_48["CA"]
print(f"CA coordinates: {ca_atom.get_vector()}")
print(f"CA B-factor:    {ca_atom.get_bfactor():.2f}")
print(f"CA occupancy:   {ca_atom.get_occupancy():.2f}")
```

**Iterating over all residues:**

```python
for model in structure:
    for chain in model:
        for residue in chain:
            # Skip water and heteroatoms
            if residue.get_id()[0] == " ":
                resname = residue.get_resname()
                resid = residue.get_id()[1]
                num_atoms = len(residue)
                print(f"Chain {chain.id} | {resname} {resid:>4d} | {num_atoms} atoms")
```

---

## 9. Working with Atomic Coordinates

### 9.1 Extracting Coordinates as NumPy Arrays

```python
import numpy as np

# Extract all CA (alpha-carbon) coordinates
ca_coords = []
ca_residues = []

for residue in structure[0]["A"].get_residues():
    if residue.get_id()[0] != " ":
        continue  # skip heteroatoms / water
    if "CA" in residue:
        ca_coords.append(residue["CA"].get_vector().get_array())
        ca_residues.append(f"{residue.get_resname()}{residue.get_id()[1]}")

ca_coords = np.array(ca_coords)
print(f"Shape: {ca_coords.shape}")  # (76, 3) — 76 residues, 3 coordinates each
print(f"First CA: {ca_residues[0]} at {ca_coords[0]}")
```

### 9.2 Calculating Distances Between Atoms

```python
from Bio.PDB import calc_angle, calc_dihedral

# Distance between two atoms
atom1 = structure[0]["A"][(" ", 1, " ")]["CA"]
atom2 = structure[0]["A"][(" ", 76, " ")]["CA"]

distance = atom1 - atom2  # Biopython overloads the minus operator!
print(f"Distance CA(1) — CA(76): {distance:.2f} Å")

# Euclidean distance with NumPy
dist_np = np.linalg.norm(ca_coords[0] - ca_coords[-1])
print(f"Same distance via NumPy:  {dist_np:.2f} Å")
```

**Distance matrix for all CA atoms:**

```python
from scipy.spatial.distance import pdist, squareform

dist_matrix = squareform(pdist(ca_coords))
print(f"Distance matrix shape: {dist_matrix.shape}")
print(f"Max CA-CA distance: {dist_matrix.max():.2f} Å")
print(f"Min non-zero CA-CA distance: {dist_matrix[dist_matrix > 0].min():.2f} Å")
```

### 9.3 Selecting Specific Residues or Chains

Biopython provides a `Select` class for fine-grained selection:

```python
from Bio.PDB import Selection

# Get all atoms in chain A
chain_a_atoms = Selection.unfold_entities(structure[0]["A"], "A")
print(f"Atoms in chain A: {len(chain_a_atoms)}")

# Get all residues across all chains
all_residues = Selection.unfold_entities(structure[0], "R")
print(f"Total residues: {len(all_residues)}")
```

You can also write your own selector:

```python
from Bio.PDB import Select

class CASelect(Select):
    """Select only alpha-carbon atoms from standard residues."""
    def accept_atom(self, atom):
        if atom.get_name() == "CA" and atom.get_parent().get_id()[0] == " ":
            return True
        return False
```

---

## 10. Basic Structure Analysis

### 10.1 Ramachandran Plot

The **Ramachandran plot** shows the backbone dihedral angles (φ, ψ) for each residue — one of the most important quality checks for a protein structure.

```python
from Bio.PDB.Polypeptide import PPBuilder
import matplotlib.pyplot as plt
import numpy as np

ppb = PPBuilder()
phi_psi = []

for pp in ppb.build_peptides(structure):
    angles = pp.get_phi_psi_list()
    for phi, psi in angles:
        if phi is not None and psi is not None:
            phi_psi.append((np.degrees(phi), np.degrees(psi)))

phi_psi = np.array(phi_psi)

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(phi_psi[:, 0], phi_psi[:, 1], s=15, alpha=0.7, c="steelblue", edgecolors="k", linewidths=0.3)
ax.set_xlabel("φ (degrees)", fontsize=14)
ax.set_ylabel("ψ (degrees)", fontsize=14)
ax.set_title("Ramachandran Plot — 1UBQ", fontsize=16)
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
ax.axhline(0, color="grey", linewidth=0.5)
ax.axvline(0, color="grey", linewidth=0.5)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("ramachandran_1UBQ.png", dpi=150)
plt.show()
```

### 10.2 B-Factor Analysis

B-factors reveal which parts of the structure are flexible vs. rigid.

```python
import matplotlib.pyplot as plt
import numpy as np

# Extract per-residue average B-factor (CA atoms)
residue_ids = []
b_factors = []

for residue in structure[0]["A"].get_residues():
    if residue.get_id()[0] != " ":
        continue
    if "CA" in residue:
        residue_ids.append(residue.get_id()[1])
        b_factors.append(residue["CA"].get_bfactor())

residue_ids = np.array(residue_ids)
b_factors = np.array(b_factors)

fig, ax = plt.subplots(figsize=(12, 4))
ax.bar(residue_ids, b_factors, color="steelblue", edgecolor="none", width=1.0)
ax.set_xlabel("Residue Number", fontsize=13)
ax.set_ylabel("B-factor (Å²)", fontsize=13)
ax.set_title("B-factor per Residue — 1UBQ (CA atoms)", fontsize=15)
ax.set_xlim(residue_ids.min() - 1, residue_ids.max() + 1)
plt.tight_layout()
plt.savefig("bfactor_1UBQ.png", dpi=150)
plt.show()
```

**Interpreting the plot:**
- **Low B-factors** (< 15 Å²): Well-ordered, rigid regions (often β-sheets, core residues).
- **High B-factors** (> 30 Å²): Flexible regions (loops, termini, surface residues).
- Ubiquitin is a compact, stable protein — expect generally low B-factors with peaks at the termini and loops.

### 10.3 Contact Maps

A **contact map** shows which residues are spatially close (typically < 8 Å between CA atoms). This reveals the protein's fold and secondary structure.

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Calculate distance matrix
dist_matrix = squareform(pdist(ca_coords))

# Create binary contact map (threshold = 8 Å)
threshold = 8.0
contact_map = (dist_matrix < threshold).astype(int)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Distance map
im0 = axes[0].imshow(dist_matrix, cmap="viridis_r", origin="lower")
axes[0].set_title("CA Distance Map", fontsize=14)
axes[0].set_xlabel("Residue Index")
axes[0].set_ylabel("Residue Index")
plt.colorbar(im0, ax=axes[0], label="Distance (Å)")

# Contact map
axes[1].imshow(contact_map, cmap="Greys", origin="lower")
axes[1].set_title(f"Contact Map (< {threshold} Å)", fontsize=14)
axes[1].set_xlabel("Residue Index")
axes[1].set_ylabel("Residue Index")

plt.tight_layout()
plt.savefig("contact_map_1UBQ.png", dpi=150)
plt.show()
```

---

## 11. 3D Visualisation with py3Dmol

`py3Dmol` provides interactive 3D molecular graphics inside Jupyter notebooks, powered by 3Dmol.js.

### 11.1 Basic Rendering

```python
import py3Dmol

# Load from the RCSB PDB directly by ID
view = py3Dmol.view(query="pdb:1UBQ", width=800, height=600)

# Cartoon representation
view.setStyle({"cartoon": {"color": "spectrum"}})
view.zoomTo()
view.show()
```

**Other representation styles:**

```python
# Stick model (all atoms)
view.setStyle({"stick": {}})

# Ball and stick
view.setStyle({"stick": {}, "sphere": {"radius": 0.3}})

# Surface representation
view.addSurface(py3Dmol.VDW, {"opacity": 0.7, "color": "white"})
```

### 11.2 Highlighting Chains, Residues, and Ligands

```python
view = py3Dmol.view(query="pdb:6LU7", width=800, height=600)

# Show the whole protein as grey cartoon
view.setStyle({"cartoon": {"color": "lightgrey"}})

# Highlight the active site residues (catalytic dyad: HIS41, CYS145)
view.addStyle(
    {"resi": [41, 145], "chain": "A"},
    {"stick": {"colorscheme": "orangeCarbon", "radius": 0.2}}
)

# Show the bound ligand (N3 inhibitor) as sticks
view.addStyle(
    {"hetflag": True},
    {"stick": {"colorscheme": "greenCarbon", "radius": 0.15}}
)

# Add labels
view.addLabel("HIS 41", {"position": {"resi": 41, "chain": "A"}, "fontSize": 12})
view.addLabel("CYS 145", {"position": {"resi": 145, "chain": "A"}, "fontSize": 12})

view.zoomTo()
view.show()
```

### 11.3 Colouring by B-Factor

```python
view = py3Dmol.view(query="pdb:1UBQ", width=800, height=600)

# Colour by B-factor using a blue-white-red gradient
view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "rwb", "min": 5, "max": 40}}})
view.zoomTo()
view.show()
```

---

## 12. Static Visualisation with Matplotlib

For publication figures where you need full control, you can plot structural features with Matplotlib.

**Backbone trace (CA atoms) coloured by B-factor:**

```python
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Use the CA coordinates and B-factors we extracted earlier
sc = ax.scatter(
    ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2],
    c=b_factors, cmap="coolwarm", s=40, edgecolors="k", linewidths=0.3
)

# Draw backbone trace
ax.plot(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2], color="grey", alpha=0.5, linewidth=1)

ax.set_xlabel("X (Å)")
ax.set_ylabel("Y (Å)")
ax.set_zlabel("Z (Å)")
ax.set_title("1UBQ Backbone Trace (CA) Coloured by B-factor", fontsize=14)

plt.colorbar(sc, ax=ax, label="B-factor (Å²)", shrink=0.6)
plt.tight_layout()
plt.savefig("backbone_trace_1UBQ.png", dpi=150)
plt.show()
```

**Residue property plot (multiple chains):**

```python
import pandas as pd
import matplotlib.pyplot as plt

# Build a DataFrame of per-residue data
data = []
for chain in structure[0]:
    for residue in chain:
        if residue.get_id()[0] != " ":
            continue
        atoms = list(residue.get_atoms())
        avg_b = np.mean([a.get_bfactor() for a in atoms])
        data.append({
            "chain": chain.id,
            "resid": residue.get_id()[1],
            "resname": residue.get_resname(),
            "num_atoms": len(atoms),
            "avg_bfactor": avg_b
        })

df = pd.DataFrame(data)
print(df.head(10))
print(f"\nTotal standard residues: {len(df)}")
print(f"Average B-factor: {df['avg_bfactor'].mean():.2f} Å²")
```

---

## 13. Writing Modified Structures Back to File

After modifying a structure (e.g. removing waters, selecting a chain), you can write it out:

**Writing a PDB file:**

```python
from Bio.PDB import PDBIO

io = PDBIO()
io.set_structure(structure)

# Write the full structure
io.save("1UBQ_output.pdb")

# Write only CA atoms using our custom selector
io.save("1UBQ_CA_only.pdb", select=CASelect())

print("Files written successfully.")
```

**Writing a mmCIF file:**

```python
from Bio.PDB.MMCIFIO import MMCIFIO

cif_io = MMCIFIO()
cif_io.set_structure(structure)
cif_io.save("1UBQ_output.cif")
```

**Selecting a single chain:**

```python
class ChainSelect(Select):
    """Select only a specific chain."""
    def __init__(self, chain_id):
        self.chain_id = chain_id
    
    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

io = PDBIO()
io.set_structure(structure)
io.save("1UBQ_chainA.pdb", select=ChainSelect("A"))
```

**Removing water molecules:**

```python
class NoWaterSelect(Select):
    """Exclude water molecules (HOH)."""
    def accept_residue(self, residue):
        return residue.get_resname() != "HOH"

io.save("1UBQ_no_water.pdb", select=NoWaterSelect())
```

---

## 14. Summary & Next Steps

### What We Covered

| Topic | Key Takeaway |
|---|---|
| **File formats** | mmCIF is the standard; legacy PDB is limited but ubiquitous |
| **PDB anatomy** | `ATOM` records contain coordinates, B-factors, occupancy in fixed columns |
| **mmCIF anatomy** | Dictionary-driven, two naming schemes (label vs. auth) |
| **Biopython parsing** | SMCRA hierarchy: Structure → Model → Chain → Residue → Atom |
| **Coordinate work** | Extract to NumPy arrays; compute distances, dihedrals |
| **Analysis** | Ramachandran plots, B-factor profiles, contact maps |
| **3D visualisation** | py3Dmol for interactive views in Jupyter |
| **Static plots** | Matplotlib/seaborn for publication-quality figures |
| **File I/O** | Write filtered structures with custom `Select` classes |

### Next Steps

- **Protein preparation for MD** → See the main [MD on Apocrita](../MD_On_Apocrita.md) tutorial for preparing your structure for GROMACS (adding hydrogens, solvation, topology generation).
- **MDAnalysis** → For analysing MD trajectories, explore the [MDAnalysis](https://www.mdanalysis.org/) library, which extends many of the concepts here to multi-frame trajectories.
- **ProDy** → For elastic network models and normal mode analysis, check out [ProDy](http://prody.csb.pitt.edu/).
- **PyMOL / ChimeraX** → For production-quality rendered images, these dedicated visualisation tools offer far more control than py3Dmol.

---

*Tutorial written for the QMUL Molecular Dynamics Teaching Project. Last updated: March 2026.*
