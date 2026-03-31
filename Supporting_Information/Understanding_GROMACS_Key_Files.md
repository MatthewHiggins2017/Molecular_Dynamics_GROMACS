# Understanding GROMACS Key Files & Force Field Models

**A Bioinformatics Reference Sheet for Molecular Dynamics with GROMACS**

---

## Table of Contents

**Part A — GROMACS Key Files**

1. [Overview: The GROMACS File Ecosystem](#1-overview-the-gromacs-file-ecosystem)
2. [Structure Files](#2-structure-files)
   - 2.1 [`.gro` — GROMACS Coordinate File](#21-gro--gromacs-coordinate-file)
   - 2.2 [`.pdb` — Protein Data Bank File (as used by GROMACS)](#22-pdb--protein-data-bank-file-as-used-by-gromacs)
3. [Topology Files](#3-topology-files)
   - 3.1 [`.top` — System Topology](#31-top--system-topology)
   - 3.2 [`.itp` — Include Topology (Molecule Definitions)](#32-itp--include-topology-molecule-definitions)
   - 3.3 [`posre.itp` — Position Restraint File](#33-posreitp--position-restraint-file)
4. [Parameter Files](#4-parameter-files)
   - 4.1 [`.mdp` — Molecular Dynamics Parameters](#41-mdp--molecular-dynamics-parameters)
5. [Run & Binary Files](#5-run--binary-files)
   - 5.1 [`.tpr` — Portable Run Input](#51-tpr--portable-run-input)
6. [Trajectory Files](#6-trajectory-files)
   - 6.1 [`.xtc` — Compressed Trajectory](#61-xtc--compressed-trajectory)
   - 6.2 [`.trr` — Full-Precision Trajectory](#62-trr--full-precision-trajectory)
   - 6.3 [`.edr` — Energy File](#63-edr--energy-file)
7. [Index Files](#7-index-files)
   - 7.1 [`.ndx` — Index File](#71-ndx--index-file)
8. [Checkpoint Files](#8-checkpoint-files)
   - 8.1 [`.cpt` — Checkpoint File](#81-cpt--checkpoint-file)
9. [Log Files](#9-log-files)
   - 9.1 [`.log` — Simulation Log](#91-log--simulation-log)

**Part B — Force Field Models**

10. [What Is a Force Field?](#10-what-is-a-force-field)
11. [The Functional Form](#11-the-functional-form)
12. [Force Field Directory Layout in GROMACS](#12-force-field-directory-layout-in-gromacs)
13. [All-Atom Force Fields](#13-all-atom-force-fields)
    - 13.1 [AMBER Family (ff99SB, ff14SB, ff19SB)](#131-amber-family-ff99sb-ff14sb-ff19sb)
    - 13.2 [CHARMM Family (CHARMM36 / CHARMM36m)](#132-charmm-family-charmm36--charmm36m)
    - 13.3 [OPLS-AA Family (OPLS-AA, OPLS-AA/M)](#133-opls-aa-family-opls-aa-opls-aam)
    - 13.4 [GROMOS Family (54A7, 54A8)](#134-gromos-family-54a7-54a8)
14. [Coarse-Grained Force Fields](#14-coarse-grained-force-fields)
    - 14.1 [MARTINI (v2 / v3)](#141-martini-v2--v3)
15. [Specialised & Polarisable Force Fields](#15-specialised--polarisable-force-fields)
    - 15.1 [AMOEBA](#151-amoeba)
    - 15.2 [Drude Polarisable Force Field](#152-drude-polarisable-force-field)
16. [Water Models](#16-water-models)
17. [Force Field Comparison Summary](#17-force-field-comparison-summary)
18. [How to Choose a Force Field](#18-how-to-choose-a-force-field)
19. [Summary & Next Steps](#19-summary--next-steps)

---

# Part A — GROMACS Key Files

---

## 1. Overview: The GROMACS File Ecosystem

A typical GROMACS workflow generates and consumes a large number of files. Understanding what each one contains — and when it is created — is essential for troubleshooting and reproducibility.

```
                    ┌──────────┐
  Input PDB ──────► │ pdb2gmx  │──────► .gro  (processed structure)
                    │          │──────► .top  (system topology)
                    │          │──────► .itp  (molecule topology + position restraints)
                    └──────────┘
                         │
         .mdp ──────►┌──────────┐
         .gro ──────►│  grompp  │──────► .tpr  (portable run input)
         .top ──────►└──────────┘
                         │
                    ┌──────────┐
         .tpr ────►│   mdrun  │──────► .xtc / .trr  (trajectory)
                    │          │──────► .edr          (energies)
                    │          │──────► .log          (log)
                    │          │──────► .cpt          (checkpoint)
                    └──────────┘
```

**Quick reference table:**

| Extension | Type | Human-Readable | Created By | Purpose |
|---|---|---|---|---|
| `.gro` | Structure | Yes | `pdb2gmx`, `editconf`, `solvate` | Atom coordinates (snapshot) |
| `.pdb` | Structure | Yes | Various | Atom coordinates (PDB format) |
| `.top` | Topology | Yes | `pdb2gmx`, manual | System-level bonded/nonbonded parameters |
| `.itp` | Topology | Yes | `pdb2gmx`, manual | Molecule-level parameters (included by `.top`) |
| `.mdp` | Parameters | Yes | User (manual) | Simulation settings |
| `.tpr` | Run input | **No** (binary) | `grompp` | Complete simulation input |
| `.xtc` | Trajectory | **No** (binary) | `mdrun` | Compressed coordinate trajectory |
| `.trr` | Trajectory | **No** (binary) | `mdrun` | Full-precision trajectory (coords + velocities + forces) |
| `.edr` | Energies | **No** (binary) | `mdrun` | Energy terms over time |
| `.ndx` | Index | Yes | `make_ndx`, manual | Custom atom/residue groups |
| `.cpt` | Checkpoint | **No** (binary) | `mdrun` | Restart information |
| `.log` | Log | Yes | `mdrun` | Performance stats, warnings, settings |

---

## 2. Structure Files

### 2.1 `.gro` — GROMACS Coordinate File

**Origin:** Generated by `gmx pdb2gmx` (from a PDB/mmCIF input), `gmx editconf`, `gmx solvate`, `gmx genion`, or `gmx trjconv`.

The `.gro` format is GROMACS's native, fixed-column-width coordinate format. It stores a **single snapshot** of the system.

**Example file (`protein_processed.gro`):**

```
Ubiquitin in water, after energy minimisation
 8SEQ
    1MET      N    1   2.734   2.443   0.261
    1MET     CA    2   2.627   2.541   0.284
    1MET      C    3   2.691   2.664   0.353
    1MET      O    4   2.789   2.646   0.426
    1MET     CB    5   2.525   2.488   0.385
    2GLN      N    6   2.643   2.782   0.323
    ...
   5.0840   4.2770   2.8950
```

**Line-by-line breakdown:**

| Line | Content | Notes |
|---|---|---|
| **1** | Title string | Free-text description |
| **2** | Total number of atoms | Single integer |
| **3 → N+2** | Atom records | Fixed-width columns (see below) |
| **Last** | Box vectors | 3 (cubic) or 9 (triclinic) values in **nm** |

**Atom record columns:**

| Columns | Width | Field | Units |
|---|---|---|---|
| 1–5 | 5 | Residue number | — |
| 6–10 | 5 | Residue name | — |
| 11–15 | 5 | Atom name | — |
| 16–20 | 5 | Atom serial number | — |
| 21–28 | 8 | X coordinate | **nm** (not Å!) |
| 29–36 | 8 | Y coordinate | nm |
| 37–44 | 8 | Z coordinate | nm |
| 45–52 | 8 | X velocity (optional) | nm/ps |
| 53–60 | 8 | Y velocity (optional) | nm/ps |
| 61–68 | 8 | Z velocity (optional) | nm/ps |

> **Critical note:** GROMACS uses **nanometres (nm)**, not Ångströms (Å). 1 nm = 10 Å. This is the #1 source of unit confusion for newcomers.

**Limitations:**
- Residue numbers wrap at 99999.
- Atom serial numbers wrap at 99999.
- Fixed precision (3 decimal places for coordinates = 0.001 nm = 0.01 Å).

### 2.2 `.pdb` — Protein Data Bank File (as used by GROMACS)

GROMACS can read and write standard PDB files (see the [Understanding Protein Files](Understanding_Protein_Files.md) tutorial for the full format breakdown). PDB files are used as:

- **Input** to `pdb2gmx` (the starting point of most workflows).
- **Output** from `trjconv -o output.pdb` for visualisation in PyMOL/ChimeraX.

GROMACS writes PDB files in Ångströms (converting internally from nm).

---

## 3. Topology Files

### 3.1 `.top` — System Topology

**Origin:** Generated by `gmx pdb2gmx` and then manually edited (e.g. after solvation with `gmx solvate` or ion addition with `gmx genion`).

The `.top` file is the **master topology** for the entire simulation system. It defines every atom, bond, angle, dihedral, and nonbonded interaction via a combination of inline definitions and `#include` directives.

**Example file (`topol.top`):**

```
;
;  Topology for ubiquitin in water
;  Generated by GROMACS pdb2gmx
;

; Force field parameters
#include "amber99sb-ildn.ff/forcefield.itp"

; Protein topology
#include "Protein_chain_A.itp"

; Position restraints (activated when -DPOSRES is passed to grompp)
#ifdef POSRES
#include "posre.itp"
#endif

; Water model topology
#include "amber99sb-ildn.ff/tip3p.itp"

; Ion parameters
#include "amber99sb-ildn.ff/ions.itp"

[ system ]
; Name
Ubiquitin in water

[ molecules ]
; Compound       #mols
Protein_chain_A     1
SOL             10234
NA                  8
CL                  6
```

**Section breakdown:**

| Section / Directive | Purpose |
|---|---|
| `#include "...forcefield.itp"` | Loads the force field atom types, nonbonded parameters, and combination rules |
| `#include "Protein_chain_A.itp"` | Loads the molecule-specific topology (bonds, angles, dihedrals, charges) |
| `#ifdef POSRES` | Conditional block — position restraints are only active if you pass `-DPOSRES` to `grompp` |
| `#include "...tip3p.itp"` | Water model definition |
| `#include "...ions.itp"` | Ion parameters (Na⁺, Cl⁻, etc.) |
| `[ system ]` | System name (cosmetic) |
| `[ molecules ]` | **Critical:** lists every molecule type and count, in the **exact order** they appear in the coordinate file |

> **Common pitfall:** The `[ molecules ]` section must match the coordinate file exactly. If you add water with `gmx solvate`, you must update the `SOL` count. `solvate` does this automatically, but `genion` requires you to update the ion and water counts.

### 3.2 `.itp` — Include Topology (Molecule Definitions)

**Origin:** Generated by `gmx pdb2gmx` for each molecule; can also be written manually or generated by tools like CGenFF, ACPYPE, or ATB for ligands.

The `.itp` file defines a single molecule type in full detail. It is `#include`-d by the `.top` file.

**Example structure (`Protein_chain_A.itp`):**

```
[ moleculetype ]
; Name            nrexcl
Protein_chain_A     3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass
     1         N3      1     MET      N      1    0.14940  14.010000
     2          H      1     MET     H1      2    0.09760   1.008000
     3          H      1     MET     H2      3    0.09760   1.008000
     4          H      1     MET     H3      4    0.09760   1.008000
     5         CT      1     MET     CA      5    0.02210  12.010000
     6         HP      1     MET     HA      6    0.11160   1.008000
    ...

[ bonds ]
;  ai    aj funct            c0            c1
    1     2     1
    1     3     1
    1     4     1
    1     5     1
    5     6     1
    5     7     1
    5    15     1
    ...

[ pairs ]
;  ai    aj funct
    1     8     1
    1    16     1
    ...

[ angles ]
;  ai    aj    ak funct            c0            c1
    2     1     3     1
    2     1     4     1
    2     1     5     1
    ...

[ dihedrals ] ; proper dihedrals
;  ai    aj    ak    al funct     c0     c1     c2     c3     c4     c5
    2     1     5     6     9
    2     1     5     7     9
    ...

[ dihedrals ] ; improper dihedrals
;  ai    aj    ak    al funct     c0     c1     c2     c3
   15     5    17    16     4
    ...
```

**Section breakdown:**

| Section | Purpose |
|---|---|
| `[ moleculetype ]` | Molecule name and number of excluded neighbours (`nrexcl` — typically 3 for proteins, meaning 1-2 and 1-3 nonbonded interactions are excluded) |
| `[ atoms ]` | Every atom: index, force field type, residue, partial charge, mass |
| `[ bonds ]` | Bonded pairs with function type (1 = harmonic) |
| `[ pairs ]` | 1-4 nonbonded pairs (special scaling for Lennard-Jones and Coulomb) |
| `[ angles ]` | Angle triplets with function type |
| `[ dihedrals ]` | Proper (torsional) and improper (out-of-plane) dihedrals |

### 3.3 `posre.itp` — Position Restraint File

**Origin:** Generated by `gmx pdb2gmx` alongside the molecule topology.

Position restraints pin heavy atoms to their initial positions during equilibration, preventing the protein from moving while the solvent relaxes.

**Example (`posre.itp`):**

```
; Position restraint file for Protein_chain_A
; Generated by pdb2gmx

[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
   5    1       1000       1000       1000
   7    1       1000       1000       1000
  15    1       1000       1000       1000
  17    1       1000       1000       1000
  ...
```

| Column | Meaning |
|---|---|
| `i` | Atom index (must match the `[ atoms ]` numbering in the `.itp`) |
| `funct` | Function type (1 = harmonic) |
| `fcx`, `fcy`, `fcz` | Force constant in each dimension (kJ mol⁻¹ nm⁻²) |

A force constant of **1000 kJ mol⁻¹ nm⁻²** is the default. This is strong enough to restrain the backbone during equilibration but can be reduced for gradual release protocols.

---

## 4. Parameter Files

### 4.1 `.mdp` — Molecular Dynamics Parameters

**Origin:** Written manually by the user (or adapted from templates).

The `.mdp` (Molecular Dynamics Parameter) file controls **everything** about how the simulation is run. It is consumed by `gmx grompp` to build the `.tpr`.

**Example — Energy Minimisation (`em.mdp`):**

```
; Energy minimisation parameters
integrator  = steep         ; Steepest descent minimisation
emtol       = 1000.0        ; Stop when max force < 1000 kJ/mol/nm
emstep      = 0.01          ; Initial step size (nm)
nsteps      = 50000         ; Max number of minimisation steps

; Neighbour searching
nstlist     = 10            ; Update neighbour list every 10 steps
cutoff-scheme = Verlet      ; Verlet cutoff scheme
ns_type     = grid          ; Grid-based neighbour searching
rcoulomb    = 1.0           ; Short-range Coulomb cutoff (nm)
rvdw        = 1.0           ; Short-range van der Waals cutoff (nm)

; Electrostatics
coulombtype = PME           ; Particle Mesh Ewald for long-range electrostatics
pme_order   = 4             ; PME interpolation order
fourierspacing = 0.16       ; Grid spacing for FFT (nm)
```

**Example — Production MD (`md.mdp`):**

```
; Production MD parameters
integrator  = md            ; Leap-frog integrator
dt          = 0.002         ; Time step: 2 fs
nsteps      = 5000000       ; Total: 10 ns (5,000,000 × 0.002 ps)

; Output control
nstxout-compressed = 5000   ; Write compressed coords every 10 ps
nstlog      = 5000          ; Write log every 10 ps
nstenergy   = 5000          ; Write energies every 10 ps

; Bond constraints
continuation = yes          ; Continuing from NVT/NPT equilibration
constraint_algorithm = lincs ; LINCS for bond constraints
constraints = h-bonds       ; Constrain all bonds involving hydrogen
lincs_iter  = 1             ; LINCS iteration order
lincs_order = 4             ; LINCS expansion order

; Neighbour searching
nstlist     = 10
cutoff-scheme = Verlet
rcoulomb    = 1.0
rvdw        = 1.0

; Electrostatics
coulombtype = PME
pme_order   = 4
fourierspacing = 0.16

; Temperature coupling
tcoupl      = V-rescale     ; Modified Berendsen thermostat
tc-grps     = Protein Non-Protein
tau_t       = 0.1     0.1   ; Time constant (ps)
ref_t       = 300     300   ; Reference temperature (K)

; Pressure coupling
pcoupl      = Parrinello-Rahman  ; Pressure coupling method
pcoupltype  = isotropic
tau_p       = 2.0           ; Time constant (ps)
ref_p       = 1.0           ; Reference pressure (bar)
compressibility = 4.5e-5    ; Water compressibility (bar⁻¹)

; Periodic boundary conditions
pbc         = xyz           ; 3D periodic boundaries

; Dispersion correction
DispCorr    = EnerPres      ; Long-range dispersion corrections for energy & pressure
```

**Key parameter groups:**

| Group | Key Parameters | What They Control |
|---|---|---|
| **Integrator** | `integrator`, `dt`, `nsteps` | Algorithm, time step, total simulation length |
| **Output** | `nstxout-compressed`, `nstenergy`, `nstlog` | How often to write trajectory/energy/log |
| **Constraints** | `constraints`, `constraint_algorithm` | Which bonds to constrain (enables 2 fs time step) |
| **Electrostatics** | `coulombtype`, `rcoulomb`, `pme_order` | How long-range electrostatics are handled |
| **van der Waals** | `rvdw`, `DispCorr` | VdW cutoff and dispersion corrections |
| **Temperature** | `tcoupl`, `ref_t`, `tau_t` | Thermostat type, target temperature, coupling strength |
| **Pressure** | `pcoupl`, `ref_p`, `tau_p` | Barostat type, target pressure, coupling strength |
| **PBC** | `pbc` | Periodic boundary conditions |

---

## 5. Run & Binary Files

### 5.1 `.tpr` — Portable Run Input

**Origin:** Generated by `gmx grompp` from the `.mdp`, `.gro`, and `.top` files.

The `.tpr` is a **binary file** that bundles everything needed for `mdrun`:
- All force field parameters (from `.top` / `.itp`)
- Simulation parameters (from `.mdp`)
- Starting coordinates (from `.gro`)
- Starting velocities (if present)

It is self-contained and portable — you can copy it to a different machine and run the simulation without any other files.

**Inspecting a `.tpr`:**

```bash
# Print a human-readable summary
gmx dump -s em.tpr | head -100

# Check simulation parameters
gmx dump -s md.tpr | grep -E "nsteps|dt|ref-t|ref-p|coulombtype"
```

---

## 6. Trajectory Files

### 6.1 `.xtc` — Compressed Trajectory

**Origin:** Written by `gmx mdrun` during the simulation.

The `.xtc` format uses lossy compression to store coordinates with reduced precision (typically 0.001 nm = 0.01 Å). It is the standard trajectory format for analysis.

| Property | Value |
|---|---|
| **Contains** | Coordinates only |
| **Precision** | ~0.001 nm (configurable) |
| **Typical size** | ~10 KB per frame per 1000 atoms |
| **Velocities** | No |
| **Forces** | No |

### 6.2 `.trr` — Full-Precision Trajectory

**Origin:** Written by `gmx mdrun` if requested (via `nstxout`, `nstvout`, `nstfout` in `.mdp`).

The `.trr` format stores coordinates, velocities, and/or forces at full (double) precision. Files are **much larger** than `.xtc`.

| Property | Value |
|---|---|
| **Contains** | Coordinates, velocities, forces |
| **Precision** | Full (double) |
| **Typical size** | ~50–100 KB per frame per 1000 atoms |
| **Use case** | When you need velocities/forces, or for restart |

### 6.3 `.edr` — Energy File

**Origin:** Written by `gmx mdrun`.

A binary file containing the time series of all energy terms (kinetic, potential, temperature, pressure, density, etc.). Analysed with `gmx energy`:

```bash
# Extract temperature and potential energy
echo "Temperature\nPotential\n\n" | gmx energy -f md.edr -o energy.xvg
```

---

## 7. Index Files

### 7.1 `.ndx` — Index File

**Origin:** Generated by `gmx make_ndx` or written manually.

Index files define custom groups of atoms for analysis or simulation settings. GROMACS has built-in default groups (`Protein`, `Backbone`, `Water`, `System`, etc.), but you often need custom groups.

**Example (`index.ndx`):**

```
[ Protein ]
   1    2    3    4    5    6    7    8    9   10
  11   12   13   14   15   16   17   18   19   20
  ...

[ Backbone ]
   1    5   15   17   28   30   ...

[ Active_Site ]
 321  322  323  324  325  640  641  642  643  644

[ Ligand ]
 1232 1233 1234 1235 1236 1237 1238
```

Each group starts with `[ GroupName ]` followed by atom indices (1-based).

```bash
# Create an index file interactively
gmx make_ndx -f system.gro -o index.ndx

# Common commands inside make_ndx:
#   r 41 145     → select residues 41 and 145
#   name 6 ActiveSite → rename group 6
#   q            → quit and save
```

---

## 8. Checkpoint Files

### 8.1 `.cpt` — Checkpoint File

**Origin:** Written periodically by `gmx mdrun` (default: every 15 minutes).

The checkpoint file contains the **complete simulation state**: coordinates, velocities, Nosé-Hoover / Parrinello-Rahman coupling variables, random number generator states, and the step number. It enables exact restart:

```bash
# Restart a simulation from a checkpoint
gmx mdrun -s md.tpr -cpi md.cpt -deffnm md
```

---

## 9. Log Files

### 9.1 `.log` — Simulation Log

**Origin:** Written by `gmx mdrun`.

The `.log` file contains:
- All simulation parameters as interpreted by GROMACS.
- Warnings and notes from `grompp` (also echoed here).
- Step-by-step energy output (if `nstlog` > 0).
- **Performance data** at the end (ns/day, core-hours, load balancing).

**Checking performance:**

```bash
# Jump to the performance section at the end
tail -30 md.log
```

Typical output:

```
               Core t (s)   Wall t (s)        (%)
       Time:     3456.000      432.000      800.0
                 (ns/day)    (hour/ns)
Performance:       20.00        1.200
```

---

# Part B — Force Field Models

---

## 10. What Is a Force Field?

A **force field** is a mathematical model that describes the potential energy of a molecular system as a function of atomic coordinates. It consists of:

1. **A functional form** — the equations that calculate energies and forces.
2. **A parameter set** — the numerical constants (bond lengths, angles, partial charges, Lennard-Jones parameters, etc.) that are fitted to experimental data and/or quantum mechanical calculations.

The force field is an *approximation*. It treats atoms as classical particles (no quantum effects), represents electrons implicitly through partial charges and van der Waals parameters, and uses simple harmonic or periodic functions for bonded interactions. Despite these simplifications, modern force fields are remarkably accurate for many biomolecular applications.

---

## 11. The Functional Form

Most biomolecular force fields share a common functional form:

$$V_{\text{total}} = V_{\text{bonded}} + V_{\text{nonbonded}}$$

### Bonded Interactions

$$V_{\text{bonded}} = \underbrace{\sum_{\text{bonds}} \frac{1}{2} k_b (r - r_0)^2}_{\text{Bond stretching}} + \underbrace{\sum_{\text{angles}} \frac{1}{2} k_\theta (\theta - \theta_0)^2}_{\text{Angle bending}} + \underbrace{\sum_{\text{dihedrals}} k_\phi [1 + \cos(n\phi - \delta)]}_{\text{Torsional rotation}} + \underbrace{\sum_{\text{impropers}} k_\xi (\xi - \xi_0)^2}_{\text{Out-of-plane}}$$

| Term | What It Models | Key Parameters |
|---|---|---|
| **Bonds** | Covalent bond stretching | $k_b$ (stiffness), $r_0$ (equilibrium length) |
| **Angles** | Bond angle bending | $k_\theta$ (stiffness), $\theta_0$ (equilibrium angle) |
| **Proper dihedrals** | Rotation around bonds | $k_\phi$ (barrier height), $n$ (periodicity), $\delta$ (phase) |
| **Improper dihedrals** | Planarity of sp² centres | $k_\xi$ (stiffness), $\xi_0$ (equilibrium angle) |

### Nonbonded Interactions

$$V_{\text{nonbonded}} = \underbrace{\sum_{i<j} 4\varepsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6 \right]}_{\text{Lennard-Jones (van der Waals)}} + \underbrace{\sum_{i<j} \frac{q_i q_j}{4\pi\varepsilon_0 r_{ij}}}_{\text{Coulomb (electrostatics)}}$$

| Term | What It Models | Key Parameters |
|---|---|---|
| **Lennard-Jones** | Van der Waals attraction + Pauli repulsion | $\varepsilon$ (well depth), $\sigma$ (contact distance) |
| **Coulomb** | Electrostatic interactions | $q_i$, $q_j$ (partial charges) |

> **1-4 interactions:** Atoms separated by exactly three bonds receive scaled nonbonded interactions. The scale factors differ by force field (e.g. AMBER: 1/2 for LJ, 1/1.2 for Coulomb; CHARMM: 1.0 for both with explicit 1-4 LJ parameters).

---

## 12. Force Field Directory Layout in GROMACS

GROMACS stores force fields as directories (typically in `$GMXLIB` or included with the GROMACS installation):

```
amber99sb-ildn.ff/
├── forcefield.itp          # Master include: defines defaults, atom types, nonbonded params
├── forcefield.doc          # Description of the force field
├── aminoacids.rtp          # Residue Topology Parameters: all standard amino acids
├── aminoacids.hdb          # Hydrogen database: rules for adding H atoms
├── aminoacids.r2b          # Residue-to-building-block mapping
├── aminoacids.arn          # Residue name aliases
├── aminoacids.n.tdb        # N-terminal modifications (patches)
├── aminoacids.c.tdb        # C-terminal modifications (patches)
├── ffbonded.itp            # All bonded parameters (bonds, angles, dihedrals)
├── ffnonbonded.itp         # All nonbonded parameters (atom types, LJ, charges)
├── tip3p.itp               # Water model topology
├── ions.itp                # Monatomic ion parameters
├── spc.itp                 # Alternative water model
├── gb.itp                  # Implicit solvent parameters (if supported)
└── watermodels.dat         # List of available water models for this FF
```

Key files:

| File | Purpose |
|---|---|
| `forcefield.itp` | Entry point: sets `[ defaults ]` (combination rules, fudge factors), `#include`s `ffnonbonded.itp` and `ffbonded.itp` |
| `ffnonbonded.itp` | `[ atomtypes ]` — atom type names, masses, charges, LJ σ/ε |
| `ffbonded.itp` | `[ bondtypes ]`, `[ angletypes ]`, `[ dihedraltypes ]` — equilibrium values and force constants |
| `aminoacids.rtp` | Templates for building each residue (atom names, charges, bonds, impropers) |
| `aminoacids.hdb` | Rules for protonation (which hydrogens to add and where) |
| `*.tdb` | Terminal patches (NH₃⁺ at N-terminus, COO⁻ at C-terminus, etc.) |

---

## 13. All-Atom Force Fields

### 13.1 AMBER Family (ff99SB, ff14SB, ff19SB)

**Origin:** Developed at the University of California (Peter Kollman's group and successors). "AMBER" stands for Assisted Model Building with Energy Refinement. The name refers to both the force field and the separate AMBER software suite, but the parameters are available in GROMACS format.

**Lineage:** ff94 → ff99 → ff99SB → ff99SB-ILDN → ff14SB → ff19SB

| Aspect | Details |
|---|---|
| **Current recommended version** | **ff19SB** (2020) |
| **Parameterisation** | Fitted to quantum mechanical gas-phase data + condensed-phase thermodynamic properties |
| **Combination rules** | Lorentz-Berthelot |
| **1-4 scaling** | LJ: 0.5, Coulomb: 1/1.2 (= 0.8333) |
| **Charge model** | RESP (Restrained Electrostatic Potential) charges from QM |
| **GROMACS port** | `amber99sb-ildn.ff` (ff99SB-ILDN built-in); ff14SB and ff19SB available via community ports |

**Advantages:**
- Excellent backbone dihedral parameters, particularly for disordered proteins and IDPs (ff19SB).
- Very well validated for protein folding and dynamics.
- Large user community; extensive literature for validation.
- RESP charges are physically grounded and well-tested.

**Disadvantages:**
- Lipid and carbohydrate parameters require separate AMBER-family FFs (Lipid17, GLYCAM).
- Older versions (ff99SB) over-stabilise α-helices.
- Transferring parameters from AMBER software format to GROMACS can introduce subtle errors if not done carefully.

**Common use cases:**
- Protein folding and stability studies
- Intrinsically disordered proteins (ff19SB + OPC water)
- Protein-ligand binding (combined with GAFF/GAFF2 for ligands)
- Free energy perturbation calculations

---

### 13.2 CHARMM Family (CHARMM36 / CHARMM36m)

**Origin:** Developed at Harvard (Martin Karplus's group) and the University of Maryland (Alexander MacKerell's group). "CHARMM" = Chemistry at HARvard Macromolecular Mechanics.

| Aspect | Details |
|---|---|
| **Current recommended version** | **CHARMM36m** (2017) |
| **Parameterisation** | Iterative fitting to QM data + condensed-phase properties (osmotic pressure, NMR data, crystal structures) |
| **Combination rules** | Lorentz-Berthelot |
| **1-4 scaling** | LJ: 1.0 (but with explicit NBFIX pairs), Coulomb: 1.0 |
| **Charge model** | Group-based charges fitted to reproduce interaction energies with water |
| **GROMACS port** | `charmm36-jul2022.ff` (official port at mackerell.umaryland.edu/charmm_ff.shtml) |

**Advantages:**
- **Universally parameterised:** proteins, nucleic acids, lipids, carbohydrates, and drug-like molecules are all parametrised *together* in a consistent framework — ideal for complex systems.
- CHARMM36m improves backbone sampling for IDPs and disordered regions.
- CHARMM-GUI (charmm-gui.org) provides an automated system builder that generates GROMACS-ready files.
- NBFIX corrections for specific ion-protein and ion-nucleic acid interactions.
- CGenFF (CHARMM General Force Field) enables automatic parameterisation of small molecules.

**Disadvantages:**
- CHARMM format uses `NBFIX` (pair-specific LJ overrides) extensively — not all software handles this correctly.
- Slightly less accurate for protein folding in some benchmarks vs. ff19SB.
- Requires CHARMM-modified TIP3P water (differs from standard TIP3P: has LJ on hydrogens).

**Common use cases:**
- Membrane protein simulations (CHARMM36 lipids are the gold standard)
- Protein-membrane-ligand systems
- Nucleic acid simulations (DNA/RNA)
- Glycoprotein modelling (CHARMM36 carbohydrates)
- Drug design with CGenFF ligands

---

### 13.3 OPLS-AA Family (OPLS-AA, OPLS-AA/M)

**Origin:** Developed at Yale by William Jorgensen. "OPLS" = Optimized Potentials for Liquid Simulations. "AA" = All-Atom.

| Aspect | Details |
|---|---|
| **Current recommended version** | **OPLS-AA/M** (2015+, ongoing updates) |
| **Parameterisation** | Fitted primarily to reproduce **liquid-state thermodynamic properties** (densities, heats of vaporisation) |
| **Combination rules** | Geometric mean for both σ and ε |
| **1-4 scaling** | LJ: 0.5, Coulomb: 0.5 |
| **Charge model** | 1.14×CM1A-LBCC (derived from semiempirical QM + empirical corrections) |
| **GROMACS port** | `oplsaa.ff` (built-in); updated ports for OPLS-AA/M available from the Jorgensen group |

**Advantages:**
- Historically strong for small organic molecules and liquid simulations.
- Geometric combination rules make it straightforward to mix parameters.
- LigParGen server provides automated ligand parameterisation.
- OPLS-AA/M significantly improved protein backbone parameters.

**Disadvantages:**
- Lipid and nucleic acid parameters are less mature than AMBER or CHARMM equivalents.
- Smaller user community compared to AMBER/CHARMM for biomolecular work.
- The built-in GROMACS `oplsaa.ff` is outdated — users must manually install newer OPLS-AA/M ports.
- Fewer automated tools for complex system setup.

**Common use cases:**
- Organic liquid simulations
- Protein-ligand binding (particularly with organic solvents)
- Materials science applications
- Host-guest chemistry
- Small-molecule solvation free energies

---

### 13.4 GROMOS Family (54A7, 54A8)

**Origin:** Developed at ETH Zürich / University of Groningen (Wilfred van Gunsteren's group). "GROMOS" = GROningen MOlecular Simulation. Note: GROMOS is a **united-atom** force field — nonpolar hydrogens (e.g. on CH₃, CH₂) are merged into the heavy atom.

| Aspect | Details |
|---|---|
| **Current recommended versions** | **54A7** (2011), **54A8** (2013) |
| **Parameterisation** | Fitted to reproduce condensed-phase thermodynamic properties (free energies of solvation and hydration) |
| **Combination rules** | Geometric mean for C₆ and C₁₂ parameters (not σ/ε) |
| **Representation** | United-atom (nonpolar H atoms are implicit) |
| **GROMACS port** | `gromos54a7.ff` (built-in) |

**Advantages:**
- **Speed:** ~3–4× fewer atoms than all-atom force fields (no nonpolar hydrogens) → significantly faster simulations.
- Well-validated free energy calculations (thermodynamic integration).
- Good performance for lipid bilayers (in united-atom representation).
- C₆/C₁₂ parameterisation directly in Lennard-Jones form — avoids σ/ε conversion issues.

**Disadvantages:**
- **United-atom approximation** loses detail at the hydrogen level — less accurate for hydrogen bonding geometries and NMR observables.
- Less actively developed than AMBER/CHARMM (development has slowed).
- Force field uses a different parameter format (C₆/C₁₂) which can cause confusion.
- Fewer community resources and tutorials compared to AMBER/CHARMM.
- Limited ligand parameterisation tools compared to CGenFF or GAFF.

**Common use cases:**
- Large-scale membrane simulations (where speed is critical)
- Free energy calculations (alchemical transformations)
- Coarse "first-pass" exploration of protein dynamics
- Historical studies and method development

---

## 14. Coarse-Grained Force Fields

### 14.1 MARTINI (v2 / v3)

**Origin:** Developed at the University of Groningen (Siewert-Jan Marrink's group). Named after the Martini cocktail — the idea being to provide a "smooth" (simplified) representation.

MARTINI maps roughly **4 heavy atoms → 1 coarse-grained (CG) bead**, reducing the system size by ~10× and enabling simulations on the **microsecond–millisecond** timescale.

| Aspect | Details |
|---|---|
| **Current version** | **MARTINI 3** (2021) |
| **Resolution** | ~4:1 mapping (4 heavy atoms per bead) |
| **Parameterisation** | Top-down: fitted to reproduce partition free energies between water and organic solvents |
| **GROMACS support** | Fully supported; martinize2 tool converts atomistic structures to CG |
| **Time step** | 20–30 fs (vs. 2 fs for all-atom) |

**MARTINI bead types (v3):**

| Type | Description | Example Mapping |
|---|---|---|
| **P** (Polar) | Hydrophilic | Water (4 water molecules → 1 bead) |
| **N** (Nonpolar) | Intermediate | Backbone amide |
| **C** (Apolar) | Hydrophobic | Alkyl chains |
| **Q** (Charged) | Charged | Lysine NH₃⁺, glutamate COO⁻ |
| **S** (Small) | 2–3 atom mapping | Ring fragments |

**Advantages:**
- **Massive speed-up:** 100–1000× faster than all-atom simulations.
- Accesses microsecond–millisecond timescales for large systems (membranes, vesicles, protein aggregation).
- Good for lipid self-assembly, membrane properties, protein-membrane interactions.
- MARTINI 3 has improved small molecule representation and water/ion behaviour.
- Widely used and actively developed.

**Disadvantages:**
- **Cannot capture atomistic detail** — no hydrogen bonds, no explicit water structure, no side chain rotamers.
- Protein secondary structure must be restrained (MARTINI does not fold proteins).
- Dielectric constant and effective time must be calibrated (the "MARTINI time" factor of ~4× is approximate).
- Free energies and kinetics should be interpreted cautiously.
- Back-mapping to atomistic resolution introduces artefacts that require further equilibration.

**Common use cases:**
- Membrane assembly and lipid phase behaviour
- Protein-membrane interactions and insertion
- Large-scale protein aggregation (amyloid, phase separation)
- Vesicle formation and fusion
- Drug delivery system design (liposomes, nanoparticles)
- Screening: identify interesting configurations, then refine with all-atom MD

---

## 15. Specialised & Polarisable Force Fields

### 15.1 AMOEBA

**Origin:** Developed by Jay Ponder (Washington University in St. Louis) and Pengyu Ren.

AMOEBA (Atomic Multipole Optimized Energetics for Biomolecular Applications) is a **polarisable force field** that uses atomic multipoles (up to quadrupoles) instead of fixed point charges, plus induced dipoles to model electronic polarisation.

| Aspect | Details |
|---|---|
| **Electrostatics** | Permanent atomic multipoles + mutual induced dipoles |
| **Polarisation** | Self-consistent induced dipoles via Thole damping |
| **VdW** | Buffered 14-7 potential (softer than 12-6 LJ) |
| **GROMACS support** | Limited — primarily used in Tinker/OpenMM |

**Advantages:**
- Most physically accurate electrostatics of any biomolecular FF.
- Excellent for ion solvation, metal coordination, and polarisation-sensitive phenomena.
- No need for NBFIX-style corrections.

**Disadvantages:**
- **5–10× more expensive** than fixed-charge force fields.
- Requires iterative self-consistent field calculation each step.
- Limited GROMACS support; primarily used in Tinker and OpenMM.
- Smaller parameter coverage for biomolecules.

**Common use cases:** Ion channels, metalloenzymes, protein-ion interactions, QM/MM validation.

### 15.2 Drude Polarisable Force Field

**Origin:** Developed by the MacKerell group (same group as CHARMM) as a polarisable extension of CHARMM.

The Drude model attaches a massless charged "Drude particle" to each polarisable atom via a harmonic spring. The displacement of the Drude particle represents induced polarisation.

| Aspect | Details |
|---|---|
| **Electrostatics** | Fixed charges + Drude oscillator particles |
| **Polarisation** | Classical Drude oscillator (charge-on-spring) |
| **GROMACS support** | Supported since GROMACS 2021 |

**Advantages:**
- Compatible with the CHARMM ecosystem.
- Captures polarisation effects at lower cost than AMOEBA.
- Parameters available for proteins, nucleic acids, lipids, and ions.
- Supported in GROMACS (unlike AMOEBA).

**Disadvantages:**
- ~2–3× more expensive than fixed-charge CHARMM36.
- Requires a shorter time step (~1 fs) or extended Lagrangian integration.
- Fewer validation studies than CHARMM36m.
- More complex system setup and equilibration protocols.

**Common use cases:** Ion transport, membrane electrostatics, systems where polarisation is critical.

---

## 16. Water Models

The choice of water model matters enormously — water typically makes up >80% of the atoms in a solvated simulation. Each force field is parameterised against a specific water model.

| Model | Sites | Force Field Compatibility | Key Properties |
|---|---|---|---|
| **TIP3P** | 3 | AMBER (standard) | Fast; underestimates viscosity, overestimates diffusion |
| **CHARMM-TIP3P** | 3 (modified) | CHARMM | Like TIP3P but with LJ on H atoms — **must** use with CHARMM |
| **SPC** | 3 | GROMOS | Simple; similar trade-offs to TIP3P |
| **SPC/E** | 3 | General | Better density/diffusion than SPC; includes self-polarisation correction |
| **TIP4P** | 4 | General | Off-atom charge site; better dielectric properties |
| **TIP4P-Ew** | 4 | AMBER (improved) | Reparameterised TIP4P for Ewald summation |
| **TIP4P/2005** | 4 | General | Best overall 4-site model; good phase diagram |
| **OPC** | 4 | AMBER (ff19SB) | Optimal Point Charge; best with ff19SB for IDPs |
| **TIP5P** | 5 | General | Two lone pairs; better tetrahedral structure, but slow |

> **Rule of thumb:** Always use the water model that your force field was parameterised with. Mixing (e.g. TIP4P with AMBER ff99SB) can give incorrect results.

---

## 17. Force Field Comparison Summary

| Force Field | Type | Speed | Proteins | Lipids | Nucleic Acids | Ligands | Best For |
|---|---|---|---|---|---|---|---|
| **AMBER ff19SB** | All-atom | ●●●○ | ★★★★★ | ★★★☆☆ ¹ | ★★★★☆ | ★★★★☆ (GAFF2) | Protein folding, IDPs |
| **CHARMM36m** | All-atom | ●●●○ | ★★★★☆ | ★★★★★ | ★★★★★ | ★★★★☆ (CGenFF) | Membranes, multi-component systems |
| **OPLS-AA/M** | All-atom | ●●●○ | ★★★★☆ | ★★☆☆☆ | ★★☆☆☆ | ★★★★☆ (LigParGen) | Organic liquids, small molecules |
| **GROMOS 54A7** | United-atom | ●●●●○ | ★★★☆☆ | ★★★☆☆ | ★★☆☆☆ | ★★★☆☆ (ATB) | Speed-critical, free energies |
| **MARTINI 3** | Coarse-grained | ●●●●● | ★★★☆☆ ² | ★★★★★ | ★★☆☆☆ | ★★★☆☆ | Large-scale, long timescales |
| **Drude** | Polarisable | ●●○○○ | ★★★★☆ | ★★★☆☆ | ★★★☆☆ | ★★☆☆☆ | Polarisation-sensitive systems |
| **AMOEBA** | Polarisable | ●○○○○ | ★★★★☆ | ★☆☆☆☆ | ★★☆☆☆ | ★★★☆☆ | Metal coordination, ions |

¹ Lipid17 is available as a separate AMBER-family FF.
² MARTINI proteins require secondary structure restraints; cannot fold.

Speed: ● = slowest, ●●●●● = fastest. Quality: ★ = poorest, ★★★★★ = best.

---

## 18. How to Choose a Force Field

Use this decision tree as a starting point:

```
Is your system primarily a membrane/lipid system?
├── Yes → CHARMM36m (or MARTINI for very large systems)
└── No
    ├── Is it a protein in solution?
    │   ├── Is the protein intrinsically disordered?
    │   │   ├── Yes → AMBER ff19SB + OPC water
    │   │   └── No → AMBER ff14SB/ff19SB or CHARMM36m (both excellent)
    │   └── Does it contain a ligand?
    │       ├── Yes → CHARMM36m + CGenFF  OR  AMBER ff19SB + GAFF2
    │       └── No → Any of the above
    ├── Is it a nucleic acid system?
    │   └── CHARMM36 (DNA/RNA well parameterised)
    ├── Do you need microsecond+ timescales or very large (>1M atom) systems?
    │   └── MARTINI 3 → then back-map and refine with all-atom
    ├── Is electronic polarisation critical?
    │   └── Drude (if GROMACS) or AMOEBA (if OpenMM/Tinker)
    └── Is it an organic liquid / materials science problem?
        └── OPLS-AA/M
```

**General rules:**
1. **Consistency is king.** Never mix parameters from different force fields.
2. **Use the recommended water model** for your chosen force field.
3. **Check the literature** for your system — has someone validated a particular FF for your protein/membrane/ligand?
4. **When in doubt, CHARMM36m or AMBER ff19SB** — they have the broadest validation and largest user bases.

---

## 19. Summary & Next Steps

### What We Covered

| Topic | Key Takeaway |
|---|---|
| **`.gro`** | GROMACS coordinate file; fixed-width, uses nanometres |
| **`.top` / `.itp`** | Topology hierarchy: `.top` includes `.itp` files; `[ molecules ]` must match coordinate file |
| **`.mdp`** | User-written simulation settings — controls every aspect of the run |
| **`.tpr`** | Binary run input — self-contained, portable |
| **`.xtc` / `.trr`** | Compressed vs. full-precision trajectories |
| **`.edr` / `.log`** | Energy data and human-readable simulation log |
| **`.ndx`** | Custom atom group definitions |
| **`.cpt`** | Checkpoint for exact restart |
| **Force fields** | Functional form is shared; parameter sets differ in philosophy and accuracy |
| **AMBER/CHARMM/OPLS/GROMOS** | All-atom (or united-atom) options with different strengths |
| **MARTINI** | Coarse-grained — fast but low-resolution |
| **Polarisable FFs** | More accurate electrostatics at higher computational cost |
| **Water models** | Must match your force field; never mix |

### Next Steps

- **Preparing a system for MD** → See the main [MD on Apocrita](../MD_On_Apocrita.md) tutorial, which walks through `pdb2gmx`, solvation, ion addition, energy minimisation, and equilibration.
- **Understanding protein structure files** → See [Understanding Protein Files](Understanding_Protein_Files.md) for a deep dive into PDB and mmCIF formats with Python parsing.
- **Writing `.mdp` files** → The [GROMACS manual](https://manual.gromacs.org/current/user-guide/mdp-options.html) has comprehensive documentation for every `.mdp` parameter.
- **Parameterising ligands** → Look into [CGenFF](https://cgenff.umaryland.edu/) (for CHARMM), [ACPYPE](https://github.com/alanwilter/acpype) / [antechamber](https://ambermd.org/antechamber/antechamber.html) (for AMBER/GAFF), or [ATB](https://atb.uq.edu.au/) (for GROMOS).

---

*Reference sheet written for the QMUL Molecular Dynamics Teaching Project. Last updated: March 2026.*
