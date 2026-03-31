# GROMACS Key Definitions

**A Comprehensive Glossary for Students of Molecular Dynamics Simulation**

Use this reference when you encounter an unfamiliar term in GROMACS documentation, tutorials, `.mdp` files, log output, or the scientific literature. Entries are grouped by topic and ordered alphabetically within each section.

> **Cross-references:** For deeper explanations of many of these concepts, see the companion tutorials:
> - [Understanding GROMACS Theory](Understanding_GROMACS_Theory.md)
> - [Understanding GROMACS Key Files & Force Field Models](Understanding_GROMACS_Key_Files.md)
> - [GROMACS Common Pitfalls & Mastery](GROMACS_Common_Pitfalls_And_Mastery.md)
> - [Understanding Protein Files](Understanding_Protein_Files.md)

---

## Table of Contents

1. [Units & Fundamental Constants](#1-units--fundamental-constants)
2. [Structural Biology & Protein Terminology](#2-structural-biology--protein-terminology)
3. [Force Field & Energy Terminology](#3-force-field--energy-terminology)
4. [Integration & Dynamics](#4-integration--dynamics)
5. [Constraints & Rigid Bodies](#5-constraints--rigid-bodies)
6. [Electrostatics & Long-Range Interactions](#6-electrostatics--long-range-interactions)
7. [Periodic Boundary Conditions & Box Geometry](#7-periodic-boundary-conditions--box-geometry)
8. [Temperature Control (Thermostats)](#8-temperature-control-thermostats)
9. [Pressure Control (Barostats)](#9-pressure-control-barostats)
10. [Statistical Mechanics & Thermodynamics](#10-statistical-mechanics--thermodynamics)
11. [Energy Minimisation](#11-energy-minimisation)
12. [Solvation, Ions & System Setup](#12-solvation-ions--system-setup)
13. [Trajectory & Analysis Terminology](#13-trajectory--analysis-terminology)
14. [Free Energy Methods](#14-free-energy-methods)
15. [Enhanced Sampling Methods](#15-enhanced-sampling-methods)
16. [Coarse-Graining](#16-coarse-graining)
17. [GROMACS-Specific Commands & Tools](#17-gromacs-specific-commands--tools)
18. [File Formats](#18-file-formats)
19. [HPC & Performance Terminology](#19-hpc--performance-terminology)
20. [Water Models](#20-water-models)
21. [Common Abbreviations](#21-common-abbreviations)

---

## 1. Units & Fundamental Constants

### Ångström (Å)
Unit of length equal to $10^{-10}$ m (0.1 nm). The traditional unit in crystallography and PDB files. **GROMACS does not use Ångströms internally — it uses nanometres.** 1 nm = 10 Å.

### Bar
Unit of pressure. 1 bar ≈ 0.987 atm ≈ 100 kPa. Standard simulation pressure is 1 bar.

### Boltzmann Constant ($k_B$)
$k_B = 1.38065 \times 10^{-23}$ J/K = 8.3145 × 10⁻³ kJ/(mol·K). Relates temperature to energy at the molecular level. The thermal energy at 300 K is $k_BT \approx 2.494$ kJ/mol.

### Dalton (Da)
Unit of atomic mass. 1 Da = 1 atomic mass unit (amu) = 1.66054 × 10⁻²⁷ kg. Carbon-12 has a mass of exactly 12 Da. GROMACS uses amu internally.

### Elementary Charge ($e$)
The charge of a proton: $e = 1.602 \times 10^{-19}$ C. Partial charges in force fields are expressed in units of $e$ (e.g. the oxygen of TIP3P water has charge −0.834 $e$).

### Femtosecond (fs)
$10^{-15}$ s. Bond vibrations occur on the ~10 fs timescale.

### Kelvin (K)
SI unit of absolute temperature. Room temperature ≈ 298–300 K. Physiological temperature = 310 K (37 °C).

### kJ/mol (Kilojoule per mole)
GROMACS's standard energy unit. To convert: 1 kcal/mol = 4.184 kJ/mol. Thermal energy $k_BT$ at 300 K ≈ 2.494 kJ/mol.

### Nanometre (nm)
$10^{-9}$ m. **The primary length unit in GROMACS.** All coordinates, cutoffs, and box dimensions are in nm. 1 nm = 10 Å.

### Nanosecond (ns)
$10^{-9}$ s. A typical "short" MD simulation runs for 10–100 ns. Microsecond ($\mu$s = $10^{-6}$ s) simulations are routine on modern GPUs.

### Picosecond (ps)
$10^{-12}$ s. **GROMACS's standard time unit.** The time step `dt = 0.002` means 2 fs = 0.002 ps. Output settings like `nstxout-compressed = 5000` with `dt = 0.002` write every 10 ps.

---

## 2. Structural Biology & Protein Terminology

### Alpha-Carbon (Cα / CA)
The central carbon atom of each amino acid, bonded to the amino group (N), carboxyl group (C), hydrogen (HA), and side chain (CB). Cα atoms define the backbone trace and are commonly used for RMSD and superposition.

### Alpha-Helix (α-helix)
A common secondary structure element: a right-handed helical coil stabilised by i → i+4 backbone hydrogen bonds (CO of residue $i$ to NH of residue $i+4$). ~3.6 residues per turn, rise of 0.54 nm per turn.

### Amino Acid
The monomer unit of proteins. 20 standard amino acids, each with a backbone (N–Cα–C=O) and a variable side chain (R group). Represented by three-letter codes (e.g. ALA, GLY, LEU) or one-letter codes (A, G, L).

### B-factor (Temperature Factor)
A measure of atomic displacement / thermal motion, reported in Å². High B-factors indicate flexible or disordered regions. Stored in the `occupancy` and `B-factor` columns of PDB `ATOM` records. Formally, $B = 8\pi^2 \langle u^2 \rangle$ where $\langle u^2 \rangle$ is the mean-square atomic displacement.

### Backbone
The repeating N–Cα–C(=O) chain of a polypeptide. Excludes side chain atoms. The backbone dihedral angles φ and ψ determine the protein's secondary structure (see **Ramachandran plot**).

### Beta-Sheet (β-sheet)
A secondary structure formed by extended polypeptide strands connected by hydrogen bonds between backbone NH and CO groups of adjacent strands. Can be parallel (same strand direction) or antiparallel (opposite direction).

### Chain
A single covalently connected polypeptide (or nucleic acid strand) within a PDB structure, identified by a chain ID (A, B, C, ...).

### Cofactor
A non-protein molecule required for a protein's biological function. Can be organic (coenzyme, e.g. NAD⁺, FAD) or inorganic (metal ion, e.g. Zn²⁺, Fe²⁺/³⁺). Cofactors appear as `HETATM` records in PDB files and usually require custom force field parameterisation.

### Conformation
A specific 3D arrangement of atoms in a molecule. Proteins can adopt many conformations; the native state is the biologically active one.

### Cryo-EM (Cryo-Electron Microscopy)
An experimental technique that determines protein structure by imaging frozen particles with an electron beam. Resolution typically 2–4 Å. Produces electron density maps from which atomic models are built.

### Disulfide Bond (S–S Bond)
A covalent bond between the sulfur atoms (Sγ) of two cysteine residues. Stabilises protein structure. Must be explicitly defined in the GROMACS topology. PDB files list these as `SSBOND` records.

### Dihedral Angle (Torsion)
The angle between two planes defined by four consecutive atoms A–B–C–D. The backbone dihedrals **phi (φ)** (C–N–Cα–C) and **psi (ψ)** (N–Cα–C–N) define secondary structure preferences.

### Heterogen / HETATM
An atom belonging to a non-standard residue. In PDB files, these are recorded as `HETATM` lines (as opposed to `ATOM` lines for standard residues). Includes ligands, cofactors, modified residues, and water molecules (HOH).

### Hydrogen Bond (H-bond)
A non-covalent interaction between a hydrogen atom bonded to an electronegative atom (donor, e.g. N–H) and another electronegative atom (acceptor, e.g. O, N). Typical energy: 5–30 kJ/mol. Key stabiliser of secondary structure and protein–water interactions. GROMACS default criteria: donor–acceptor distance ≤ 0.35 nm, angle ≤ 30°.

### Intrinsically Disordered Protein (IDP)
A protein (or protein region) that lacks a stable 3D structure under physiological conditions. IDPs are functional and involved in signalling, regulation, and binding. They are notoriously difficult to simulate — force fields often make them too compact.

### Ligand
A molecule that binds to a protein, typically at a specific binding site. Can be a drug, substrate, inhibitor, or signalling molecule. In GROMACS, ligands require custom parameterisation (e.g. via CGenFF, GAFF, ACPYPE).

### Loop
A stretch of polypeptide connecting regular secondary structure elements (helices and strands). Loops are often flexible, have high B-factors, and are frequently missing in crystal structures.

### N-terminus / C-terminus
The two ends of a polypeptide chain. The N-terminus has a free amino group (typically protonated to NH₃⁺ at physiological pH); the C-terminus has a free carboxyl group (deprotonated to COO⁻). These terminal groups are "patched" by GROMACS's `pdb2gmx` using rules in the `.tdb` files.

### NMR (Nuclear Magnetic Resonance)
An experimental technique that determines protein structure in solution by measuring distances between atomic nuclei. Produces an **ensemble** of structures (typically 10–20 models). NMR structures in PDB files contain multiple `MODEL` records.

### Occupancy
In a PDB file, the fraction of time an atom is observed at a given position (0 to 1). An occupancy of 1.00 means the atom is always in this position. Values < 1 indicate alternate conformations (see **Alternate Conformation**).

### PDB (Protein Data Bank)
The worldwide repository for experimentally determined 3D structures of biological macromolecules ([rcsb.org](https://www.rcsb.org)). Each entry has a 4-character ID (e.g. `1UBQ`, `6LU7`).

### PDB ID
A 4-character alphanumeric code uniquely identifying a structure in the PDB (e.g. `1UBQ` for ubiquitin).

### Phi (φ) / Psi (ψ)
The two main backbone dihedral angles. φ is defined by atoms C(i-1)–N(i)–Cα(i)–C(i). ψ is defined by N(i)–Cα(i)–C(i)–N(i+1). Together they determine the secondary structure of each residue and are plotted on the Ramachandran plot.

### Protonation State
The number of hydrogen atoms (protons) on a titratable group at a given pH. Histidine, aspartate, glutamate, cysteine, lysine, and tyrosine all have pH-dependent protonation states. Must be assigned before simulation.

### Ramachandran Plot
A scatter plot of backbone dihedral angles φ vs. ψ for all residues. Allowed regions correspond to common secondary structures: α-helix (φ ≈ −60°, ψ ≈ −45°), β-sheet (φ ≈ −120°, ψ ≈ +130°). Used as a structure validation tool.

### Residue
A single amino acid unit within a polypeptide chain. Identified by a residue name (e.g. ALA), residue number (e.g. 48), and chain ID (e.g. A). In GROMACS topology, residues are defined in `.rtp` (Residue Topology Parameter) files.

### Resolution
A measure of the detail visible in an experimentally determined structure, in Ångströms. Lower is better: sub-1.0 Å = atomic resolution, 1.5–2.5 Å = typical high-resolution, 3.0–4.0 Å = medium, >4.0 Å = low. Only applies to diffraction-based and cryo-EM methods, not NMR.

### Salt Bridge
An electrostatic interaction between a positively charged residue (Arg, Lys, His⁺) and a negatively charged residue (Asp, Glu) at close distance (<0.4 nm). Important for protein stability and function.

### Secondary Structure
The local 3D arrangement of the polypeptide backbone, classified into α-helices, β-sheets, turns, coils, etc. Assigned from coordinates using algorithms like DSSP or STRIDE.

### Side Chain
The variable "R group" that distinguishes the 20 amino acids. Attached to the Cα atom. Determines the chemical properties of each residue (hydrophobic, polar, charged, aromatic, etc.).

### Solvent-Accessible Surface Area (SASA)
The surface area of a protein that is accessible to a solvent molecule (modelled as a sphere, typically 0.14 nm radius). Computed by `gmx sasa`. Changes in SASA can indicate folding, binding, or conformational change.

### Turn
A short secondary structure element (3–5 residues) that reverses the direction of the polypeptide chain.

### X-ray Crystallography
The dominant experimental technique for protein structure determination. A protein crystal is exposed to X-rays, producing a diffraction pattern from which electron density is calculated. Individual atom positions are fitted into the density. Produces a single model (not an ensemble).

---

## 3. Force Field & Energy Terminology

### Angle Bending
A bonded interaction between three atoms (A–B–C) that resists deviation from the equilibrium bond angle $\theta_0$. Modelled as a harmonic potential: $V(\theta) = \frac{1}{2}k_\theta(\theta - \theta_0)^2$.

### AMBER
A family of force fields developed at the University of California (Assisted Model Building with Energy Refinement). Key versions: ff99SB, ff14SB, ff19SB. Uses RESP charges, Lorentz-Berthelot combining rules. Distributed both as the AMBER software and as GROMACS-compatible parameter files.

### Atom Type
A label that encodes the chemical nature and hybridisation of an atom within a force field (e.g. `CT` = sp³ carbon, `N3` = protonated nitrogen in AMBER). Atom types determine LJ parameters and which bonded parameters apply. Defined in `ffnonbonded.itp`.

### B-spline Interpolation
A mathematical method used in PME to distribute point charges onto a regular grid. The `pme_order` parameter (typically 4 = cubic B-spline) controls the order of interpolation.

### Bond Stretching
The simplest bonded interaction: resistance to changing a covalent bond length. Modelled as $V(r) = \frac{1}{2}k_b(r - r_0)^2$ (harmonic spring). The force constant $k_b$ is very large (200,000–500,000 kJ/mol/nm²) because covalent bonds are stiff.

### CHARMM
Chemistry at HARvard Macromolecular Mechanics. A force field family with strong coverage of proteins, lipids, nucleic acids, and carbohydrates. Current version: CHARMM36m. Uses Lorentz-Berthelot combining rules and NBFIX pair-specific corrections.

### Charge Group
A group of atoms whose partial charges sum to an integer (or near-integer). Used in older GROMACS cutoff schemes to reduce artefacts at the cutoff boundary. The Verlet cutoff scheme (standard since GROMACS 5.0) uses atom-based cutoffs and makes charge groups obsolete.

### Combination Rules (Combining Rules)
Rules for computing Lennard-Jones parameters between unlike atom types from individual atom parameters. **Lorentz-Berthelot:** $\sigma_{ij} = (\sigma_i + \sigma_j)/2$, $\varepsilon_{ij} = \sqrt{\varepsilon_i \varepsilon_j}$. **Geometric:** both are geometric means. Set by `comb-rule` in `[ defaults ]`.

### Coulomb Interaction
The electrostatic interaction between two partial charges: $V = q_i q_j / (4\pi\varepsilon_0 r)$. Decays as $1/r$ — very long-ranged, requiring special treatment (see **PME**).

### Dihedral (Torsion)
A bonded interaction between four atoms (A–B–C–D). **Proper dihedrals** model rotation around the B–C bond: $V(\phi) = k_\phi[1 + \cos(n\phi - \delta)]$. **Improper dihedrals** maintain planarity or chirality: $V(\xi) = \frac{1}{2}k_\xi(\xi - \xi_0)^2$.

### Dispersion
Attractive van der Waals interactions arising from correlated fluctuations in electron clouds (London dispersion forces). Scales as $r^{-6}$. Captured by the attractive term of the Lennard-Jones potential.

### Dispersion Correction (DispCorr)
An analytical correction for the Lennard-Jones energy and pressure beyond the cutoff distance. Applied when `DispCorr = EnerPres` in the `.mdp` file. Assumes a uniform atom distribution beyond the cutoff.

### Electrostatic Potential (ESP)
The electric potential at a point in space due to the charge distribution of a molecule. Force field partial charges (especially RESP charges in AMBER) are fitted to reproduce the quantum mechanical ESP around a molecule.

### Force Constant
The stiffness parameter in a harmonic potential. Higher force constants = stiffer interaction. Units: kJ mol⁻¹ nm⁻² (bonds), kJ mol⁻¹ rad⁻² (angles).

### Force Field
A complete mathematical model for molecular potential energy, consisting of a functional form (equations) and a parameter set (numerical values). In GROMACS, a force field is stored as a directory (e.g. `amber99sb-ildn.ff/`) containing `.itp`, `.rtp`, `.hdb`, and other files.

### Fudge Factors (fudgeLJ, fudgeQQ)
Scale factors applied to 1-4 nonbonded interactions. Defined in `[ defaults ]` of the force field. AMBER: fudgeLJ=0.5, fudgeQQ=0.8333. OPLS: fudgeLJ=0.5, fudgeQQ=0.5. CHARMM: fudgeLJ=1.0, fudgeQQ=1.0.

### GROMOS
GROningen MOlecular Simulation. A united-atom force field family (nonpolar hydrogens are implicit). Current versions: 54A7, 54A8. Uses C₆/C₁₂ parameters directly rather than σ/ε.

### Harmonic Potential
$V(x) = \frac{1}{2}k(x - x_0)^2$. Used for bonds, angles, improper dihedrals, and position restraints. The force is $F = -k(x - x_0)$ — linear restoring force, like a spring (Hooke's law).

### Improper Dihedral
A dihedral angle used to maintain planarity (e.g. peptide bond, aromatic ring) or chirality (e.g. Cα centre). Not a true torsion around a bond — instead defined by four atoms arranged to detect out-of-plane distortion.

### Lennard-Jones (LJ) Potential
$V(r) = 4\varepsilon[(\sigma/r)^{12} - (\sigma/r)^6]$. Models van der Waals interactions: the $r^{-12}$ term is short-range repulsion (Pauli exclusion); the $r^{-6}$ term is long-range attraction (London dispersion). Parameters: $\varepsilon$ (well depth, kJ/mol), $\sigma$ (contact distance, nm).

### Morse Potential
$V(r) = D_e[1 - e^{-a(r-r_0)}]^2$. A more physically realistic bond potential that correctly models bond dissociation. Rarely used in biomolecular force fields because it is more expensive and bond breaking is not the goal of classical MD.

### NBFIX
A CHARMM-specific mechanism for overriding combining rules with pair-specific Lennard-Jones parameters. Used extensively for ion–protein and ion–nucleic acid interactions where standard combining rules give incorrect results.

### 1-4 Interactions (One-Four Interactions)
Nonbonded interactions between atoms separated by exactly three bonds. These are handled specially (with scaled fudge factors) rather than using full nonbonded parameters, because the intervening bonded terms already partially account for the interaction.

### nrexcl (Number of Excluded Neighbours)
The `nrexcl` parameter in `[ moleculetype ]` specifies how many bonded neighbours are excluded from nonbonded calculations. `nrexcl = 3` (standard for proteins) means 1-2 and 1-3 nonbonded interactions are fully excluded; 1-4 interactions use fudge factors.

### OPLS-AA
Optimized Potentials for Liquid Simulations, All-Atom. A force field developed by William Jorgensen (Yale). Uses geometric combining rules for both σ and ε. Current version: OPLS-AA/M.

### Parameterisation
The process of determining numerical values for force field parameters (charges, LJ parameters, force constants) by fitting to experimental data and/or quantum mechanical calculations.

### Partial Charge
The effective fractional electric charge assigned to an atom in a force field. Not a physical observable — it is a model parameter designed to reproduce the molecular electrostatic potential. Always in units of elementary charge $e$. Must sum to the total molecular charge.

### Pauli Repulsion
The sharp increase in energy when electron clouds overlap, arising from the Pauli exclusion principle of quantum mechanics. Approximated by the $r^{-12}$ term in the Lennard-Jones potential.

### Potential Energy Surface (PES)
The hypersurface of energy as a function of all atomic coordinates: $V(\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_N)$. Minima are stable conformations; saddle points are transition states. The PES has an astronomically large number of local minima for a protein.

### RESP Charges (Restrained Electrostatic Potential)
The charge-fitting method used in AMBER force fields. Partial charges are determined by fitting to the quantum mechanical electrostatic potential while applying a restraint that keeps chemically equivalent atoms similar.

### Proper Dihedral
A torsional potential modelling rotation around a covalent bond. Uses a cosine functional form: $V(\phi) = k[1 + \cos(n\phi - \delta)]$, where $n$ is the multiplicity (number of minima per full rotation) and $\delta$ is the phase.

### United Atom
A representation where nonpolar hydrogen atoms (e.g. on CH₃, CH₂) are merged into their parent heavy atom, reducing the number of interaction sites. Used by GROMOS and some older force fields. Trades accuracy for speed.

### Van der Waals (VdW) Interaction
The collective term for weak, short-range attractive forces between atoms/molecules, arising from fluctuating and permanent dipoles. Modelled by the Lennard-Jones potential in most force fields.

---

## 4. Integration & Dynamics

### Acceleration
$\mathbf{a}_i = \mathbf{F}_i / m_i$. The acceleration of atom $i$ is the force divided by its mass (Newton's second law).

### Degrees of Freedom ($N_f$)
The number of independent variables needed to describe the state of the system. For $N$ atoms in 3D with $N_c$ constraints and COM removal: $N_f = 3N - N_c - 3$. Temperature is computed from $N_f$.

### Deterministic
A simulation where the outcome is fully determined by the initial conditions and parameters. Nosé-Hoover and Verlet integration are deterministic. Stochastic thermostats (V-rescale, Langevin) include random elements.

### dt (Time Step)
The discrete time increment used in numerical integration. In GROMACS, specified in picoseconds in the `.mdp` file: `dt = 0.002` = 2 fs. Must resolve the fastest motion in the system.

### Equation of Motion
Newton's second law applied to each atom: $m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i = -\nabla_i V$. The equations of motion are integrated numerically using the leap-frog algorithm.

### Hamiltonian
The total energy function of the system: $H = T + V$ (kinetic + potential). In a microcanonical (NVE) simulation, the Hamiltonian is conserved. With thermostats/barostats, an extended Hamiltonian is conserved.

### Hydrogen Mass Repartitioning (HMR)
A technique that increases hydrogen atom masses (to ~3–4 Da) by redistributing mass from their bonded heavy atoms. Slows down H-vibrations, allowing a 4 fs time step with `constraints = h-bonds`. The total mass of each residue is unchanged.

### Integrator
The numerical algorithm that advances positions and velocities through time. In GROMACS: `md` = leap-frog, `md-vv` = velocity Verlet, `steep` = steepest descent (minimisation), `cg` = conjugate gradient (minimisation), `sd` = stochastic/Langevin dynamics.

### Kinetic Energy
$E_{\text{kin}} = \frac{1}{2}\sum_i m_i v_i^2$. The energy of motion. Related to temperature via the equipartition theorem.

### Langevin Dynamics
A stochastic integrator that adds friction and random forces to Newton's equations: $m\ddot{\mathbf{r}} = \mathbf{F} - \gamma m \dot{\mathbf{r}} + \boldsymbol{\eta}(t)$. Acts as an implicit thermostat. Used via `integrator = sd` in GROMACS. Useful for implicit-solvent simulations or enhanced conformational sampling.

### Leap-Frog Algorithm
The default GROMACS integrator. Positions and velocities are evaluated at alternating half-steps: $\mathbf{v}(t+\frac{\Delta t}{2}) = \mathbf{v}(t-\frac{\Delta t}{2}) + \mathbf{a}(t)\Delta t$; $\mathbf{r}(t+\Delta t) = \mathbf{r}(t) + \mathbf{v}(t+\frac{\Delta t}{2})\Delta t$. Symplectic, time-reversible, requires one force evaluation per step.

### Maxwell-Boltzmann Distribution
The probability distribution of velocities at a given temperature: $P(v) \propto v^2 \exp(-mv^2/2k_BT)$. Initial velocities in MD are drawn from this distribution using `gen_vel = yes`.

### nsteps
The total number of integration steps. Total simulation time = `nsteps × dt`. E.g. `nsteps = 5000000`, `dt = 0.002` → 10,000 ps = 10 ns.

### Phase Space
The abstract multi-dimensional space whose coordinates are all positions and momenta of all atoms: $(\mathbf{r}_1, \ldots, \mathbf{r}_N, \mathbf{p}_1, \ldots, \mathbf{p}_N)$. Has $6N$ dimensions. Each point represents a complete microstate.

### Stochastic
Involving random elements. Stochastic thermostats (V-rescale, Langevin) add random perturbations to ensure correct canonical sampling. Stochastic dynamics `integrator = sd` includes friction and noise.

### Symplectic
A property of integrators that conserve phase-space volume (Liouville's theorem). Symplectic integrators (Verlet, leap-frog) exhibit excellent long-term energy conservation with no systematic drift. Non-symplectic methods (e.g. simple Euler) accumulate energy errors over time.

### Time-Reversibility
A property of Newtonian mechanics: reversing all velocities at time $t$ causes the system to retrace its trajectory. The leap-frog integrator preserves this property. Some thermostats (Berendsen, Langevin with friction) break time-reversibility.

### Trajectory
A time-ordered series of system snapshots (coordinates, and optionally velocities and forces). Stored in `.xtc` (compressed) or `.trr` (full-precision) files.

### Velocity Verlet
A variant of the Verlet integrator that evaluates positions and velocities at the same time $t$ (unlike leap-frog, which staggers them by $\Delta t/2$). Used via `integrator = md-vv`. Slightly more expensive per step but provides synchronised positions and velocities.

### Verlet Algorithm
The foundational integration scheme: $\mathbf{r}(t+\Delta t) = 2\mathbf{r}(t) - \mathbf{r}(t-\Delta t) + \mathbf{a}(t)\Delta t^2$. Position error $\mathcal{O}(\Delta t^4)$. Leap-frog and velocity Verlet are mathematically equivalent formulations.

---

## 5. Constraints & Rigid Bodies

### Constraint
A holonomic constraint fixes a degree of freedom (e.g. a bond length) to a constant value, removing it from the dynamics. Reduces the number of degrees of freedom $N_f$, enabling a larger time step.

### h-bonds (Hydrogen Bond Constraints)
`.mdp` setting: `constraints = h-bonds`. Constrains all bonds involving hydrogen atoms to their equilibrium length. This is the standard setting for 2 fs time steps with the leap-frog integrator.

### All-Bonds
`.mdp` setting: `constraints = all-bonds`. Constrains all bond lengths (not just X–H). Allows time steps of up to ~4–5 fs but removes some physical dynamics. Less commonly used.

### LINCS (Linear Constraint Solver)
GROMACS's default constraint algorithm for all bonds except water. Works by projecting unconstrained positions back onto the constraint surface. Fast, parallel, and handles branched molecules. Cannot handle closed constraint loops (rings).

### SETTLE
An analytical constraint algorithm specifically for rigid three-site water molecules (e.g. TIP3P, SPC). Solves the constraint equations exactly in one step — no iteration needed. Very efficient since water dominates most simulation systems.

### SHAKE
An older iterative constraint algorithm that can handle closed loops (unlike LINCS). Slower and less parallel than LINCS. GROMACS falls back to SHAKE for ring molecules with all bonds constrained.

---

## 6. Electrostatics & Long-Range Interactions

### Coulomb Cutoff
The real-space distance beyond which direct electrostatic pair interactions are not computed. Set by `rcoulomb` (typically 1.0–1.2 nm). With PME, the short-range part is computed directly within the cutoff; the long-range part is computed in reciprocal space.

### Dielectric Constant ($\varepsilon_r$)
A material property describing how strongly it screens electric fields. Vacuum: $\varepsilon_r = 1$. Water at 298 K: $\varepsilon_r \approx 78$. In explicit water MD, $\varepsilon_r = 1$ is used (the screening effect arises naturally from the water molecules).

### Ewald Summation
A method for computing the electrostatic energy of a periodic system by splitting it into a short-range real-space sum and a long-range reciprocal-space sum, each converging rapidly. The basis for PME.

### Fourier Spacing
The `.mdp` parameter `fourierspacing` (nm) controls the PME reciprocal-space grid density. Typical value: 0.12–0.16 nm. Smaller = more grid points = more accurate but slower.

### Particle Mesh Ewald (PME)
The standard method for computing long-range electrostatics in GROMACS (`coulombtype = PME`). Charges are interpolated onto a grid, the electrostatic potential is solved via 3D FFT, and forces are interpolated back to atoms. Complexity: $\mathcal{O}(N \log N)$.

### Reaction Field
An alternative to PME for treating long-range electrostatics (`coulombtype = reaction-field`). Assumes the region beyond the cutoff is a uniform dielectric. Less accurate than PME for heterogeneous systems (e.g. proteins in water) but can be faster.

### Reciprocal Space
The Fourier-transformed representation of real space. In PME, the smooth part of the electrostatic potential is computed in reciprocal space using FFTs. Complementary to the real-space direct sum.

### rvdw
The van der Waals (Lennard-Jones) cutoff distance in nm. Interactions beyond this distance are set to zero (with a potential/force switch near the cutoff, or a dispersion correction applied). Typically 1.0–1.2 nm.

### Screening
The weakening of electrostatic interactions by intervening charges (e.g. solvent molecules, ions). In Ewald summation, point charges are "screened" with Gaussian distributions to make the real-space sum converge faster.

### Switch Function
A mathematical function that smoothly brings an interaction to zero at the cutoff distance, avoiding the discontinuity of a hard cutoff. Used for VdW interactions in some force field setups (e.g. CHARMM recommends a force-switch for LJ: `vdwtype = Cut-off`, `vdw-modifier = Force-switch`, `rvdw-switch = 1.0`).

---

## 7. Periodic Boundary Conditions & Box Geometry

### Box Vectors
The three vectors $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$ that define the simulation box. For a cubic box, $\mathbf{a} = (L, 0, 0)$, $\mathbf{b} = (0, L, 0)$, $\mathbf{c} = (0, 0, L)$. The last line of a `.gro` file contains these vectors (3 values for orthogonal boxes, 9 for triclinic).

### Cubic Box
A rectangular box with $a = b = c$ and all angles 90°. Simple but uses the most volume (least efficient for globular proteins).

### Minimum Image Convention
In a periodic system, each atom interacts with the **nearest image** of every other atom (which may be in a neighbouring replica box rather than the primary box). Requires that the cutoff is smaller than half the smallest box dimension.

### Periodic Boundary Conditions (PBC)
The simulation box is replicated infinitely in all directions. Atoms that leave one face re-enter on the opposite face. Eliminates surface effects, maintains constant density, and enables Ewald summation. Set by `pbc = xyz` in the `.mdp` file.

### Rhombic Dodecahedron
A 12-faced polyhedron that is the most efficient box shape for globular proteins (~71% of the volume of the circumscribing cube). Set with `gmx editconf -bt dodecahedron`.

### Self-Interaction
An artefact where a molecule interacts with its own periodic image through the box boundary. Indicates the box is too small. Causes incorrect forces and structural properties.

### Truncated Octahedron
A 14-faced polyhedron (~77% of the cubic volume). The most space-efficient Bravais lattice for spherical solutes. Set with `gmx editconf -bt octahedron`.

---

## 8. Temperature Control (Thermostats)

### Berendsen Thermostat
First-order temperature coupling: $dT/dt = (T_0 - T)/\tau_T$. Rescales velocities to exponentially relax toward $T_0$. Fast convergence but **does not generate a correct canonical ensemble** — suppresses kinetic energy fluctuations. Use only for equilibration.

### Coupling Group (tc-grps)
A set of atoms assigned to the same thermostat. Typically `tc-grps = Protein Non-Protein` (separate thermostats for protein and solvent). Each group is independently coupled to the target temperature. Groups must cover **all** atoms in the system.

### Equipartition Theorem
Each quadratic degree of freedom contributes $\frac{1}{2}k_BT$ to the average energy: $\langle \frac{1}{2}mv^2 \rangle = \frac{1}{2}k_BT$ per component. Connects temperature to kinetic energy. Only valid for classical systems at equilibrium.

### Flying Ice Cube
A notorious artefact of the Berendsen thermostat where kinetic energy gradually accumulates in centre-of-mass translational motion at the expense of internal (vibrational, rotational) modes. The system looks cold internally while the global kinetic energy is correct. The molecule drifts ("flies") while its atoms freeze. Mitigated by removing COM motion (`comm-mode = Linear`).

### gen_vel / gen_seed / gen_temp
`.mdp` parameters for generating initial velocities. `gen_vel = yes` draws velocities from a Maxwell-Boltzmann distribution at temperature `gen_temp` using random seed `gen_seed`. Different seeds give independent trajectories for replica simulations.

### Nosé-Hoover Thermostat
An extended-system thermostat that adds a friction variable $\xi$ to the equations of motion. Generates a rigorous canonical (NVT) ensemble. Deterministic and time-reversible. Can exhibit oscillatory relaxation — not ideal for equilibration. In GROMACS: `tcoupl = nose-hoover`.

### tau_t (Temperature Coupling Time Constant)
Controls how strongly the thermostat is coupled. Smaller $\tau_T$ = tighter coupling (faster response but more perturbation to dynamics). Typical values: 0.1 ps for V-rescale, 1.0–2.0 ps for Nosé-Hoover. Units: ps.

### Temperature
In MD, the instantaneous temperature is defined as $T = 2E_{\text{kin}}/(N_f k_B)$, where $N_f$ is the number of degrees of freedom. It fluctuates from step to step — this is normal. The time-averaged temperature should equal the target.

### V-rescale (Velocity Rescale)
A stochastic modification of the Berendsen thermostat (Bussi et al., 2007) that produces a **correct canonical ensemble** while retaining Berendsen's fast relaxation. The **recommended thermostat** for most GROMACS workflows. In GROMACS: `tcoupl = V-rescale`.

---

## 9. Pressure Control (Barostats)

### Berendsen Barostat
First-order pressure coupling: rescales box and coordinates to relax pressure toward $P_0$. Fast, monotonic convergence but **does not generate a correct NPT ensemble**. Use only for equilibration.

### Compressibility ($\beta$)
The isothermal compressibility of the medium: $\beta = -\frac{1}{V}\frac{\partial V}{\partial P}$. For liquid water at 300 K: $\beta \approx 4.5 \times 10^{-5}$ bar⁻¹. Set in the `.mdp` file for pressure coupling.

### Density
Mass per unit volume (kg/m³). For a water-solvated protein system at 300 K and 1 bar, the equilibrium density should be ~993–1005 kg/m³ (depending on the water model). Monitored during NPT equilibration using `gmx energy`.

### Isotropic Pressure Coupling
All box dimensions are scaled equally (`pcoupltype = isotropic`). Appropriate for globular proteins in solution.

### MTTK
Martyna-Tuckerman-Tobias-Klein barostat. A velocity-Verlet-based barostat that generates a correct NPT ensemble. Used with `integrator = md-vv`. Less common than Parrinello-Rahman.

### Parrinello-Rahman Barostat
The box vectors are treated as dynamical variables with their own equations of motion. Generates a **correct NPT ensemble** and allows anisotropic box deformation. The standard choice for production MD. Can oscillate if started far from equilibrium — always equilibrate with Berendsen first. In GROMACS: `pcoupl = Parrinello-Rahman`.

### Pressure
In MD, the instantaneous pressure is computed from the virial: $P = \frac{1}{V}[Nk_BT + \frac{1}{3}\sum_{i<j}\mathbf{r}_{ij}\cdot\mathbf{F}_{ij}]$. Extremely noisy (fluctuations of ±100–500 bar are normal at 1 bar). Only the time-averaged pressure is meaningful.

### Semi-Isotropic Pressure Coupling
$x$/$y$ dimensions scale together but independently from $z` (`pcoupltype = semiisotropic`). Essential for membrane simulations where the membrane plane (xy) and its normal (z) have different compressibility.

### Surface Tension Coupling
A variant of pressure coupling where a target surface tension (γ) is maintained for the xy plane. Used in membrane simulations to investigate area-dependent properties: `pcoupltype = surface-tension`.

### tau_p (Pressure Coupling Time Constant)
Controls barostat strength. Typical values: 1.0–2.0 ps for Berendsen, 5.0 ps for Parrinello-Rahman. Larger values = gentler coupling, less perturbation to dynamics.

### Virial
The sum of $\mathbf{r}_{ij} \cdot \mathbf{F}_{ij}$ over all interacting pairs. Quantifies the contribution of interatomic forces to pressure. Dominates the pressure in a condensed-phase system. Computed internally by GROMACS at each step.

---

## 10. Statistical Mechanics & Thermodynamics

### Autocorrelation Function
$C(\tau) = \langle A(t) A(t + \tau) \rangle / \langle A^2 \rangle$. Measures how a property correlates with itself at a later time. Used to determine correlation times (how long the system "remembers" its state) and to compute transport coefficients.

### Block Averaging
A technique for estimating error bars from correlated time series data. The trajectory is split into blocks; the variance of block averages gives the standard error. Block size should exceed the correlation time.

### Boltzmann Distribution
The probability of a system being in a microstate with energy $E$: $P(E) \propto e^{-E/k_BT}$. Sampled correctly by the canonical ensemble. The foundation of equilibrium statistical mechanics.

### Canonical Ensemble (NVT)
Fixed number of particles (N), volume (V), and temperature (T). Generated by correct thermostats (V-rescale, Nosé-Hoover). Energy fluctuates. The partition function is $Z = \sum_i e^{-E_i/k_BT}$.

### Correlation Time
The time over which a property remains significantly correlated with its own past. For a property with correlation time $\tau_c$, you need $t_{\text{sim}} \gg \tau_c$ for adequate sampling. Effective independent samples ≈ $t_{\text{sim}} / (2\tau_c)$.

### Ensemble
A (conceptually infinite) collection of microstates consistent with the same macroscopic constraints. In MD, different ensembles are generated by different combinations of thermostats and barostats (NVE, NVT, NPT).

### Entropy
$S = -k_B \sum_i P_i \ln P_i$ (Gibbs entropy). A measure of the number of accessible microstates. Cannot be directly computed from a single MD trajectory. Contributes to free energy: $G = H - TS$.

### Ergodic Hypothesis
The assumption that a system will visit all accessible microstates given sufficient time, so that time averages equal ensemble averages. Fails when the simulation is too short to cross energy barriers between states — the fundamental **sampling problem** in MD.

### Fluctuation
The instantaneous deviation of a property from its average: $\delta A = A - \langle A \rangle$. Fluctuations are physical and related to response functions: $\langle (\delta E)^2 \rangle_{NVT} = k_BT^2 C_V$ (heat capacity from energy fluctuations), $\langle (\delta V)^2 \rangle_{NPT} = k_BT V \beta$ (compressibility from volume fluctuations).

### Free Energy ($G$, $F$)
The thermodynamic potential that determines spontaneous processes. Gibbs free energy: $G = H - TS$ (constant $P$). Helmholtz free energy: $F = U - TS$ (constant $V$). **Cannot be computed directly from MD** — only free energy *differences* are accessible via perturbation, integration, or other methods.

### Isothermal-Isobaric Ensemble (NPT)
Fixed N, pressure (P), and temperature (T). The standard ensemble for biomolecular simulation, matching experimental conditions. Volume fluctuates. Generated by combining a thermostat and a barostat.

### Microcanonical Ensemble (NVE)
Fixed N, V, and total energy (E). What you get with no thermostat or barostat — pure Newtonian dynamics. Temperature fluctuates. Useful for testing integrator accuracy (energy should be conserved).

### Microstate
A single point in phase space — a complete specification of all atomic positions and momenta. Macroscopic observables are averages over microstates.

### Partition Function ($Z$)
$Z = \sum_i e^{-E_i/k_BT}$ (NVT). A sum over all possible microstates. Contains all thermodynamic information ($F = -k_BT \ln Z$) but is practically impossible to compute for a macroscopic system. Free energy *differences* bypass the need for absolute $Z$.

### Sampling
The process of exploring phase space (different configurations) during a simulation. **Adequate sampling** means the simulation has visited all relevant microstates enough times to compute reliable averages. The sampling problem is the most fundamental limitation of MD.

---

## 11. Energy Minimisation

### Conjugate Gradient
An optimisation algorithm that uses information from the previous step to choose a better search direction (conjugate to the previous gradient). Converges faster than steepest descent near a minimum but is less robust for severely overlapping atoms. In GROMACS: `integrator = cg`.

### emstep
The initial step size (nm) for steepest descent minimisation. Typically 0.01 nm. GROMACS adjusts this adaptively during minimisation.

### emtol
The maximum force tolerance for convergence (kJ/mol/nm). Minimisation stops when the maximum force on any atom falls below `emtol`. Typical: 1,000 kJ/mol/nm (standard) or 500 (strict).

### Global Minimum
The lowest-energy point on the entire potential energy surface. For proteins, this is (arguably) the native state. Energy minimisation **cannot** find the global minimum — it always converges to the nearest **local minimum**.

### Local Minimum
A point where the energy is lower than all immediately neighbouring configurations, but not necessarily the lowest energy overall. Proteins have astronomically many local minima. Standard minimisation always finds a local minimum near the starting structure.

### Steepest Descent
The simplest optimisation algorithm: move each atom in the direction of the force (downhill in energy). Robust but slow to converge near a minimum (zig-zag path). In GROMACS: `integrator = steep`.

---

## 12. Solvation, Ions & System Setup

### Counter-Ion
An ion added to neutralise the net charge of the system (e.g. Na⁺ to balance negative protein charge, or Cl⁻ to balance positive charge). Required for PME.

### editconf
GROMACS tool (`gmx editconf`) for editing structure files. Primary use: defining the simulation box shape and size (`-bt cubic/dodecahedron/octahedron`, `-d 1.0` for edge distance).

### genion
GROMACS tool (`gmx genion`) that replaces solvent molecules with ions. Used to neutralise the system (`-neutral`) and optionally add salt (`-conc 0.15` for 0.15 M). Updates the topology `[ molecules ]` automatically.

### Ionic Strength
The concentration of dissolved ions, affecting electrostatic screening (Debye-Hückel theory). Physiological ionic strength: ~0.15 M NaCl. Set with `gmx genion -conc 0.15`.

### pdb2gmx
The GROMACS tool that converts a PDB file to GROMACS format. It:
1. Selects the force field and water model.
2. Assigns atom types and partial charges.
3. Adds hydrogen atoms.
4. Assigns protonation states (interactively for HIS).
5. Generates `.gro` (structure), `.top` (topology), and `.itp` (molecule definition + position restraints) files.

### Solvation
The process of adding explicit water molecules around the solute to fill the simulation box. Performed by `gmx solvate`. The solvent equilibrates around the protein during NVT/NPT equilibration.

### solvate
GROMACS tool (`gmx solvate`) that fills the simulation box with solvent molecules (typically water). Reads a pre-equilibrated water box and tiles/trims it to fill the available space. Updates the `[ molecules ]` section of the topology.

---

## 13. Trajectory & Analysis Terminology

### Cluster Analysis
Grouping trajectory frames into structurally similar sets. Common algorithms: GROMOS clustering (rmsd-based), k-means, hierarchical. Used to identify distinct conformational states. In GROMACS: `gmx cluster`.

### Contact Map
A 2D matrix showing which residue pairs are spatially close (typically Cα–Cα distance < 8 Å). Diagonal = trivial self-contacts. Off-diagonal patterns reveal the protein fold, secondary structure, and tertiary contacts.

### DSSP
Define Secondary Structure of Proteins. An algorithm (and program) that assigns secondary structure types (H = helix, E = strand, T = turn, etc.) from atomic coordinates based on hydrogen bonding patterns. Used by `gmx do_dssp` (requires external DSSP installation).

### Frame
A single snapshot of the system at one point in time within a trajectory. Each frame contains all atomic coordinates (and optionally velocities/forces). Writing frequency is controlled by `nstxout-compressed`.

### Mean Squared Displacement (MSD)
$\text{MSD}(\tau) = \langle |\mathbf{r}(t+\tau) - \mathbf{r}(t)|^2 \rangle$. Used to compute the diffusion coefficient $D$ via the Einstein relation: $D = \lim_{\tau \to \infty} \text{MSD}(\tau) / (6\tau)$. Computed by `gmx msd`.

### Principal Component Analysis (PCA)
A dimensionality reduction technique that identifies the dominant modes of motion (largest-amplitude collective fluctuations) from the covariance matrix of atomic coordinates. The first few principal components (eigenvectors) often capture biologically relevant motions. In GROMACS: `gmx covar` (build covariance matrix) → `gmx anaeig` (project and analyse).

### Radius of Gyration ($R_g$)
A measure of the compactness of a protein: $R_g = \sqrt{\frac{\sum_i m_i |\mathbf{r}_i - \mathbf{r}_{\text{COM}}|^2}{\sum_i m_i}}$. Increases upon unfolding. Computed by `gmx gyrate`.

### RMSD (Root Mean Square Deviation)
$\text{RMSD} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}|\mathbf{r}_i - \mathbf{r}_i^{\text{ref}}|^2}$. Measures the average deviation of atomic positions from a reference structure after optimal superposition (least-squares fitting). Commonly used to track structural stability over time. Computed by `gmx rms`.

### RMSF (Root Mean Square Fluctuation)
$\text{RMSF}_i = \sqrt{\langle |\mathbf{r}_i - \langle \mathbf{r}_i \rangle|^2 \rangle}$. The time-averaged fluctuation of each atom (or residue) around its mean position. A per-residue measure of flexibility — analogous to the crystallographic B-factor. Computed by `gmx rmsf`.

### Superposition (Fitting / Alignment)
The process of rotating and translating one structure to minimise the RMSD relative to a reference. Essential before computing RMSD or RMSF to remove rigid-body translation and rotation. Typically performed on backbone or Cα atoms.

### trjconv
The GROMACS tool for trajectory conversion and manipulation (`gmx trjconv`). Essential for: fixing PBC artefacts (`-pbc whole/mol/nojump`), centring molecules (`-center`), extracting time ranges (`-b`, `-e`), extracting specific atom groups, and format conversion.

### XVG File
A plain-text data file produced by GROMACS analysis tools. Contains columns of data (e.g. time vs. RMSD). Includes Grace/xmgrace header lines (starting with `@`). Can be plotted with xmgrace, or loaded into Python/gnuplot by skipping comment lines.

---

## 14. Free Energy Methods

### Alchemical Transformation
A computational technique where a molecule is gradually "mutated" into another (e.g. drug A → drug B, or solvated ligand → nothing) by introducing a coupling parameter $\lambda$ that interpolates between the two states. Used to compute binding and solvation free energies.

### BAR (Bennett Acceptance Ratio)
The most statistically efficient method for computing free energy differences from overlapping forward and reverse perturbation data. Implemented in GROMACS via `gmx bar`.

### Coupling Parameter ($\lambda$)
A variable that interpolates between two thermodynamic states in alchemical free energy calculations. $\lambda = 0$: initial state. $\lambda = 1$: final state. The simulation is run at multiple $\lambda$ values.

### FEP (Free Energy Perturbation)
$\Delta G = -k_BT \ln \langle e^{-\Delta V/k_BT} \rangle_A$. Computes the free energy difference between state A and state B by exponential averaging of energy differences on trajectories at state A. Requires good overlap between states.

### MM/PBSA
Molecular Mechanics / Poisson-Boltzmann Surface Area. An approximate endpoint method for estimating binding free energies from MD snapshots. Treats the binding energy as: $\Delta G_{\text{bind}} \approx \Delta E_{\text{MM}} + \Delta G_{\text{solv}} - T\Delta S$. Fast but approximate — the entropy term is often neglected or estimated crudely.

### MM/GBSA
Similar to MM/PBSA but uses the Generalised Born (GB) model for solvation instead of Poisson-Boltzmann. Faster but less accurate.

### Soft-Core Potential
A modified Lennard-Jones potential used during alchemical transformations to prevent singularities when atoms are appearing or disappearing. The $r^{-12}$ repulsion is smoothed at short distances to avoid "particle overlap" at intermediate $\lambda$ values.

### Thermodynamic Cycle
A closed cycle of state transformations used to compute quantities that are difficult to measure directly. Example: to compute the relative binding free energy of two ligands, compute the alchemical transformation of one into the other both in solvent and in the protein.

### Thermodynamic Integration (TI)
$\Delta G = \int_0^1 \langle \partial V/\partial \lambda \rangle_\lambda \, d\lambda$. Computes free energy differences by integrating the ensemble average of $\partial V/\partial \lambda$ over the coupling parameter. Requires simulations at multiple $\lambda$ values.

---

## 15. Enhanced Sampling Methods

### Accelerated MD (aMD / GaMD)
Methods that add a boost potential to fill energy wells, allowing the system to escape kinetic traps more easily. Gaussian accelerated MD (GaMD) gives a Gaussian-shaped boost, allowing reweighting to recover the unbiased ensemble.

### Collective Variable (CV)
A reduced-dimension coordinate that describes a slow process of interest (e.g. a distance, angle, RMSD, coordination number). Used in metadynamics, umbrella sampling, and other biasing methods to drive sampling along relevant degrees of freedom.

### Metadynamics
An enhanced sampling method that iteratively adds Gaussian bias potentials along one or more collective variables, gradually filling energy minima and encouraging the system to explore new regions. Converges to the negative of the free energy surface. Implemented via the **PLUMED** plugin.

### PLUMED
An open-source plugin for free energy calculations and enhanced sampling that interfaces with GROMACS and other MD engines. Provides metadynamics, steered MD, umbrella sampling, and analysis of collective variables. Must be compiled with GROMACS.

### Potential of Mean Force (PMF)
The free energy as a function of one or more collective variables: $W(\xi) = -k_BT \ln P(\xi)$, where $P(\xi)$ is the probability distribution along $\xi$. The "landscape" that a reaction coordinate traverses.

### Replica Exchange MD (REMD / REST2)
Multiple copies (replicas) of the system are simulated at different temperatures (temperature REMD) or with different Hamiltonians (Hamiltonian REMD / REST2). Periodically, configurations are swapped between replicas based on a Metropolis criterion. High-temperature replicas cross barriers easily; swaps allow low-temperature replicas to access these configurations.

### Steered MD (SMD)
An external force is applied to pull the system along a collective variable (e.g. pulling a ligand out of a binding site). Produces non-equilibrium work trajectories from which free energies can be estimated via the Jarzynski equality.

### Umbrella Sampling
A biasing method where multiple simulations (windows) are run with harmonic restraints along a collective variable. Each window samples a narrow region. The biased distributions are unbiased and combined using **WHAM** (Weighted Histogram Analysis Method) to reconstruct the full free energy profile (PMF). In GROMACS: `gmx wham`.

### WHAM (Weighted Histogram Analysis Method)
A statistical method for combining data from multiple umbrella sampling windows to reconstruct an unbiased free energy profile. In GROMACS: `gmx wham`.

---

## 16. Coarse-Graining

### All-Atom (AA)
A molecular representation where every atom (including hydrogens) is treated explicitly. Maximum detail but highest computational cost.

### Back-Mapping
The process of converting a coarse-grained structure back to atomistic resolution. Involves placing atoms within CG beads and re-equilibrating. Tools: `backward.py` (MARTINI), `initram.sh`. The back-mapped structure requires energy minimisation and short equilibration.

### Bead
A single interaction site in a coarse-grained model that represents a group of atoms (typically 2–6 heavy atoms). Heavier and smoother than atoms, enabling larger time steps and faster simulations.

### Bottom-Up Coarse-Graining
Parameterisation of CG models from atomistic simulations (reference data). Methods include iterative Boltzmann inversion (IBI), force matching, and relative entropy minimisation. Produces models that reproduce atomistic structural properties.

### Mapping
The scheme that defines which atoms belong to which CG bead. MARTINI uses roughly a 4:1 heavy atom-to-bead mapping. The mapping determines the resolution and transferability of the CG model.

### MARTINI
A widely used coarse-grained force field for biomolecular simulations. Approximately 4 heavy atoms per bead. Parameterised top-down (fitted to partition free energies). Current version: MARTINI 3. Excellent for membrane systems and large-scale dynamics.

### Top-Down Coarse-Graining
Parameterisation of CG models from macroscopic experimental data (e.g. partition coefficients, surface tensions, densities). MARTINI uses this approach. Produces models with correct thermodynamic properties.

### United Atom (UA)
An intermediate representation where nonpolar hydrogens (CH, CH₂, CH₃) are subsumed into their parent heavy atom. Fewer sites than all-atom but more than CG. Used by GROMOS.

---

## 17. GROMACS-Specific Commands & Tools

### gmx bar
Compute free energy differences using the Bennett Acceptance Ratio method from `.xvg` energy difference files.

### gmx check
Verify the integrity of trajectory (`.xtc`, `.trr`) and topology (`.tpr`) files. Reports number of frames, time range, and any detected errors.

### gmx cluster
Perform conformational clustering on a trajectory to identify distinct structural states.

### gmx covar
Build the covariance matrix of atomic fluctuations for Principal Component Analysis.

### gmx anaeig
Analyse the eigenvectors from PCA: project trajectories onto principal components, compute cosine content (a diagnostic for insufficient sampling — high cosine content ≈ random diffusion).

### gmx do_dssp
Compute secondary structure assignment over the trajectory using the DSSP algorithm.

### gmx dump
Print a human-readable representation of binary GROMACS files (`.tpr`, `.cpt`, `.edr`). Essential for debugging: `gmx dump -s md.tpr | head -100`.

### gmx editconf
Edit structure files: change box shape/size, centre the molecule, convert formats, translate/rotate coordinates.

### gmx energy
Extract thermodynamic quantities (temperature, pressure, density, potential energy, kinetic energy, etc.) from the `.edr` energy file and write them as `.xvg` time series.

### gmx genion
Replace selected solvent molecules with ions. Used to neutralise charge and/or add salt.

### gmx grompp
The GROMACS pre-processor. Combines `.mdp` (parameters), `.gro` (coordinates), `.top` (topology), and optionally `.ndx` (index) into a `.tpr` (portable run input) file. Performs extensive error checking — **always read its output carefully**.

### gmx gyrate
Compute the radius of gyration over time.

### gmx hbond
Analyse hydrogen bonds between selected groups (number, occupancy, lifetime).

### gmx make_ndx
Create custom index files (`.ndx`) by defining atom groups interactively or programmatically.

### gmx mdrun
The MD engine — runs the simulation defined in the `.tpr` file. Supports MPI parallelism, GPU offloading, and checkpoint restart.

### gmx mindist
Compute the minimum distance between two groups over time. Also reports contact numbers within a cutoff.

### gmx msd
Compute mean squared displacement and derive diffusion coefficients.

### gmx pdb2gmx
Convert a PDB file to GROMACS format, selecting force field and water model, assigning atom types and charges, adding hydrogens, and generating topology files.

### gmx rms
Compute RMSD from a reference structure over time.

### gmx rmsf
Compute per-atom or per-residue root mean square fluctuation (RMSF).

### gmx sasa
Compute solvent-accessible surface area over time.

### gmx solvate
Fill the simulation box with solvent molecules.

### gmx trjcat
Concatenate multiple trajectory files into one. Use `-settime` to adjust time offsets.

### gmx trjconv
Convert and manipulate trajectories: fix PBC, centre molecules, extract frames, change output format, select atom groups.

### gmx wham
Perform Weighted Histogram Analysis Method to compute a potential of mean force from umbrella sampling data.

### grompp warnings and notes
`grompp` output is classified as:
- **NOTE:** Informational — usually safe to ignore.
- **WARNING:** Potential issue — read carefully, may or may not be a problem.
- **ERROR:** Fatal problem — the `.tpr` will not be generated until it is fixed.

Always read the full `grompp` output, especially the warnings.

---

## 18. File Formats

### `.cif` (mmCIF / PDBx)
The current standard format for macromolecular structures. Dictionary-driven, no column-width limits, supports large complexes. Master format of the PDB since 2014.

### `.cpt` (Checkpoint)
Binary file containing the full simulation state (coordinates, velocities, coupling variables, RNG state, step number). Written periodically by `mdrun` for exact restart.

### `.edr` (Energy)
Binary file containing all energy terms, temperature, pressure, density, and other properties as time series. Read by `gmx energy`.

### `.gro` (GROMACS Coordinate)
Fixed-width text format for atomic coordinates. Uses nm. Contains: title, atom count, atom records (residue number, name, atom name, serial number, x, y, z, [vx, vy, vz]), and box vectors. Limited to 5-digit residue and atom numbers.

### `.itp` (Include Topology)
Text file defining a single molecule type: `[moleculetype]`, `[atoms]`, `[bonds]`, `[pairs]`, `[angles]`, `[dihedrals]`. Included by the `.top` file via `#include`.

### `.log` (Simulation Log)
Human-readable log file from `mdrun`. Contains simulation parameters, step-by-step energy output, performance statistics, and any warnings.

### `.mdp` (Molecular Dynamics Parameters)
User-written text file specifying all simulation parameters: integrator, time step, cutoffs, thermostat, barostat, output frequencies, etc. Consumed by `grompp`.

### `.ndx` (Index)
Text file defining named groups of atom indices (`[GroupName]` followed by atom numbers). Used to specify custom selections for analysis or coupling groups.

### `.pdb` (Protein Data Bank)
Fixed-column-width text format (80 chars/line) for atomic coordinates. Uses Ångströms. The traditional (but retiring) interchange format for macromolecular structures. Max 99,999 atoms, 62 chains.

### `.top` (System Topology)
The master topology file. Contains `#include` directives for force field parameters and molecule definitions, plus the `[system]` and `[molecules]` sections. Defines the *complete* interacting system.

### `.tpr` (Portable Run Input)
Binary file containing everything needed for `mdrun`: force field parameters, coordinates, velocities, and simulation settings. Self-contained and portable. Created by `grompp`.

### `.trr` (Full-Precision Trajectory)
Binary trajectory format storing coordinates, velocities, and/or forces at double precision. Large files. Used when velocities or forces are needed (e.g. for restart or force analysis).

### `.xtc` (Compressed Trajectory)
Binary trajectory format storing coordinates only, with lossy compression (~0.001 nm precision). The standard format for routine analysis. ~10× smaller than `.trr`.

### `.xvg` (XmGrace Data)
Text data file with Grace/xmgrace formatting headers. Columns of numerical data (e.g. time vs. observable). Produced by most GROMACS analysis tools.

---

## 19. HPC & Performance Terminology

### Checkpoint Restart
Resuming a simulation from a `.cpt` file: `gmx mdrun -cpi md.cpt`. Exact continuation of the run — same trajectory as if it had never been interrupted.

### Core-Hours
A measure of computational cost: number of CPU cores × wall-clock hours. Used for resource allocation on HPC clusters.

### Domain Decomposition (DD)
GROMACS's spatial parallelisation method. The simulation box is divided into rectangular sub-domains, each assigned to an MPI rank. Atoms near domain boundaries exchange information (halo exchange). The decomposition adjusts dynamically for load balancing.

### GPU Offloading
Moving computationally expensive tasks from the CPU to a GPU. In GROMACS: nonbonded forces (`-nb gpu`), PME (`-pme gpu`), bonded forces (`-bonded gpu`), and the update/constraint step (`-update gpu`).

### Load Balancing
The process of distributing work evenly across processors. GROMACS performs dynamic load balancing by shifting domain boundaries. The simulation log reports the load imbalance percentage.

### MPI (Message Passing Interface)
A standard for inter-process communication in parallel computing. GROMACS uses MPI for multi-node parallelism. Each MPI rank handles one domain decomposition cell.

### ns/day
Nanoseconds of simulated time per wall-clock day. The standard performance metric for MD. Higher is better. Check at the end of the `.log` file.

### ntmpi / ntomp
`-ntmpi`: number of MPI thread ranks (inter-process parallelism, typically = number of GPUs). `-ntomp`: number of OpenMP threads per rank (intra-process parallelism, for CPU work). Optimal ratio depends on hardware.

### OpenMP
A shared-memory parallelisation standard. GROMACS uses OpenMP for multi-threaded work within each MPI rank (e.g. parallel force reduction). Controlled by `-ntomp`.

### Thread-MPI
A lightweight MPI implementation used by GROMACS when the real MPI library is not available. Allows multi-rank execution on a single machine without installing a full MPI library.

### Wall Time
The actual elapsed time (hours, minutes) of a job on a cluster. Distinguished from CPU time (core-hours). Cluster jobs have maximum wall-time limits.

---

## 20. Water Models

### CHARMM-TIP3P (TIPS3P)
A modified TIP3P model with Lennard-Jones interactions on hydrogen atoms. **Must** be used with the CHARMM force field. Slightly different properties from standard TIP3P.

### OPC (Optimal Point Charge)
A 4-site water model with an off-atom charge site. Recommended with AMBER ff19SB, especially for intrinsically disordered proteins. Produces better water dielectric and diffusion properties than TIP3P.

### SPC (Simple Point Charge)
A 3-site water model with a pure O–H bond angle of 109.47° (tetrahedral). Used with GROMOS force fields. Simple and fast but less accurate than SPC/E.

### SPC/E (Extended Simple Point Charge)
SPC with a self-polarisation correction. Better density, diffusion, and dielectric properties than SPC. Compatible with many force fields.

### TIP3P (Transferable Intermolecular Potential, 3 Points)
The most widely used water model. Three interaction sites: O (LJ + charge) and two H (charge only). Fast to compute but overestimates diffusion and underestimates viscosity. Default for most AMBER simulations.

### TIP4P
A 4-site model with an off-atom negative charge site (M-site) along the H–O–H bisector. Better electrostatic properties than TIP3P but ~30% slower (extra site). Available as TIP4P, TIP4P-Ew (for Ewald), TIP4P/2005 (best overall phase diagram).

### TIP4P/2005
A reparameterised TIP4P model that reproduces the water phase diagram, temperature of maximum density, and many other properties better than any other rigid, non-polarisable model. The recommended 4-site model for general use.

### TIP5P
A 5-site model with two lone-pair sites. Better tetrahedral structure and orientation properties but 60–80% slower than TIP3P. Rarely used in biomolecular work.

---

## 21. Common Abbreviations

| Abbreviation | Full Form |
|---|---|
| **AA** | All-Atom |
| **aMD** | Accelerated Molecular Dynamics |
| **BAR** | Bennett Acceptance Ratio |
| **CG** | Coarse-Grained |
| **COM** | Centre of Mass |
| **cryo-EM** | Cryo-Electron Microscopy |
| **CV** | Collective Variable / Heat Capacity |
| **DD** | Domain Decomposition |
| **DFT** | Density Functional Theory |
| **DSSP** | Define Secondary Structure of Proteins |
| **EM** | Energy Minimisation |
| **ESP** | Electrostatic Potential |
| **FEP** | Free Energy Perturbation |
| **FF** | Force Field |
| **FFT** | Fast Fourier Transform |
| **GAFF** | Generalised AMBER Force Field (for small molecules) |
| **GaMD** | Gaussian Accelerated Molecular Dynamics |
| **GPU** | Graphics Processing Unit |
| **GROMACS** | GROningen MAchine for Chemical Simulations |
| **HMR** | Hydrogen Mass Repartitioning |
| **HPC** | High-Performance Computing |
| **IDP** | Intrinsically Disordered Protein |
| **LJ** | Lennard-Jones |
| **LINCS** | Linear Constraint Solver |
| **MARTINI** | (Not an acronym — named after the cocktail) |
| **MD** | Molecular Dynamics |
| **MM/PBSA** | Molecular Mechanics / Poisson-Boltzmann Surface Area |
| **MMCIF** | Macromolecular Crystallographic Information Framework |
| **MPI** | Message Passing Interface |
| **MSD** | Mean Squared Displacement |
| **NMR** | Nuclear Magnetic Resonance |
| **NPT** | Isothermal-Isobaric Ensemble (constant N, P, T) |
| **NVE** | Microcanonical Ensemble (constant N, V, E) |
| **NVT** | Canonical Ensemble (constant N, V, T) |
| **PBC** | Periodic Boundary Conditions |
| **PCA** | Principal Component Analysis |
| **PDB** | Protein Data Bank |
| **PES** | Potential Energy Surface |
| **PME** | Particle Mesh Ewald |
| **PMF** | Potential of Mean Force |
| **PR** | Parrinello-Rahman (barostat) |
| **QM** | Quantum Mechanics |
| **QM/MM** | Quantum Mechanics / Molecular Mechanics |
| **REMD** | Replica Exchange Molecular Dynamics |
| **RESP** | Restrained Electrostatic Potential |
| **RMSD** | Root Mean Square Deviation |
| **RMSF** | Root Mean Square Fluctuation |
| **SASA** | Solvent-Accessible Surface Area |
| **SETTLE** | (Analytical constraint solver for water) |
| **SHAKE** | (Iterative constraint solver — not an acronym) |
| **SMD** | Steered Molecular Dynamics |
| **SPC** | Simple Point Charge |
| **TI** | Thermodynamic Integration |
| **TIP3P** | Transferable Intermolecular Potential, 3 Points |
| **UA** | United Atom |
| **VdW** | Van der Waals |
| **VMD** | Visual Molecular Dynamics (visualisation software) |
| **WHAM** | Weighted Histogram Analysis Method |

---

## Quick Reference: Units in GROMACS

| Quantity | Unit | Symbol |
|---|---|---|
| Length | nanometre | nm |
| Time | picosecond | ps |
| Mass | atomic mass unit | amu (Da) |
| Charge | elementary charge | $e$ |
| Energy | kJ/mol | kJ mol⁻¹ |
| Force | kJ/(mol·nm) | kJ mol⁻¹ nm⁻¹ |
| Temperature | Kelvin | K |
| Pressure | bar | bar |
| Velocity | nm/ps | nm ps⁻¹ |
| Force constant (bond) | kJ/(mol·nm²) | kJ mol⁻¹ nm⁻² |
| Force constant (angle) | kJ/(mol·rad²) | kJ mol⁻¹ rad⁻² |

---

*Glossary compiled for the QMUL Molecular Dynamics Teaching Project. Last updated: March 2026.*