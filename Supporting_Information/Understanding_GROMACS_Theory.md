# Understanding GROMACS Theory

**The Biophysics and Computational Science Behind Molecular Dynamics Simulation**

---

## Table of Contents

1. [Why Molecular Dynamics?](#1-why-molecular-dynamics)
2. [The Central Idea: Classical Mechanics on Atoms](#2-the-central-idea-classical-mechanics-on-atoms)
3. [From Quantum Mechanics to Classical Approximation](#3-from-quantum-mechanics-to-classical-approximation)
   - 3.1 [The Born-Oppenheimer Approximation](#31-the-born-oppenheimer-approximation)
   - 3.2 [From Electrons to Point Charges](#32-from-electrons-to-point-charges)
   - 3.3 [What We Gain and What We Lose](#33-what-we-gain-and-what-we-lose)
4. [The Potential Energy Function (Force Field)](#4-the-potential-energy-function-force-field)
   - 4.1 [Bonded Interactions in Detail](#41-bonded-interactions-in-detail)
   - 4.2 [Nonbonded Interactions in Detail](#42-nonbonded-interactions-in-detail)
   - 4.3 [Combining Rules](#43-combining-rules)
5. [Newton's Equations of Motion](#5-newtons-equations-of-motion)
6. [Numerical Integration: The Leap-Frog Algorithm](#6-numerical-integration-the-leap-frog-algorithm)
   - 6.1 [Why Not Just Solve the Equations Analytically?](#61-why-not-just-solve-the-equations-analytically)
   - 6.2 [The Verlet Family of Integrators](#62-the-verlet-family-of-integrators)
   - 6.3 [Leap-Frog in Practice](#63-leap-frog-in-practice)
   - 6.4 [Time Step Selection](#64-time-step-selection)
7. [Constraint Algorithms: LINCS and SETTLE](#7-constraint-algorithms-lincs-and-settle)
   - 7.1 [Why Constrain Bonds?](#71-why-constrain-bonds)
   - 7.2 [LINCS (Linear Constraint Solver)](#72-lincs-linear-constraint-solver)
   - 7.3 [SETTLE (for Water)](#73-settle-for-water)
8. [Handling Long-Range Electrostatics: Particle Mesh Ewald](#8-handling-long-range-electrostatics-particle-mesh-ewald)
   - 8.1 [The Problem with Coulomb's Law](#81-the-problem-with-coulombs-law)
   - 8.2 [Ewald Summation — The Idea](#82-ewald-summation--the-idea)
   - 8.3 [Particle Mesh Ewald (PME) — Making It Fast](#83-particle-mesh-ewald-pme--making-it-fast)
   - 8.4 [PME Parameters in GROMACS](#84-pme-parameters-in-gromacs)
9. [Periodic Boundary Conditions](#9-periodic-boundary-conditions)
   - 9.1 [The Minimum Image Convention](#91-the-minimum-image-convention)
   - 9.2 [Box Shapes](#92-box-shapes)
   - 9.3 [Artefacts and Limitations](#93-artefacts-and-limitations)
10. [Neighbour Searching: The Verlet Cutoff Scheme](#10-neighbour-searching-the-verlet-cutoff-scheme)
11. [Temperature Control (Thermostats)](#11-temperature-control-thermostats)
    - 11.1 [Why Control Temperature?](#111-why-control-temperature)
    - 11.2 [The Equipartition Theorem and Kinetic Energy](#112-the-equipartition-theorem-and-kinetic-energy)
    - 11.3 [Berendsen Thermostat](#113-berendsen-thermostat)
    - 11.4 [Velocity-Rescale (V-rescale) Thermostat](#114-velocity-rescale-v-rescale-thermostat)
    - 11.5 [Nosé-Hoover Thermostat](#115-nosé-hoover-thermostat)
    - 11.6 [Which Thermostat Should I Use?](#116-which-thermostat-should-i-use)
12. [Pressure Control (Barostats)](#12-pressure-control-barostats)
    - 12.1 [Why Control Pressure?](#121-why-control-pressure)
    - 12.2 [The Virial and Pressure](#122-the-virial-and-pressure)
    - 12.3 [Berendsen Barostat](#123-berendsen-barostat)
    - 12.4 [Parrinello-Rahman Barostat](#124-parrinello-rahman-barostat)
    - 12.5 [Which Barostat Should I Use?](#125-which-barostat-should-i-use)
13. [Statistical Mechanics: Connecting Simulation to Experiment](#13-statistical-mechanics-connecting-simulation-to-experiment)
    - 13.1 [Ensembles](#131-ensembles)
    - 13.2 [The Ergodic Hypothesis](#132-the-ergodic-hypothesis)
    - 13.3 [Free Energy and Sampling](#133-free-energy-and-sampling)
14. [Energy Minimisation](#14-energy-minimisation)
    - 14.1 [Steepest Descent](#141-steepest-descent)
    - 14.2 [Conjugate Gradient](#142-conjugate-gradient)
    - 14.3 [Why Minimise Before MD?](#143-why-minimise-before-md)
15. [Equilibration: NVT and NPT](#15-equilibration-nvt-and-npt)
16. [Key Assumptions of Molecular Dynamics](#16-key-assumptions-of-molecular-dynamics)
17. [Known Limitations and When MD Can Fail](#17-known-limitations-and-when-md-can-fail)
18. [Summary Table of Key Algorithms in GROMACS](#18-summary-table-of-key-algorithms-in-gromacs)
19. [Glossary](#19-glossary)
20. [Further Reading & Next Steps](#20-further-reading--next-steps)

---

## 1. Why Molecular Dynamics?

Proteins are not static sculptures. They breathe, flex, bind, and unfold on timescales from femtoseconds (bond vibrations) to seconds (folding). Experimental methods — X-ray crystallography, cryo-EM, NMR — give us snapshots or ensemble averages, but rarely the full dynamic movie.

**Molecular Dynamics (MD) simulation** fills this gap. It computes how every atom in a system moves over time by solving Newton's equations of motion, given a model for the forces between atoms (the force field). The result is a **trajectory** — a time-ordered series of atomic coordinates — from which we can extract:

- Conformational changes and protein flexibility
- Binding free energies and pathways
- Membrane behaviour and lipid–protein interactions
- Thermodynamic properties (entropy, heat capacity)
- Transport properties (diffusion coefficients, viscosity)

GROMACS (GROningen MAchine for Chemical Simulations) is one of the fastest and most widely used MD engines, particularly strong for biomolecular systems.

---

## 2. The Central Idea: Classical Mechanics on Atoms

At its core, MD is strikingly simple:

1. **Define a potential energy function** $V(\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_N)$ that tells you the energy of the system for any configuration of $N$ atoms.
2. **Calculate the force** on each atom: $\mathbf{F}_i = -\nabla_i V$
3. **Integrate Newton's second law** to get new positions and velocities: $\mathbf{F}_i = m_i \mathbf{a}_i$
4. **Advance time** by a small step $\Delta t$ and repeat.

That is the entire algorithm. Everything else — thermostats, barostats, PME, constraints — is refinement to make this procedure accurate, stable, and physically meaningful.

```
┌─────────────────────────────────────────────────┐
│                    MD Loop                       │
│                                                  │
│  ┌───────────────┐                               │
│  │ Calculate      │                               │
│  │ Forces F(t)    │◄─── Potential energy function │
│  └───────┬───────┘                               │
│          │                                        │
│          ▼                                        │
│  ┌───────────────┐                               │
│  │ Update         │                               │
│  │ velocities     │◄─── F = ma                   │
│  │ v(t+½Δt)      │                               │
│  └───────┬───────┘                               │
│          │                                        │
│          ▼                                        │
│  ┌───────────────┐                               │
│  │ Update         │                               │
│  │ positions      │◄─── r(t+Δt) = r(t) + v·Δt   │
│  │ r(t+Δt)       │                               │
│  └───────┬───────┘                               │
│          │                                        │
│          ▼                                        │
│  ┌───────────────┐                               │
│  │ Apply          │                               │
│  │ constraints,   │                               │
│  │ thermostat,    │                               │
│  │ barostat       │                               │
│  └───────┬───────┘                               │
│          │                                        │
│          ▼                                        │
│     t = t + Δt ──────► Repeat                    │
│                                                  │
└─────────────────────────────────────────────────┘
```

---

## 3. From Quantum Mechanics to Classical Approximation

### 3.1 The Born-Oppenheimer Approximation

In reality, molecules are governed by **quantum mechanics**. The full quantum description requires solving the time-dependent Schrödinger equation for all electrons and nuclei simultaneously — computationally impossible for anything larger than a few atoms.

The **Born-Oppenheimer approximation** (1927) separates the problem into two parts based on a simple physical observation: **electrons are ~1,836× lighter than protons** and therefore move much faster. From the perspective of the electrons, the nuclei are essentially frozen. From the perspective of the nuclei, the electrons instantaneously adjust to any new nuclear configuration.

This gives us a two-step approach:

1. **Solve for the electrons** at fixed nuclear positions → this gives the electronic energy $E_{\text{elec}}(\mathbf{R})$ as a function of nuclear coordinates $\mathbf{R}$.
2. **Treat the nuclei** as moving on the **potential energy surface (PES)** defined by $E_{\text{elec}}(\mathbf{R})$.

The PES is the landscape of energy as a function of all atomic positions. Minima correspond to stable structures; saddle points correspond to transition states; the topology of the PES determines all the thermodynamics and kinetics of the system.

### 3.2 From Electrons to Point Charges

Even with the Born-Oppenheimer approximation, computing $E_{\text{elec}}(\mathbf{R})$ from quantum mechanics is expensive. Ab initio MD (AIMD) does this and is limited to ~100–1,000 atoms for picoseconds.

Classical MD makes a further, much more drastic approximation: **replace the quantum electronic energy surface with an analytical function** — the force field. This means:

- **Electrons disappear entirely.** Their effects are captured implicitly through:
  - Fixed partial charges on atoms (electrostatics)
  - Lennard-Jones parameters (van der Waals / Pauli repulsion)
  - Bonded terms (bonds, angles, dihedrals)
- **Atoms become point particles** obeying Newton's classical mechanics.
- **No bond breaking or formation** — the connectivity (which atoms are bonded) is fixed for the entire simulation.

### 3.3 What We Gain and What We Lose

| We Gain | We Lose |
|---|---|
| Simulations of 10⁵–10⁷ atoms | Quantum effects (tunnelling, zero-point energy) |
| Microsecond–millisecond timescales | Chemical reactions (bond breaking/forming) |
| Direct connection to thermodynamics via statistical mechanics | Electronic excited states, charge transfer |
| Fast force evaluation (~ns per step for 100k atoms) | Explicit electron correlation effects |
| Systematic, reproducible dynamics | Accurate treatment of metal coordination |

> **Key insight for bioinformaticians:** MD cannot simulate enzyme catalysis (bond breaking), electron transfer, or photochemistry. For these, you need QM/MM (quantum mechanics/molecular mechanics) hybrid methods, which treat the active site quantum mechanically and the surroundings classically.

---

## 4. The Potential Energy Function (Force Field)

The total potential energy of the system is:

$$V_{\text{total}} = V_{\text{bonded}} + V_{\text{nonbonded}}$$

> For the full force field comparison (AMBER, CHARMM, OPLS, GROMOS, MARTINI), see the companion tutorial [Understanding GROMACS Key Files & Force Field Models](Understanding_GROMACS_Key_Files.md). Here we focus on the physics of each term.

### 4.1 Bonded Interactions in Detail

#### Bond Stretching (Harmonic Potential)

$$V_{\text{bond}}(r) = \frac{1}{2} k_b (r - r_0)^2$$

- $r$ = current bond length
- $r_0$ = equilibrium bond length (e.g. 0.153 nm for C–C, 0.100 nm for O–H)
- $k_b$ = force constant (stiffness; typically 200,000–500,000 kJ mol⁻¹ nm⁻²)

**Physics:** This is Hooke's law — the same physics as a spring. The harmonic approximation is valid for small deviations from equilibrium (which is almost always the case in MD, because bonds vibrate by only ~0.01 nm around their resting length).

**Limitation:** A harmonic potential cannot model bond dissociation. The energy increases forever as the bond stretches, rather than plateauing (as a real bond does before breaking). This is why classical MD cannot simulate chemical reactions. The Morse potential $V(r) = D_e[1 - e^{-a(r-r_0)}]^2$ is more realistic but rarely used in biomolecular force fields because it is slower to compute and bond breaking is not the goal.

#### Angle Bending

$$V_{\text{angle}}(\theta) = \frac{1}{2} k_\theta (\theta - \theta_0)^2$$

- $\theta$ = current bond angle (formed by three atoms A–B–C)
- $\theta_0$ = equilibrium angle (e.g. 109.5° for sp³ carbon, 120° for sp²)
- $k_\theta$ = force constant (typically 200–800 kJ mol⁻¹ rad⁻²)

**Physics:** Again harmonic. Angle deformations are softer than bond stretching (lower force constants) and contribute to the overall flexibility of the molecule.

#### Proper Dihedral Angles (Torsions)

$$V_{\text{dihedral}}(\phi) = k_\phi [1 + \cos(n\phi - \delta)]$$

- $\phi$ = torsion angle (the angle between planes ABC and BCD for four atoms A–B–C–D)
- $k_\phi$ = barrier height
- $n$ = periodicity (multiplicity: how many minima per full 360° rotation)
- $\delta$ = phase offset

**Physics:** Unlike bonds and angles, dihedral rotations involve relatively **low energy barriers** (typically 1–20 kJ mol⁻¹). These are the primary degrees of freedom that determine protein conformation. The cosine form naturally captures the periodicity of rotation — e.g. $n=3$ for a sp³–sp³ bond (three-fold symmetry with gauche+, trans, gauche⁻ minima).

**Key insight:** The backbone dihedral angles φ (phi: C–N–Cα–C) and ψ (psi: N–Cα–C–N) are the variables plotted on a **Ramachandran plot**. The accuracy of dihedral parameters directly determines how well a force field reproduces secondary structure preferences — this is why force field developers spend enormous effort refitting these terms.

#### Improper Dihedrals

$$V_{\text{improper}}(\xi) = \frac{1}{2} k_\xi (\xi - \xi_0)^2$$

**Physics:** These maintain **planarity** of sp² centres (e.g. the peptide bond, aromatic rings) and **chirality** of chiral centres. Without improper dihedrals, the peptide bond could pucker out of plane and aromatic rings could become non-planar.

### 4.2 Nonbonded Interactions in Detail

Nonbonded interactions are computed for all pairs of atoms that are not directly bonded (or separated by only two bonds). They dominate the computational cost of MD (~90% of CPU time).

#### Lennard-Jones (LJ) Potential

$$V_{\text{LJ}}(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]$$

- $r$ = distance between atoms $i$ and $j$
- $\varepsilon$ = well depth (strength of attraction; typically 0.1–2.0 kJ/mol)
- $\sigma$ = distance at which $V = 0$ (effective atomic "size"; typically 0.25–0.4 nm)

**The physics, term by term:**

| Term | Sign | Physical Origin |
|---|---|---|
| $(\sigma/r)^{12}$ | Repulsive | **Pauli exclusion:** electron clouds cannot overlap. Rises steeply at short range. |
| $-(\sigma/r)^{6}$ | Attractive | **London dispersion forces:** correlated fluctuations in electron clouds create transient dipole-dipole attractions. The $r^{-6}$ dependence is exact from quantum perturbation theory. |

The LJ potential has a minimum at $r_{\min} = 2^{1/6}\sigma \approx 1.122\sigma$, where the energy is $-\varepsilon$.

**Limitation:** The $r^{-12}$ repulsion is not physically rigorous — it is chosen purely for computational convenience (it is the square of $r^{-6}$, allowing fast computation). The real repulsion decays roughly as $e^{-r}$, but the error is small at biomolecular temperatures because few atoms approach the hard repulsive wall.

**In practice:** LJ interactions are short-ranged (significant only within ~1 nm) and are truncated at a cutoff distance $r_c$ (typically 1.0–1.2 nm in GROMACS). A long-range dispersion correction (`DispCorr = EnerPres`) is applied to account for the tail beyond the cutoff.

#### Coulomb (Electrostatic) Interaction

$$V_{\text{Coulomb}}(r) = \frac{q_i q_j}{4\pi\varepsilon_0 r}$$

- $q_i$, $q_j$ = partial atomic charges (in units of elementary charge $e$)
- $\varepsilon_0$ = permittivity of free space
- $r$ = interatomic distance

**Physics:** Electrostatics govern hydrogen bonding, salt bridges, protein-water interactions, and the overall electrostatic environment of the active site. Unlike LJ, the Coulomb interaction decays as $1/r$ — **very slowly** — which means it cannot simply be truncated at a cutoff without serious errors. This is why Particle Mesh Ewald (Section 8) is essential.

**Where the charges come from:** Partial charges are NOT ionic charges. The partial charge on a carbon atom might be +0.02 or −0.18, reflecting the uneven sharing of electron density in covalent bonds. They are determined by:

- **RESP fitting** (AMBER): fit charges to reproduce the quantum mechanical electrostatic potential (ESP) around the molecule.
- **Charge optimisation** (CHARMM): iteratively fit charges to reproduce experimental interaction energies with water and QM dipole moments.

> **Dielectric constant:** In vacuum, $\varepsilon_r = 1$. In classical MD, we do not use a dielectric constant (water's dielectric effect emerges naturally from explicit water molecules). Implicit solvent models use $\varepsilon_r \approx 78$ for water, but this loses the atomistic detail of solvation.

#### 1-4 Interactions

Atoms separated by exactly three bonds (1-4 pairs) receive **scaled** nonbonded interactions rather than full or zero:

| Force Field | 1-4 LJ scaling | 1-4 Coulomb scaling |
|---|---|---|
| AMBER | 0.5 | 0.8333 (= 1/1.2) |
| CHARMM | 1.0 (with NBFIX) | 1.0 |
| OPLS-AA | 0.5 | 0.5 |
| GROMOS | 1.0 | 1.0 |

These scale factors (`fudgeLJ` and `fudgeQQ` in GROMACS) are set in the `[ defaults ]` section of the force field and **must not be changed** — they are an integral part of the parameterisation.

### 4.3 Combining Rules

When computing LJ interactions between atoms of different types, we need to combine their individual $\sigma$ and $\varepsilon$ values. Two schemes are common:

**Lorentz-Berthelot** (AMBER, CHARMM):

$$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2} \qquad \varepsilon_{ij} = \sqrt{\varepsilon_i \cdot \varepsilon_j}$$

**Geometric** (OPLS):

$$\sigma_{ij} = \sqrt{\sigma_i \cdot \sigma_j} \qquad \varepsilon_{ij} = \sqrt{\varepsilon_i \cdot \varepsilon_j}$$

These are defined in the `[ defaults ]` section of the GROMACS topology:
- `comb-rule = 2` → Lorentz-Berthelot
- `comb-rule = 3` → Geometric

> **Warning:** Mixing parameters from force fields that use different combining rules is physically wrong and will give incorrect results.

---

## 5. Newton's Equations of Motion

For each atom $i$ with mass $m_i$ at position $\mathbf{r}_i$:

$$\mathbf{F}_i = m_i \mathbf{a}_i = m_i \frac{d^2\mathbf{r}_i}{dt^2}$$

The force is the negative gradient of the potential:

$$\mathbf{F}_i = -\frac{\partial V}{\partial \mathbf{r}_i}$$

This is a system of $3N$ coupled second-order ordinary differential equations (three spatial dimensions × $N$ atoms). For a typical solvated protein system with 50,000 atoms, that is 150,000 coupled equations. They cannot be solved analytically — we must integrate them numerically.

---

## 6. Numerical Integration: The Leap-Frog Algorithm

### 6.1 Why Not Just Solve the Equations Analytically?

Newton's equations for two interacting particles (the two-body problem) have exact analytical solutions. But for three or more particles (the N-body problem), no general closed-form solution exists. With 50,000+ atoms, numerical integration is the only option.

The key requirements for a good integrator are:

1. **Accuracy:** the trajectory should closely approximate the true dynamics.
2. **Stability:** the total energy should not drift over millions of steps.
3. **Time-reversibility:** running the simulation backward should retrace the same path (a property of Newtonian mechanics).
4. **Symplecticity:** the integrator should preserve the geometric structure of Hamiltonian mechanics (phase-space volume). This guarantees long-term energy conservation.
5. **Efficiency:** the force calculation is expensive, so the integrator should require only **one force evaluation per time step**.

### 6.2 The Verlet Family of Integrators

The standard family of integrators used in MD is based on the **Verlet algorithm** (1967). Start from a Taylor expansion of the position around time $t$:

$$\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t)\Delta t + \frac{1}{2}\mathbf{a}(t)\Delta t^2 + \mathcal{O}(\Delta t^3)$$

$$\mathbf{r}(t - \Delta t) = \mathbf{r}(t) - \mathbf{v}(t)\Delta t + \frac{1}{2}\mathbf{a}(t)\Delta t^2 - \mathcal{O}(\Delta t^3)$$

Adding these two equations:

$$\mathbf{r}(t + \Delta t) = 2\mathbf{r}(t) - \mathbf{r}(t - \Delta t) + \mathbf{a}(t)\Delta t^2 + \mathcal{O}(\Delta t^4)$$

The position error is $\mathcal{O}(\Delta t^4)$ — remarkably accurate. Velocities don't appear explicitly, but can be obtained from:

$$\mathbf{v}(t) = \frac{\mathbf{r}(t + \Delta t) - \mathbf{r}(t - \Delta t)}{2\Delta t}$$

### 6.3 Leap-Frog in Practice

GROMACS uses the **leap-frog** variant, where velocities and positions are evaluated at **alternating half-steps** (they "leap" over each other):

$$\mathbf{v}\left(t + \frac{\Delta t}{2}\right) = \mathbf{v}\left(t - \frac{\Delta t}{2}\right) + \frac{\mathbf{F}(t)}{m}\Delta t$$

$$\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}\left(t + \frac{\Delta t}{2}\right)\Delta t$$

**Step-by-step for one MD iteration:**

```
1. Compute forces F(t) from current positions r(t)
2. Update velocities:  v(t+½Δt) = v(t-½Δt) + [F(t)/m] × Δt
3. Update positions:   r(t+Δt)  = r(t) + v(t+½Δt) × Δt
4. Apply constraints (LINCS/SETTLE)
5. Apply thermostat/barostat
6. t → t + Δt, goto 1
```

**Why leap-frog?**
- Symplectic → excellent long-term energy conservation.
- Time-reversible → satisfies a fundamental property of Newtonian mechanics.
- Only **one force evaluation per time step** → as efficient as possible.
- The velocity at time $t$ (if needed) can be approximated: $\mathbf{v}(t) \approx \frac{1}{2}[\mathbf{v}(t-\frac{\Delta t}{2}) + \mathbf{v}(t+\frac{\Delta t}{2})]$

### 6.4 Time Step Selection

The time step $\Delta t$ must be small enough to resolve the fastest motion in the system. If it is too large, atoms will "overshoot" and the simulation will become unstable (energy will diverge).

| Motion | Period | Frequency |
|---|---|---|
| O–H bond vibration | ~10 fs | ~3,500 cm⁻¹ |
| C–H bond vibration | ~10 fs | ~3,000 cm⁻¹ |
| C–C bond vibration | ~20 fs | ~1,000 cm⁻¹ |
| Bond angle vibration | ~20–50 fs | 500–1,500 cm⁻¹ |
| Dihedral rotation | ~100–10,000 fs | variable |

**Rule of thumb:** $\Delta t$ should be ≤ 1/10 of the shortest oscillation period.

- **Without constraints:** $\Delta t$ = 1 fs (limited by O–H vibrations at ~10 fs)
- **With H-bond constraints (LINCS):** $\Delta t$ = 2 fs (O–H vibrations frozen; next fastest motion is ~20 fs)
- **With all-bond constraints:** $\Delta t$ = 4-5 fs (rarely used; reduces accuracy)
- **With hydrogen mass repartitioning (HMR):** $\Delta t$ = 4 fs (heavier H atoms slow down the fastest vibrations)
- **Coarse-grained (MARTINI):** $\Delta t$ = 20–30 fs (beads are heavier and smoother)

> **Practical impact:** Doubling the time step from 1 fs to 2 fs (by constraining H-bonds) halves the number of force evaluations needed for the same physical time, effectively doubling simulation speed with negligible loss of accuracy. This is why `constraints = h-bonds` is the standard setting.

---

## 7. Constraint Algorithms: LINCS and SETTLE

### 7.1 Why Constrain Bonds?

Bond vibrations (especially X–H) are the fastest motions in a biomolecular system. They are also, from a biological perspective, the least interesting — the *length* of a C–H bond barely affects protein function. By **constraining** these bonds to their equilibrium length, we:

1. Remove the fastest degrees of freedom → allow a larger time step.
2. Are more physically accurate for hydrogen bonds anyway (at 300 K, quantum zero-point energy keeps H bonds near their equilibrium length, but the classical harmonic potential would allow unrealistic stretching).

### 7.2 LINCS (Linear Constraint Solver)

**LINCS** (Hess et al., 1997) is GROMACS's default constraint algorithm for all bonds except water (which uses SETTLE).

**How it works (simplified):**

1. The integrator proposes unconstrained new positions $\mathbf{r}'$.
2. LINCS projects these positions back onto the constraint surface (i.e., corrects bond lengths to their target values).
3. This is done iteratively: first a linear correction, then a non-linear rotation correction.

The key parameter is `lincs_iter` (number of correction iterations, default 1) and `lincs_order` (expansion order, default 4). For most protein simulations, the defaults are fine. Increase `lincs_order` to 6–8 for very stiff systems or if you see LINCS warnings.

**Limitation:** LINCS cannot handle closed loops of constraints (e.g. ring molecules with all bonds constrained). For this, GROMACS falls back to the older SHAKE algorithm.

### 7.3 SETTLE (for Water)

**SETTLE** (Miyamoto & Kollman, 1992) is an analytical constraint algorithm specifically designed for rigid three-site water molecules (e.g. TIP3P, SPC). It solves the constraint equations exactly in one step — no iteration needed.

Since water typically makes up >80% of atoms in a solvated system, SETTLE's speed is critical for overall performance.

---

## 8. Handling Long-Range Electrostatics: Particle Mesh Ewald

### 8.1 The Problem with Coulomb's Law

The Coulomb potential decays as $1/r$. For a finite cutoff $r_c$, the energy of interactions beyond $r_c$ is:

$$\int_{r_c}^{\infty} \frac{1}{r} \cdot 4\pi r^2 \rho \, dr \propto \int_{r_c}^{\infty} r \, dr = \infty$$

In 3D, this integral **diverges**. Simply truncating electrostatics at a cutoff causes:

- Large, discontinuous forces at the cutoff boundary.
- Artificial ordering (molecules align to avoid the cutoff surface).
- Serious errors in structural and thermodynamic properties.

This is unlike Lennard-Jones interactions, where the $r^{-6}$ attraction decays fast enough that truncation is acceptable with a small analytically-computed correction.

**Bottom line:** You *must* handle long-range electrostatics properly. In GROMACS, this means **Particle Mesh Ewald (PME)**.

### 8.2 Ewald Summation — The Idea

The Ewald method (1921, originally for ionic crystals) splits the Coulomb sum into two parts that each converge rapidly:

$$V_{\text{Coulomb}} = \underbrace{V_{\text{short-range}}}_{
\text{Real space}} + \underbrace{V_{\text{long-range}}}_{\text{Reciprocal space}} - \underbrace{V_{\text{self}}}_{\text{Self-correction}}$$

**The trick:** surround each point charge with a Gaussian charge distribution of equal magnitude but opposite sign. This "screens" the charge, making the resulting interaction short-ranged (converges rapidly in real space). Then, to correct for the screening, add back the Gaussian distributions in **reciprocal (Fourier) space**, where they converge rapidly due to the smooth nature of Gaussians.

**Physically:**
1. **Real space sum:** Each charge interacts with nearby charges through a screened (erfc) potential. Computed directly for all pairs within a cutoff. Cost: $\mathcal{O}(N)$.
2. **Reciprocal space sum:** The smooth charge density is represented on a grid and the electrostatic potential is computed using a Fast Fourier Transform. Cost: $\mathcal{O}(N \log N)$.
3. **Self-energy correction:** Remove the artificial interaction of each Gaussian with itself.

### 8.3 Particle Mesh Ewald (PME) — Making It Fast

Classical Ewald summation has a cost of $\mathcal{O}(N^{3/2})$. **PME** (Darden et al., 1993) accelerates the reciprocal-space sum by:

1. **Interpolating** charges onto a regular 3D grid (using B-spline interpolation of order `pme_order`, typically 4).
2. Computing the electrostatic potential on the grid using a **3D FFT** (Fast Fourier Transform).
3. **Interpolating** forces back from the grid to the atoms.

This reduces the reciprocal-space cost to $\mathcal{O}(N \log N)$, making it practical for systems of >100,000 atoms.

```
                    Charge Distribution
                          │
              ┌───────────┴───────────┐
              │                       │
       Screened charges         Smooth Gaussians
       (real space)             (reciprocal space)
              │                       │
       Direct pair sum          Interpolate → Grid
       within cutoff            3D FFT
       O(N)                     Solve Poisson eq.
              │                 Inverse FFT
              │                 Interpolate → Forces
              │                 O(N log N)
              │                       │
              └───────────┬───────────┘
                          │
                   Total Coulomb Force
```

### 8.4 PME Parameters in GROMACS

```
coulombtype    = PME
rcoulomb       = 1.0       ; Real-space cutoff (nm)
pme_order      = 4         ; B-spline interpolation order
fourierspacing = 0.16      ; Grid spacing (nm) — determines grid size
```

| Parameter | Effect of Increasing | Typical Value |
|---|---|---|
| `rcoulomb` | More real-space work, less reciprocal-space work | 1.0–1.2 nm |
| `pme_order` | More accurate interpolation, slightly more expensive | 4 (cubic) |
| `fourierspacing` | Finer grid → more accurate, more expensive FFT | 0.12–0.16 nm |

> **GPU acceleration:** GROMACS can offload PME to a GPU, which is a major performance advantage. The real-space and reciprocal-space calculations can run on the GPU and CPU simultaneously.

---

## 9. Periodic Boundary Conditions

### 9.1 The Minimum Image Convention

In a real experiment, a protein is surrounded by effectively infinite solvent. In a simulation, we can only afford a finite number of water molecules. **Periodic boundary conditions (PBC)** solve this by replicating the simulation box infinitely in all directions:

```
 ┌─────────┬─────────┬─────────┐
 │ image   │ image   │ image   │
 │  (-1,1) │  (0,1)  │  (1,1)  │
 ├─────────┼─────────┼─────────┤
 │ image   │ PRIMARY │ image   │
 │ (-1,0)  │  BOX    │  (1,0)  │
 ├─────────┼─────────┼─────────┤
 │ image   │ image   │ image   │
 │ (-1,-1) │ (0,-1)  │ (1,-1)  │
 └─────────┴─────────┴─────────┘
```

When an atom leaves the box on one side, it re-enters on the opposite side. When computing forces, each atom interacts with the **nearest image** of every other atom (the minimum image convention).

PBC is essential for:
- Eliminating surface effects (no vacuum boundary)
- Enabling the use of Ewald summation (which requires periodicity)
- Maintaining constant density

### 9.2 Box Shapes

| Box Shape | GROMACS `editconf -bt` | Volume Efficiency | Use Case |
|---|---|---|---|
| Cubic | `cubic` | 100% (reference) | Simple, general purpose |
| Rhombic dodecahedron | `dodecahedron` | 71% of cubic | Globular proteins (saves ~29% water) |
| Truncated octahedron | `octahedron` | 77% of cubic | Similar savings to dodecahedron |

**Volume efficiency matters:** Using a dodecahedron box for a roughly spherical protein means ~29% fewer water molecules → ~29% faster simulation.

**Box sizing rule:** The minimum distance from the protein to the box edge should be at least half the nonbonded cutoff, plus a safety margin. Standard practice is **1.0–1.2 nm** from the solute to the nearest box face.

### 9.3 Artefacts and Limitations

**Self-interaction artefact:** If the box is too small, the protein can interact with its own periodic image. This is non-physical. For a cutoff of 1.0 nm, the minimum box vector must be > 2.0 nm, but in practice you need much more (typically the protein diameter + 2–2.4 nm).

**Artificial periodicity:** PBC imposes translational symmetry on the system. This can affect:
- Long-range correlations (the system effectively has crystalline periodicity)
- Free energies of charged species (monopole interactions between periodic images require finite-size corrections)
- Large conformational changes (if the protein unfolds, it may interact with itself across the boundary)

> **For charged systems:** GROMACS requires a net-neutral simulation box for PME. If your protein has a net charge, counterions (Na⁺, Cl⁻) must be added.

---

## 10. Neighbour Searching: The Verlet Cutoff Scheme

Computing all $N(N-1)/2$ pairwise interactions every step would be $\mathcal{O}(N^2)$ — prohibitively expensive. Since LJ interactions are negligible beyond ~1 nm, we can use a **neighbour list**: a list of atom pairs within a cutoff distance, recomputed every `nstlist` steps.

GROMACS uses the **Verlet cutoff scheme:**

1. Build a neighbour list with a **buffer region** (slightly larger than the actual cutoff).
2. Use this list for `nstlist` steps (typically 10–50).
3. The buffer ensures that atoms which drift into the cutoff during those steps are still included.

```
   ┌──────────────────────────────────┐
   │          Buffer zone             │
   │   ┌──────────────────────┐      │
   │   │    Cutoff sphere     │      │
   │   │                      │      │
   │   │     Atom i ●         │      │
   │   │                      │      │
   │   └──────────────────────┘      │
   └──────────────────────────────────┘
   
   r_list = r_cutoff + buffer
```

GROMACS automatically calculates the optimal buffer size based on `nstlist`, the time step, and the target energy drift (<0.005 kJ/mol/ps/atom).

**Spatial decomposition:** GROMACS divides the simulation box into a grid of cells. Only atoms in neighbouring cells need to be checked as potential neighbours. This reduces the search from $\mathcal{O}(N^2)$ to $\mathcal{O}(N)$.

---

## 11. Temperature Control (Thermostats)

### 11.1 Why Control Temperature?

Newton's equations conserve **total energy** (kinetic + potential), generating a trajectory in the **microcanonical (NVE) ensemble**. But experiments are done at constant **temperature** (the canonical, NVT, or isothermal-isobaric, NPT, ensemble).

To simulate at a target temperature $T_0$, we need a thermostat: a mechanism that adds or removes kinetic energy to maintain the desired temperature.

### 11.2 The Equipartition Theorem and Kinetic Energy

The **equipartition theorem** from statistical mechanics states that each quadratic degree of freedom contributes $\frac{1}{2}k_BT$ to the average energy:

$$\langle E_{\text{kinetic}} \rangle = \frac{1}{2}Nf \cdot k_BT$$

where $Nf$ is the number of degrees of freedom (typically $3N - N_c - 3$, with $N_c$ constraints and 3 subtracted for centre-of-mass translation).

The instantaneous temperature is defined from the kinetic energy:

$$T(t) = \frac{2 E_{\text{kin}}(t)}{N_f k_B} = \frac{\sum_i m_i v_i^2}{N_f k_B}$$

The thermostat adjusts velocities so that the **time-averaged** temperature matches $T_0$.

### 11.3 Berendsen Thermostat

$$\frac{dT}{dt} = \frac{T_0 - T}{\tau_T}$$

The Berendsen thermostat (1984) rescales velocities at each step to exponentially relax the temperature toward $T_0$ with a time constant $\tau_T$.

**Scale factor:**

$$\lambda = \sqrt{1 + \frac{\Delta t}{\tau_T}\left(\frac{T_0}{T(t)} - 1\right)}$$

All velocities are multiplied by $\lambda$.

| Advantage | Disadvantage |
|---|---|
| Very efficient at reaching $T_0$ quickly | **Does not generate a correct canonical (NVT) ensemble** |
| Never overshoots | Suppresses fluctuations in kinetic energy |
| Simple and robust | The "flying ice cube" effect: kinetic energy can accumulate in low-frequency modes |

> **Critical limitation:** Berendsen thermostat suppresses energy fluctuations, meaning properties that depend on fluctuations (e.g. heat capacity $C_V = \langle (\delta E)^2 \rangle / k_BT^2$) will be wrong. Use it *only* for equilibration, never for production.

### 11.4 Velocity-Rescale (V-rescale) Thermostat

The V-rescale thermostat (Bussi et al., 2007) is a modified Berendsen thermostat that adds a stochastic term to the rescaling:

$$dK = (K_0 - K)\frac{dt}{\tau_T} + 2\sqrt{\frac{K K_0}{N_f}} \frac{dW}{\sqrt{\tau_T}}$$

where $K$ is kinetic energy, $K_0 = \frac{1}{2}N_f k_B T_0$, and $dW$ is a Wiener noise term.

| Advantage | Disadvantage |
|---|---|
| **Generates a correct canonical ensemble** | Slightly more complex than Berendsen |
| As efficient as Berendsen for equilibration | Requires a stochastic random number generator |
| Correct kinetic energy fluctuations | |
| The **recommended thermostat** in GROMACS | |

### 11.5 Nosé-Hoover Thermostat

The Nosé-Hoover thermostat (Nosé 1984, Hoover 1985) takes a fundamentally different approach: it extends the equations of motion with an additional degree of freedom $\xi$ (a "heat bath" variable) that acts as a friction:

$$\frac{d\mathbf{v}_i}{dt} = \frac{\mathbf{F}_i}{m_i} - \xi(t) \mathbf{v}_i$$

$$\frac{d\xi}{dt} = \frac{1}{Q}\left[\sum_i m_i v_i^2 - N_f k_B T_0\right]$$

$Q$ is the "mass" of the thermostat (related to $\tau_T$). The friction $\xi$ increases when the system is too hot and decreases when too cold.

| Advantage | Disadvantage |
|---|---|
| **Generates a rigorous canonical ensemble** | Can exhibit oscillatory relaxation (temperature oscillates before converging) |
| Deterministic (no random numbers) | Not efficient for equilibration (slow convergence) |
| Well-suited for production MD | Can couple to resonant motions in small systems |
| Time-reversible | Requires a chain of thermostats (Nosé-Hoover chains) for ergodicity in some cases |

### 11.6 Which Thermostat Should I Use?

| Stage | Recommended Thermostat | Why |
|---|---|---|
| **Energy minimisation** | None | Not an MD simulation |
| **Equilibration (NVT, NPT)** | **V-rescale** | Fast convergence to target T; correct ensemble |
| **Production MD** | **V-rescale** or **Nosé-Hoover** | Both give correct NVT ensemble |

> **Practical recommendation:** V-rescale (`tcoupl = V-rescale`) is the standard choice in modern GROMACS workflows for both equilibration and production. It combines Berendsen's fast relaxation with correct canonical sampling.

---

## 12. Pressure Control (Barostats)

### 12.1 Why Control Pressure?

Experiments are typically performed at constant (atmospheric) pressure, not constant volume. To reproduce this, we simulate in the **NPT ensemble** (constant Number of particles, Pressure, and Temperature). The barostat adjusts the simulation **box volume** to maintain the target pressure.

### 12.2 The Virial and Pressure

The instantaneous pressure in a simulation is computed from the **virial**:

$$P = \frac{1}{V}\left[Nk_BT + \frac{1}{3}\sum_{i<j} \mathbf{r}_{ij} \cdot \mathbf{F}_{ij}\right]$$

The first term is the ideal gas pressure (from kinetic energy). The second term is the virial — the contribution from interatomic forces. For condensed phases, the virial term is large and negative (attractive forces), and the two terms nearly cancel, making the **instantaneous pressure extremely noisy**. Fluctuations of ±100–500 bar are normal and expected in a biomolecular simulation at 1 bar.

### 12.3 Berendsen Barostat

Analogous to the Berendsen thermostat, the Berendsen barostat rescales the box vectors and atomic coordinates:

$$\mu = \left[1 - \frac{\beta \Delta t}{\tau_p}(P_0 - P(t))\right]^{1/3}$$

All positions and box vectors are scaled by $\mu$.

| Advantage | Disadvantage |
|---|---|
| Very efficient at converging the density | **Does not generate a correct NPT ensemble** |
| No oscillations | Suppresses volume fluctuations |
| Robust for equilibration | Volume-dependent properties will be wrong |

> **Use Berendsen only for equilibration**, never for production MD or any analysis that depends on volume fluctuations (compressibility, thermal expansion, etc.).

### 12.4 Parrinello-Rahman Barostat

The Parrinello-Rahman barostat (1981) treats the box vectors as dynamical variables with their own equations of motion:

$$\frac{d^2\mathbf{h}}{dt^2} = V W^{-1}(\mathbf{P} - P_0 \mathbf{I})$$

where $\mathbf{h}$ is the box matrix, $W$ is the "mass" of the box (related to $\tau_p$ and the compressibility), $\mathbf{P}$ is the pressure tensor, and $P_0$ is the target pressure.

| Advantage | Disadvantage |
|---|---|
| **Generates a correct NPT ensemble** | Can exhibit box oscillations if started far from equilibrium |
| Allows anisotropic box deformation | Should not be used from a non-equilibrated starting point |
| Widely used and validated | Requires a well-equilibrated density (use Berendsen first) |

### 12.5 Which Barostat Should I Use?

| Stage | Recommended Barostat | Why |
|---|---|---|
| **NPT Equilibration** | **Berendsen** | Fast density convergence without oscillation |
| **Production MD** | **Parrinello-Rahman** | Correct NPT ensemble |
| **Free energy calculations** | **Parrinello-Rahman** or **MTTK** | Must have correct volume fluctuations |

---

## 13. Statistical Mechanics: Connecting Simulation to Experiment

### 13.1 Ensembles

An **ensemble** is a collection of all possible microstates (atomic configurations) consistent with a set of macroscopic constraints. MD generates a time series of microstates, and we extract observable properties by averaging over this trajectory.

| Ensemble | Fixed Quantities | Fluctuating Quantities | GROMACS Realisation |
|---|---|---|---|
| **NVE** (Microcanonical) | N, V, E | T, P | No thermostat, no barostat |
| **NVT** (Canonical) | N, V, T | E, P | Thermostat only |
| **NPT** (Isothermal-isobaric) | N, P, T | E, V | Thermostat + barostat |

Most biomolecular MD is done in NPT, matching experimental conditions (constant temperature and pressure, with the system free to expand or contract).

### 13.2 The Ergodic Hypothesis

The **ergodic hypothesis** states that, given sufficient time, a system will visit all accessible microstates, and the time average of any property will equal the ensemble average:

$$\langle A \rangle_{\text{time}} = \lim_{T \to \infty} \frac{1}{T}\int_0^T A(t) \, dt \approx \langle A \rangle_{\text{ensemble}}$$

This is the fundamental assumption that connects MD simulation (which generates a time series) to statistical mechanics (which works with ensembles).

**In practice, ergodicity often fails:** If there are high energy barriers between conformational states, the simulation may be "trapped" in one state for the entire simulation time. This is the **sampling problem** — the most fundamental challenge in MD.

**Examples of sampling failure:**
- A protein that folds on a millisecond timescale cannot be fully sampled in a microsecond simulation.
- A ligand may not dissociate/reassociate within the simulation time.
- Side chain rotamer transitions in buried residues may be too slow to sample.

**Enhanced sampling methods** (replica exchange / REMD, metadynamics, accelerated MD, steered MD) are designed to overcome these barriers by biasing the simulation to explore more of the energy landscape.

### 13.3 Free Energy and Sampling

The **Helmholtz free energy** (NVT) and **Gibbs free energy** (NPT) cannot be directly computed from a single MD trajectory because they depend on the *partition function* (a sum over **all** microstates, not just those visited):

$$G = -k_BT \ln Z$$

However, **free energy *differences*** between two states can be computed using:

- **Free energy perturbation (FEP):** $\Delta G = -k_BT \ln \langle e^{-\Delta V/k_BT} \rangle_A$
- **Thermodynamic integration (TI):** $\Delta G = \int_0^1 \langle \frac{\partial V(\lambda)}{\partial \lambda} \rangle_\lambda d\lambda$
- **Bennett acceptance ratio (BAR):** optimal combination of forward and reverse perturbation data.

These methods are used in GROMACS for computing binding free energies, solvation free energies, and the effect of mutations.

---

## 14. Energy Minimisation

Before running MD, the starting structure (from X-ray/cryo-EM/homology model) typically has steric clashes, unrealistic bond lengths, or bad contacts with solvent. Energy minimisation (EM) relaxes these.

### 14.1 Steepest Descent

The simplest optimisation algorithm. At each step, move each atom in the direction of the force (negative gradient):

$$\mathbf{r}_{i}^{(n+1)} = \mathbf{r}_{i}^{(n)} + \frac{\mathbf{F}_i^{(n)}}{|\mathbf{F}_{\max}^{(n)}|} h_n$$

where $h_n$ is the step size. If the energy decreases, the step size is increased; if it increases, the step is rejected and $h_n$ is reduced.

| Advantage | Disadvantage |
|---|---|
| Extremely robust — always converges | Slow near the minimum (zig-zag behaviour) |
| Handles very bad starting structures | Cannot escape local minima |
| Simple implementation | Not suitable for precise optimisation |

### 14.2 Conjugate Gradient

Uses information from the previous step to choose a better search direction (conjugate to the previous gradient). Converges much faster than steepest descent near a minimum, but is **less robust** for very bad starting structures.

**Typical protocol:** Run steepest descent to remove the worst clashes, then switch to conjugate gradient for further refinement.

### 14.3 Why Minimise Before MD?

1. **Remove steric clashes** that would cause enormous forces and crash the simulation.
2. **Relax added hydrogen atoms** (which are placed geometrically by `pdb2gmx` and may overlap).
3. **Relax solvent** around the protein.
4. Bring the system to a nearby local energy minimum so that the initial velocities (assigned from a Maxwell-Boltzmann distribution) don't cause instability.

> **Important:** Energy minimisation finds a **local minimum**, not the global minimum. The energy landscape of a protein has an astronomically large number of local minima (the "Levinthal paradox"). EM simply removes initial bad contacts — it does not "fold" the protein or find the native state.

---

## 15. Equilibration: NVT and NPT

After energy minimisation, the system must be **equilibrated** — brought to thermal equilibrium at the target temperature and pressure before production data can be collected.

**Standard two-phase protocol:**

### Phase 1: NVT Equilibration (Heating)

```
integrator = md
dt         = 0.002
nsteps     = 50000    ; 100 ps
tcoupl     = V-rescale
ref_t      = 300
pcoupl     = no       ; No pressure coupling yet
define     = -DPOSRES ; Position restraints on protein
```

**Purpose:** Heat the system from 0 K (or wherever the velocities start) to 300 K, while restraining the protein. The solvent molecules rearrange around the restrained protein. Monitor the temperature — it should stabilise at 300 K.

### Phase 2: NPT Equilibration (Density relaxation)

```
integrator   = md
dt           = 0.002
nsteps       = 50000    ; 100 ps
tcoupl       = V-rescale
ref_t        = 300
pcoupl       = Berendsen    ; Fast density convergence
ref_p        = 1.0
define       = -DPOSRES     ; Position restraints still on
```

**Purpose:** Now allow the box volume to adjust to achieve the correct density at 1 bar. The density should converge to ~1,000 kg/m³ (for a water-solvated protein). Position restraints remain on the protein.

### Why Position Restraints During Equilibration?

The crystal structure was determined at cryogenic temperature or in a crystal lattice. Throwing it into a box of water at 300 K and letting it move freely would cause:

- Rapid, unrealistic conformational changes as the protein adjusts to solvent.
- Potential unfolding of marginally stable regions.
- Solvent not yet optimally packed around the protein.

Position restraints allow the solvent to relax while the protein stays near its starting conformation. After equilibration, restraints are removed for production MD.

---

## 16. Key Assumptions of Molecular Dynamics

Every MD simulation rests on a stack of assumptions. Understanding these is essential for interpreting results and knowing when they might fail.

| # | Assumption | What It Means | When It Might Fail |
|---|---|---|---|
| **1** | **Born-Oppenheimer approximation** | Electrons adjust instantaneously to nuclear motion | Near conical intersections, in photochemistry, or with very fast proton transfer |
| **2** | **Classical nuclei** | Atoms obey Newton's laws, not quantum mechanics | For light atoms (H) at low temperatures — quantum tunnelling and zero-point energy are ignored |
| **3** | **Fixed bonding topology** | Bonds cannot break or form | Any chemical reaction, proton transfer, or disulfide bond formation/breakage |
| **4** | **Point charges and fixed partial charges** | No electronic polarisation | Environments with very different dielectric (e.g. membrane interior vs. water), ion binding, metal coordination |
| **5** | **Pairwise additive potential** | Total energy is a sum of pair interactions | Many-body effects (polarisation, charge transfer) are averaged out; fails for metal clusters, polarisable media |
| **6** | **Harmonic bonds and angles** | Small deviations from equilibrium only | Near bond dissociation (never happens in standard MD, but means the model is inherently local) |
| **7** | **Effective pair LJ parameters** | The $r^{-12}$ repulsion is empirical | Works well in practice but is not physically rigorous |
| **8** | **Periodic boundary conditions** | The system is infinite and periodic | Finite-size effects for charged systems, long-range correlations, very large conformational changes |
| **9** | **Ergodic sampling** | The trajectory visits all relevant states | Fails for processes slower than the simulation timescale (folding, slow binding) |
| **10** | **Classical statistical mechanics** | Energy is continuous and freely exchangeable | Quantum effects in vibrational modes (especially X–H stretching); heat capacity is overestimated by classical MD |

---

## 17. Known Limitations and When MD Can Fail

### 17.1 The Timescale Problem

The most fundamental limitation. The 2 fs time step (set by bond vibrations) means that even a microsecond simulation requires 500 million steps. Many biological processes happen on timescales that are currently inaccessible:

| Process | Typical Timescale | Accessible by MD? (2026) |
|---|---|---|
| Bond vibration | 10 fs | ✅ Trivially |
| Water reorientation | 1–10 ps | ✅ |
| Side chain rotation | 10 ps – 10 ns | ✅ |
| Loop motion | 1 ns – 1 μs | ✅ (mostly) |
| Domain motion | 100 ns – 1 ms | ⚠️ Borderline |
| Protein folding | 1 μs – 10 s | ❌ For most proteins |
| Ligand binding (kon) | 1 μs – 1 ms | ⚠️ With enhanced sampling |
| Large conformational change | 1 μs – 1 s | ❌ Usually |

**Mitigation:** Enhanced sampling methods (REMD, metadynamics, accelerated MD, Gaussian accelerated MD), longer trajectories on GPU clusters, or coarse-grained models.

### 17.2 Force Field Accuracy

The force field is an approximation. Known issues include:

| Problem | Details |
|---|---|
| **Secondary structure bias** | Older force fields (e.g. AMBER ff99) over-stabilise α-helices. Modern FFs (ff19SB, CHARMM36m) have improved but are not perfect. |
| **Intrinsically disordered proteins** | IDPs are exquisitely sensitive to the balance between protein–protein and protein–water interactions. Many FFs make IDPs too compact. ff19SB + OPC and CHARMM36m + TIP4P-D are recent improvements. |
| **Salt bridges** | Some FFs overestimate the stability of salt bridges in solution; ion parameters remain an active area of research. |
| **Protein–RNA/DNA interfaces** | Cross-parameterisation between protein and nucleic acid FFs is not always well validated. |
| **Membrane systems** | Lipid force fields (CHARMM36 lipids, Slipids, Lipid17) differ in their predictions of membrane properties; area per lipid, order parameters, and phase behaviour vary. |
| **Metals and cofactors** | Classical FFs struggle with metal coordination (electronic effects). Specialised approaches like the bonded model, nonbonded model, or the cationic dummy atom model are needed for metalloenzymes. |
| **Transferability** | Parameters fitted to one type of system may not work well for another (e.g. parameters for proteins in water may not be accurate for proteins at an interface). |

### 17.3 Finite-Size Effects

With periodic boundary conditions:

- **Charged solutes** interact with their periodic images, leading to finite-size errors in solvation free energies. Corrections (e.g. the Rocklin correction) exist but add complexity.
- **Long-range order** is artificially imposed — the simulation behaves as if the system is a crystal of repeating units.
- **Conformational changes** that make the protein larger than the box will cause self-interaction artefacts.

### 17.4 Water Model Limitations

| Issue | Details |
|---|---|
| **TIP3P** | Too fast dynamics (diffusion ~2× too high), low viscosity. Structural properties are reasonable. |
| **No quantum effects** | Classical water cannot reproduce the anomalous properties of water (density maximum at 4°C) with a simple model. |
| **No polarisation** | Fixed-charge water models cannot respond to the local electrostatic environment. This matters near charged surfaces and ions. |
| **Phase behaviour** | Most 3-site models do not correctly predict ice phases, boiling points, or the temperature of maximum density. TIP4P/2005 is the best for phase diagram. |

### 17.5 Treatment of Electrostatics

| Issue | Details |
|---|---|
| **PME requires periodicity** | Cannot be used for non-periodic systems (gas phase, droplets) without modification. |
| **Tinfoil boundary conditions** | Standard Ewald/PME uses conducting (tinfoil) boundary conditions at infinity, which may not represent the physical boundary. |
| **Charge neutrality** | PME requires net-zero charge. A uniform background charge is implicitly applied, which introduces artefacts for charged systems. |

### 17.6 Sampling and Convergence

**How do you know your simulation is converged?** This is one of the hardest questions in MD. Common diagnostics:

- **RMSD plateau:** The RMSD from the starting structure should level off.
- **Property convergence:** The running average of observables (radius of gyration, SASA, number of H-bonds) should stabilise.
- **Block averaging:** Divide the trajectory into blocks and check that block averages are consistent.
- **Multiple independent simulations:** The gold standard. Run 3–5 independent simulations from different starting velocities and check that they converge to the same distribution.

> **Golden rule:** A single short MD trajectory is an anecdote, not a statistic. Always run replicas.

### 17.7 Reproducibility and Chaos

MD trajectories are **chaotic** in the mathematical sense: tiny perturbations in initial conditions (e.g. the last decimal place of one coordinate) lead to completely different trajectories after a few picoseconds. This is **not a bug** — it is a fundamental property of classical many-body dynamics (Lyapunov instability).

**Implications:**
- Two "identical" simulations with different random seeds will diverge rapidly.
- This is fine — what matters is that the **ensemble properties** (averages and distributions) are the same.
- Single-trajectory observations ("the loop moved at 50 ns") are not reproducible and should not be over-interpreted.

---

## 18. Summary Table of Key Algorithms in GROMACS

| Component | Algorithm | GROMACS Setting | Purpose |
|---|---|---|---|
| **Integration** | Leap-frog | `integrator = md` | Advance positions and velocities |
| **Constraints (bonds)** | LINCS | `constraint_algorithm = lincs` | Freeze bond lengths → 2 fs time step |
| **Constraints (water)** | SETTLE | Automatic for rigid water | Rigid water geometry |
| **Electrostatics** | PME | `coulombtype = PME` | Accurate long-range Coulomb |
| **Neighbour search** | Verlet | `cutoff-scheme = Verlet` | Efficient pair list with buffer |
| **Thermostat** | V-rescale | `tcoupl = V-rescale` | Correct NVT ensemble |
| **Thermostat** | Nosé-Hoover | `tcoupl = nose-hoover` | Alternative correct NVT |
| **Barostat (equil.)** | Berendsen | `pcoupl = Berendsen` | Fast density convergence |
| **Barostat (prod.)** | Parrinello-Rahman | `pcoupl = Parrinello-Rahman` | Correct NPT ensemble |
| **Minimisation** | Steepest descent | `integrator = steep` | Remove bad contacts |
| **Minimisation** | Conjugate gradient | `integrator = cg` | Refine near minimum |

---

## 19. Glossary

| Term | Definition |
|---|---|
| **Å (Ångström)** | Unit of length: 1 Å = 0.1 nm = 10⁻¹⁰ m |
| **B-factor** | Crystallographic temperature factor; reflects atomic displacement/flexibility |
| **Barostat** | Algorithm that controls pressure by adjusting box volume |
| **Born-Oppenheimer** | Approximation separating nuclear and electronic motion |
| **Cutoff** | Distance beyond which nonbonded interactions are truncated or switched |
| **Dihedral** | Torsion angle defined by four atoms |
| **Ensemble** | Statistical mechanical collection of microstates under fixed macroscopic conditions |
| **Ergodic** | Property of a system that visits all microstates given sufficient time |
| **Force field** | Mathematical model + parameter set defining interatomic interactions |
| **Lennard-Jones** | Model potential for van der Waals interactions ($r^{-12}$ repulsion + $r^{-6}$ attraction) |
| **LINCS** | Linear Constraint Solver for bond length constraints |
| **Leap-frog** | Symplectic integration algorithm used by GROMACS |
| **nm (nanometre)** | Unit of length: 1 nm = 10 Å. GROMACS internal unit. |
| **PBC** | Periodic Boundary Conditions |
| **PES** | Potential Energy Surface — energy as a function of all atomic coordinates |
| **PME** | Particle Mesh Ewald — algorithm for long-range electrostatics |
| **RMSD** | Root Mean Square Deviation of atomic positions from a reference |
| **Symplectic** | Property of integrators that conserve phase-space volume (ensures energy stability) |
| **Thermostat** | Algorithm that controls temperature by adjusting kinetic energy |
| **Virial** | Contribution of interatomic forces to pressure |

---

## 20. Further Reading & Next Steps

### Within This Tutorial Series

- [Understanding Protein Files](Understanding_Protein_Files.md) — PDB/mmCIF format details and Python manipulation.
- [Understanding GROMACS Key Files & Force Field Models](Understanding_GROMACS_Key_Files.md) — File format reference and force field comparison.
- [MD on Apocrita](../MD_On_Apocrita.md) — Practical step-by-step guide to running GROMACS on the QMUL HPC cluster.

### Textbooks

- **Frenkel & Smit**, *Understanding Molecular Simulation* (Academic Press) — The definitive textbook on the statistical mechanics and algorithms of simulation.
- **Leach**, *Molecular Modelling: Principles and Applications* (Pearson) — Excellent broad introduction for bioscientists.
- **Allen & Tildesley**, *Computer Simulation of Liquids* (Oxford) — Classic reference for simulation methodology.

### Key Papers

- Verlet (1967) — Original Verlet integrator. *Phys. Rev.* 159, 98.
- Hess et al. (1997) — LINCS algorithm. *J. Comput. Chem.* 18, 1463.
- Darden et al. (1993) — PME. *J. Chem. Phys.* 98, 10089.
- Bussi et al. (2007) — V-rescale thermostat. *J. Chem. Phys.* 126, 014101.
- Parrinello & Rahman (1981) — PR barostat. *J. Appl. Phys.* 52, 7182.

### GROMACS Documentation

- [GROMACS Reference Manual](https://manual.gromacs.org/) — Comprehensive documentation of all algorithms, options, and file formats.
- [GROMACS Tutorials by Justin Lemkul](http://www.mdtutorials.com/) — Excellent practical walkthroughs.

---

*Tutorial written for the QMUL Molecular Dynamics Teaching Project. Last updated: March 2026.*