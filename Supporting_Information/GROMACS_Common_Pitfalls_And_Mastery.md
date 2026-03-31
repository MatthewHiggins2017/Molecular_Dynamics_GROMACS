# GROMACS Common Pitfalls and Mastery
**Avoiding the Traps and Building Deep Expertise in Molecular Dynamics Simulation**

---

## Table of Contents

**Part A — Common Pitfalls**

1. [System Preparation Pitfalls](#1-system-preparation-pitfalls)
   - 1.1 [Unit Confusion: Ångströms vs. Nanometres](#11-unit-confusion-ångströms-vs-nanometres)
   - 1.2 [Missing or Incorrect Protonation States](#12-missing-or-incorrect-protonation-states)
   - 1.3 [Missing Residues and Incomplete Structures](#13-missing-residues-and-incomplete-structures)
   - 1.4 [Forgetting Disulfide Bonds](#14-forgetting-disulfide-bonds)
   - 1.5 [Wrong Force Field for Your System](#15-wrong-force-field-for-your-system)
   - 1.6 [Mismatched Water Model](#16-mismatched-water-model)
   - 1.7 [Box Too Small](#17-box-too-small)
   - 1.8 [Net Charge Not Neutralised](#18-net-charge-not-neutralised)
   - 1.9 [Topology `[ molecules ]` Mismatch](#19-topology--molecules--mismatch)
2. [Energy Minimisation Pitfalls](#2-energy-minimisation-pitfalls)
   - 2.1 [Skipping Minimisation Entirely](#21-skipping-minimisation-entirely)
   - 2.2 [Maximum Force Not Converging](#22-maximum-force-not-converging)
   - 2.3 [Minimisation "Succeeds" but System Is Wrong](#23-minimisation-succeeds-but-system-is-wrong)
3. [Equilibration Pitfalls](#3-equilibration-pitfalls)
   - 3.1 [Skipping or Rushing Equilibration](#31-skipping-or-rushing-equilibration)
   - 3.2 [Using Parrinello-Rahman from the Start](#32-using-parrinello-rahman-from-the-start)
   - 3.3 [Forgetting Position Restraints](#33-forgetting-position-restraints)
   - 3.4 [Temperature Not Stabilising](#34-temperature-not-stabilising)
   - 3.5 [Density Not Converging in NPT](#35-density-not-converging-in-npt)
4. [Production MD Pitfalls](#4-production-md-pitfalls)
   - 4.1 [Berendsen Thermostat/Barostat in Production](#41-berendsen-thermostatbarostat-in-production)
   - 4.2 [LINCS Warnings and Crashes](#42-lincs-warnings-and-crashes)
   - 4.3 [Atom Out of Domain / Coordinate Reset Errors](#43-atom-out-of-domain--coordinate-reset-errors)
   - 4.4 [Energy Drift in Long Simulations](#44-energy-drift-in-long-simulations)
   - 4.5 [Insufficient Output Frequency](#45-insufficient-output-frequency)
   - 4.6 [Wall-Clock Timeouts on HPC Clusters](#46-wall-clock-timeouts-on-hpc-clusters)
5. [Analysis Pitfalls](#5-analysis-pitfalls)
   - 5.1 [Analysing Non-Equilibrated Data](#51-analysing-non-equilibrated-data)
   - 5.2 [PBC Artefacts in Trajectories](#52-pbc-artefacts-in-trajectories)
   - 5.3 [Fitting and Alignment Errors](#53-fitting-and-alignment-errors)
   - 5.4 [Drawing Conclusions from a Single Trajectory](#54-drawing-conclusions-from-a-single-trajectory)
   - 5.5 [RMSD Misinterpretation](#55-rmsd-misinterpretation)
   - 5.6 [Misusing Hydrogen Bond Cutoffs](#56-misusing-hydrogen-bond-cutoffs)

**Part B — Mastering GROMACS**

6. [Level 1: Competent User](#6-level-1-competent-user)
7. [Level 2: Proficient Practitioner](#7-level-2-proficient-practitioner)
8. [Level 3: Expert / Power User](#8-level-3-expert--power-user)
9. [Level 4: Developer-Level Mastery](#9-level-4-developer-level-mastery)
10. [Essential Command-Line Mastery](#10-essential-command-line-mastery)
11. [Building a Robust Simulation Workflow](#11-building-a-robust-simulation-workflow)
12. [Validation Checklist: The "Pre-Flight Check"](#12-validation-checklist-the-pre-flight-check)
13. [Performance Tuning](#13-performance-tuning)
14. [Troubleshooting Decision Tree](#14-troubleshooting-decision-tree)
15. [Resources for Continued Learning](#15-resources-for-continued-learning)

---

# Part A — Common Pitfalls

---

## 1. System Preparation Pitfalls

System preparation is where the vast majority of simulation errors originate. A flawed starting system produces flawed results — and GROMACS will happily run a garbage simulation to completion without complaint.

### 1.1 Unit Confusion: Ångströms vs. Nanometres

**The mistake:** PDB files use **Ångströms** (Å). GROMACS uses **nanometres** (nm) internally. New users frequently confuse the two, especially when manually editing `.gro` files or setting cutoffs.

**Consequences:**
- Setting `rcoulomb = 10.0` when you meant 1.0 nm (10 Å) → electrostatics computed over 10 nm, catastrophically slow.
- Editing a `.gro` file with Å coordinates → the protein is 10× too large and the simulation immediately crashes.

**How to avoid it:**

```
1 nm = 10 Å

PDB file:   27.340 Å  →  GROMACS .gro: 2.734 nm
Cutoff:     1.0 nm  =  10 Å  (standard setting)
Box size:   5.0 nm  =  50 Å
```

> **Rule:** If a number in a `.gro` file looks like it should be in a PDB file, it's probably in the wrong units. `.gro` coordinates for a typical protein are in the range 0–10 nm, not 0–100.

### 1.2 Missing or Incorrect Protonation States

**The mistake:** Using the PDB file directly without considering the protonation state of titratable residues at your simulation pH.

**The problem:** PDB files from X-ray crystallography **rarely include hydrogen atoms** and never explicitly tell you the protonation state. At physiological pH (7.4):

| Residue | pKa | State at pH 7.4 | Easily Overlooked? |
|---|---|---|---|
| ASP | 3.7 | Deprotonated (ASP⁻) | No |
| GLU | 4.1 | Deprotonated (GLU⁻) | No |
| **HIS** | **6.0** | **Ambiguous** — could be neutral (HIE/HID) or charged (HIP) | **YES** |
| CYS | 8.3 | Usually protonated (neutral) unless in disulfide | Sometimes |
| LYS | 10.5 | Protonated (LYS⁺) | No |
| TYR | 10.1 | Protonated (neutral) | No |
| N-terminus | 8.0 | Protonated (NH₃⁺) | No |
| C-terminus | 3.1 | Deprotonated (COO⁻) | No |

**Histidine is the most dangerous:** At pH 7.4, histidine can be:
- **HIE (HIS-ε):** Proton on Nε — the most common single-proton tautomer.
- **HID (HIS-δ):** Proton on Nδ — less common but occurs in specific environments.
- **HIP (HIS⁺):** Both protons present, positively charged — common in active sites and at lower pH.

`pdb2gmx` will ask you to choose for each histidine — **never accept the default blindly**.

**How to avoid it:**
1. Use **H++ (http://newbiophysics.cs.vt.edu/H++)** or **PropKa** to predict pKa values for your specific structure.
2. Inspect the local environment of each histidine manually — is it coordinating a metal ion? Near a charged residue? Buried or solvent-exposed?
3. Document your protonation state choices.

### 1.3 Missing Residues and Incomplete Structures

**The mistake:** Using a PDB structure with missing loops, terminal residues, or unresolved regions without modelling them.

**The problem:** X-ray and cryo-EM structures frequently have:
- **Missing loops:** Flexible regions that are disordered in the crystal/particle.
- **Missing terminal residues:** N- and C-termini are often disordered.
- **Missing side chains:** Surface residues with poor electron density.

If you feed a structure with missing residues to `pdb2gmx`, it will either crash or generate a disconnected topology (a chain break where there should be a peptide bond).

**How to avoid it:**
1. Check for `REMARK 465` (missing residues) and `REMARK 470` (missing atoms) in the PDB file.
2. Model missing loops with **MODELLER**, **SWISS-MODEL**, or **AlphaFold**.
3. Use `pdb2gmx -missing` to allow missing residues (only if you intentionally want chain breaks).

### 1.4 Forgetting Disulfide Bonds

**The mistake:** Running `pdb2gmx` without specifying or verifying disulfide bonds.

**The problem:** `pdb2gmx` detects disulfide bonds based on the Sγ–Sγ distance (< 0.3 nm), but:
- If your structure has imperfect geometry, it may miss a real disulfide.
- If you modified the structure (e.g. after loop modelling), cysteine positions may have shifted.

Missing disulfides will cause the protein to unfold in regions that should be constrained.

**How to avoid it:**
1. Check the PDB `SSBOND` records.
2. Verify that `pdb2gmx` reports the expected disulfides in its output.
3. Use `-ss yes` to force automatic detection, or `-ss no` and provide them manually.

### 1.5 Wrong Force Field for Your System

**The mistake:** Choosing a force field without considering what your system contains.

**Common errors:**
| Scenario | Wrong Choice | Better Choice |
|---|---|---|
| Protein in membrane | OPLS-AA (no mature lipid FF) | CHARMM36m (integrated lipid parameters) |
| Protein + small molecule ligand | Any FF without ligand parameters | CHARMM36m + CGenFF or AMBER + GAFF2 |
| Intrinsically disordered protein | AMBER ff99SB (helix bias) | AMBER ff19SB + OPC |
| Very large system (>1M atoms) | Any all-atom FF | MARTINI 3 (coarse-grained) |

**How to avoid it:** See the [force field decision tree](Understanding_GROMACS_Key_Files.md#18-how-to-choose-a-force-field) in the companion tutorial. Make the decision deliberately and document it.

### 1.6 Mismatched Water Model

**The mistake:** Using a water model that was not parameterised with your force field.

**Why it matters:** Force field partial charges and Lennard-Jones parameters are **fitted against a specific water model**. Using TIP4P water with AMBER ff14SB (parameterised with TIP3P) will produce incorrect solvation energies, protein–water interactions, and structural properties.

| Force Field | Correct Water Model |
|---|---|
| AMBER ff14SB | TIP3P |
| AMBER ff19SB | OPC (or TIP3P) |
| CHARMM36m | CHARMM-modified TIP3P (TIPS3P) |
| OPLS-AA/M | TIP3P or TIP4P |
| GROMOS 54A7 | SPC or SPC/E |

**How to avoid it:** `pdb2gmx` will list the available water models for your chosen force field — pick one from that list.

### 1.7 Box Too Small

**The mistake:** Making the simulation box barely fit around the protein.

**The problem:** If any atom can "see" its own periodic image within the nonbonded cutoff, you get self-interaction. This corrupts forces, energies, and structural properties. GROMACS will print a warning but **will not stop the simulation**.

**Typical warning:**
```
WARNING: Atoms are too close to the box boundaries.
         This likely means your simulation box is too small.
```

**How to avoid it:**
```bash
# Ensure at least 1.0 nm between solute and box edge
gmx editconf -f protein.gro -o boxed.gro -c -d 1.0 -bt dodecahedron
```

**Rule of thumb:** The minimum distance from any solute atom to the box wall should be ≥ `rcoulomb` (typically 1.0 nm). Using **1.2 nm** provides a safety margin. For proteins that undergo large conformational changes, use **1.5 nm or more**.

### 1.8 Net Charge Not Neutralised

**The mistake:** Running a simulation with a net charge on the system.

**The problem:** PME (Particle Mesh Ewald) requires an electrostatically neutral simulation box. A net charge causes:
- A divergent Ewald energy (corrected internally by a uniform background charge, but this introduces artefacts).
- Incorrect electrostatic forces, especially near the solute.

**How to avoid it:**
```bash
# Add counter-ions to neutralise
gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral

# Or add a specific salt concentration (e.g. 0.15 M NaCl)
gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```

Always check the total charge reported by `grompp` — it should be exactly 0 (or very close due to rounding).

### 1.9 Topology `[ molecules ]` Mismatch

**The mistake:** The `[ molecules ]` section of the `.top` file does not match the coordinate file.

**When this happens:**
- You ran `gmx solvate` but forgot to update the water count.
- You ran `gmx genion` (which replaces some water molecules with ions) but didn't update both the water and ion counts.
- You manually edited the coordinate file without updating the topology.

**The error:**
```
Fatal error:
number of coordinates in coordinate file does not match topology
```

**How to avoid it:**
- `gmx solvate` updates the topology automatically (if you pass `-p topol.top`).
- `gmx genion` also updates it, but **always double-check** the `[ molecules ]` section after every step.
- The order in `[ molecules ]` must match the order of atoms in the `.gro` file exactly.

---

## 2. Energy Minimisation Pitfalls

### 2.1 Skipping Minimisation Entirely

**The mistake:** Going straight from solvation to MD without minimisation.

**What happens:** The initial structure has steric clashes (atoms overlapping) that produce enormous forces (>10⁶ kJ/mol/nm). The first MD step sends atoms flying at unrealistic velocities, leading to "atom out of domain" errors or immediate crashes.

**Fix:** Always run steepest descent minimisation:
```
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000
```

### 2.2 Maximum Force Not Converging

**The mistake:** The maximum force remains extremely high (>10⁵ kJ/mol/nm) after minimisation.

**Common causes:**
| Cause | Diagnosis | Fix |
|---|---|---|
| Overlapping atoms | Visualise in PyMOL/VMD | Remove clashing atoms, rebuild region |
| Missing parameters | Check `grompp` warnings | Add missing parameters or fix topology |
| Very bad initial structure | Check PDB quality | Use a validated, cleaned structure |
| Charged residue next to vacuum | Check box edges | Increase box size |

**Target:** Maximum force < 1,000 kJ/mol/nm (strict) or < 500 kJ/mol/nm (ideal). Values > 10,000 indicate a problem. Energy minimisation does not need to reach `emtol` — check the visual structure and the force magnitude.

### 2.3 Minimisation "Succeeds" but System Is Wrong

**The mistake:** Energy minimisation converges (low forces) but the system has a fundamental error (wrong protonation, missing bonds, wrong force field).

**Why this is dangerous:** Minimisation simply relaxes the *given* system. It does not fix:
- Missing residues
- Wrong atom types
- Incorrect charges
- Mismatched topology/coordinates

**How to avoid it:** Always visually inspect the minimised structure. Load it in a molecular viewer and check:
- Does the protein look intact?
- Are water molecules filling the box uniformly?
- Are ions present and reasonable placed?
- Are there no disconnected fragments?

---

## 3. Equilibration Pitfalls

### 3.1 Skipping or Rushing Equilibration

**The mistake:** Running only a few picoseconds of equilibration, or skipping directly to production MD.

**Why equilibration matters:**
- Solvent must reorganise around the protein (solvation shells form).
- Box volume must adjust to the correct density.
- Temperature gradients must dissipate.
- Position restraints allow gradual relaxation.

**Minimum recommended equilibration:**
| Phase | Duration | What to Check |
|---|---|---|
| NVT | 100–500 ps | Temperature stable at target |
| NPT | 100–500 ps | Density stable at ~1,000 kg/m³ |
| Unrestrained NPT (optional) | 1–5 ns | RMSD levels off |

### 3.2 Using Parrinello-Rahman from the Start

**The mistake:** Using Parrinello-Rahman barostat for NPT equilibration from a non-equilibrated state.

**What happens:** Parrinello-Rahman is a "flying box" barostat — it gives the box volume its own inertia. If the density is far from equilibrium, the box oscillates wildly, potentially crashing the simulation.

**Fix:** Always use **Berendsen** for NPT equilibration (fast, monotonic convergence), then switch to **Parrinello-Rahman** for production:

```
; NPT equilibration
pcoupl = Berendsen

; Production MD
pcoupl = Parrinello-Rahman
```

### 3.3 Forgetting Position Restraints

**The mistake:** Not applying position restraints during NVT/NPT equilibration by forgetting the `-DPOSRES` flag.

**The fix:**
```
; In your .mdp file for equilibration
define = -DPOSRES
```

Without this, the `#ifdef POSRES` block in the topology is inactive and the protein is free to move during equilibration — potentially unfolding before the solvent is properly relaxed.

### 3.4 Temperature Not Stabilising

**The mistake:** Temperature oscillates or drifts during NVT equilibration.

**Common causes:**
| Cause | Fix |
|---|---|
| Bad initial velocities | Let `pdb2gmx` or `grompp` generate them (`gen_vel = yes`) |
| `tau_t` too large or too small | Use 0.1 ps for V-rescale |
| Wrong `tc-grps` | Ensure the groups exist and cover all atoms |
| System not minimised | Go back and minimise |

**Diagnostic:**
```bash
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
xmgrace temperature.xvg
```

The temperature should stabilise within the first 10–50 ps and fluctuate around the target (e.g. 300 ± 5 K for a medium-sized system).

### 3.5 Density Not Converging in NPT

**The mistake:** Box volume does not stabilise during NPT equilibration.

**Common causes:**
| Cause | Fix |
|---|---|
| `tau_p` too large | Use 2.0 ps for Berendsen, 5.0 ps for PR |
| Wrong compressibility | Use 4.5e-5 bar⁻¹ for liquid water |
| System not equilibrated in NVT first | Run proper NVT before NPT |
| Bad topology (wrong molecule counts) | Check `[ molecules ]` |

**Target density:** A solvated protein system at 300 K and 1 bar should converge to ~993–1,005 kg/m³, depending on the water model.

---

## 4. Production MD Pitfalls

### 4.1 Berendsen Thermostat/Barostat in Production

**The mistake:** Using Berendsen coupling in production simulations.

**Why it's wrong:** Berendsen does not generate the correct statistical mechanical ensemble. It suppresses fluctuations in kinetic energy (thermostat) and volume (barostat), which means:
- Heat capacity calculations will be wrong.
- Volume compressibility will be wrong.
- Any property that depends on fluctuations (many free energy methods, thermal expansion) will be incorrect.

**Fix:**
```
; Production MD settings
tcoupl = V-rescale       ; NOT Berendsen
pcoupl = Parrinello-Rahman  ; NOT Berendsen
```

> **This is the single most common mistake in published MD studies.** Reviewers increasingly check for it.

### 4.2 LINCS Warnings and Crashes

**The mistake:** Ignoring or not understanding LINCS constraint warnings.

**Typical message:**
```
Step 12345, time 24.690: LINCS WARNING
relative constraint deviation after LINCS:
rms 0.000345, max 0.0156
```

**What it means:** LINCS failed to fully satisfy a bond constraint — a bond deviated from its target length by more than the tolerance.

**Escalation:**
- **Warnings:** Occasional warnings are usually harmless. Many per step → problem.
- **Crash:** `LINCS ERROR` → bonds are so far from equilibrium that LINCS cannot correct them. Simulation terminates.

**Common causes and fixes:**
| Cause | Fix |
|---|---|
| Bad initial structure | Better minimisation (more steps, or switch to double precision) |
| Time step too large | Reduce to 1 fs temporarily |
| Trapped atom (clashing) | Visualise the frame before the crash: `gmx trjconv -dump <time>` |
| Missing parameters | Check `grompp` for NOTE/WARNING about missing dihedrals or pairs |
| Hot atoms from restraint release | Use a gradual restraint-release protocol |

### 4.3 Atom Out of Domain / Coordinate Reset Errors

**The typical error:**
```
Fatal error:
step 5678: atom 1234 in domain decomposition cell X Y Z has coordinate value > 1e10
```

**What it means:** An atom has been launched to an absurd position (essentially infinity). This is a symptom, not a cause.

**Root causes:**
1. **Steric clash** not resolved by minimisation.
2. **Missing bonded interaction** (e.g. wrong topology for a ligand).
3. **Unstable simulation** (too large a time step, wrong constraints).
4. **Sudden release of position restraints** without gradual equilibration.

**Fix:** Find the last good frame, visualise the problematic atom, and diagnose:
```bash
# Extract the frame just before the crash
gmx trjconv -s md.tpr -f md.xtc -o crash_frame.pdb -dump <time_before_crash>
```

### 4.4 Energy Drift in Long Simulations

**The mistake:** Not monitoring total energy conservation over long simulations.

**What it means:** If the conserved energy (Hamiltonian) drifts upward over time, something is wrong with the simulation parameters.

**Common causes:**
| Cause | Fix |
|---|---|
| `nstlist` too large for the buffer | Let GROMACS auto-tune (set `nstlist = 10`, Verlet scheme) |
| Cutoff too short | Use ≥ 1.0 nm for `rcoulomb` and `rvdw` |
| No dispersion correction | Add `DispCorr = EnerPres` |
| Buggy force field port | Check the force field source against the original |

**Diagnostic:**
```bash
echo "Conserved-En." | gmx energy -f md.edr -o conserved.xvg
```

Acceptable drift: < 0.01 kJ/mol/atom over the entire simulation.

### 4.5 Insufficient Output Frequency

**The mistake:** Writing trajectory frames too infrequently (or too frequently).

| Problem | Setting | Consequence |
|---|---|---|
| **Too infrequent** (`nstxout-compressed = 500000` → 1 ns between frames) | Miss fast events, poor time resolution for analysis | Can't calculate autocorrelation functions, miss transient contacts |
| **Too frequent** (`nstxout-compressed = 50` → 0.1 ps between frames) | Enormous trajectory files, slow I/O | Fills disk, slows simulation due to writing |

**Recommended starting points:**
| Output | Setting | Effective Interval | File |
|---|---|---|---|
| Compressed coords | `nstxout-compressed = 5000` | 10 ps | `.xtc` |
| Energy | `nstenergy = 5000` | 10 ps | `.edr` |
| Log | `nstlog = 5000` | 10 ps | `.log` |
| Full-precision (if needed) | `nstxout = 0`, `nstvout = 0` | Off | `.trr` |

Adjust based on your analysis needs. For fast processes (e.g. water dynamics), write more frequently for a short sub-trajectory.

### 4.6 Wall-Clock Timeouts on HPC Clusters

**The mistake:** Submitting a GROMACS job without accounting for cluster wall-time limits.

**What happens:** The job is killed mid-step. If no checkpoint was written recently, you lose hours of work. Restarting from the last checkpoint wastes the computation since the last write.

**Fix:** Use the `-maxh` flag to tell GROMACS to stop cleanly before the wall-time limit:
```bash
# If your cluster wall time is 48 hours, stop at 47.5 hours
gmx mdrun -deffnm md -maxh 47.5
```

GROMACS will write a checkpoint and stop gracefully, allowing a clean restart:
```bash
gmx mdrun -deffnm md -cpi md.cpt -maxh 47.5
```

Also reduce the checkpoint interval for long runs:
```bash
gmx mdrun -deffnm md -cpt 15   # Checkpoint every 15 minutes (default)
```

---

## 5. Analysis Pitfalls

### 5.1 Analysing Non-Equilibrated Data

**The mistake:** Including the equilibration period in your production analysis.

**Why it's wrong:** The first part of a trajectory is biased by the initial (non-equilibrium) conditions. RMSD, radius of gyration, and other properties are still drifting. Including this data skews averages and inflates error bars.

**How to avoid it:**
1. Plot RMSD vs. time. Identify when it plateaus.
2. Discard at minimum the first 10–20% of the trajectory (or more if the RMSD hasn't stabilised).
3. Use the `-b` flag in GROMACS analysis tools:

```bash
# Skip the first 10 ns
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -b 10000  # time in ps
```

### 5.2 PBC Artefacts in Trajectories

**The mistake:** Analysing or visualising a raw trajectory without fixing periodic boundary condition artefacts.

**What you see:** The protein appears to "jump" across the box, split into fragments, or have water molecules flying around disconnected from the bulk.

**This is not a bug** — it's how PBC works. Atoms that leave one side of the box re-enter on the opposite side. The raw trajectory stores coordinates within the primary box.

**Fix with `gmx trjconv`:**
```bash
# Step 1: Make the protein whole (reconnect broken chains)
echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole

# Step 2: Centre the protein in the box
echo "Protein System" | gmx trjconv -s md.tpr -f md_whole.xtc -o md_center.xtc -pbc mol -center

# Step 3: (Optional) Remove box jumps for multi-molecule systems
echo "System" | gmx trjconv -s md.tpr -f md_center.xtc -o md_nojump.xtc -pbc nojump
```

> **Critical:** Always fix PBC *before* computing RMSD, distances, contacts, or visualising the trajectory.

### 5.3 Fitting and Alignment Errors

**The mistake:** Computing RMSD without proper least-squares fitting, or fitting to the wrong group.

**Common errors:**
| Mistake | Consequence | Correct Approach |
|---|---|---|
| No fitting at all | RMSD reflects translation + rotation + internal motion | Always fit (least-squares superposition) |
| Fit to all atoms including side chains | Noisy fit dominated by flexible loops | Fit to **backbone** or **Cα atoms** |
| Fit to one chain, compute RMSD of another | Measures relative motion, not internal variation | Be explicit about fit group vs. RMSD group |

```bash
# Fit to backbone, compute RMSD of backbone
echo "Backbone Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd_bb.xvg
```

### 5.4 Drawing Conclusions from a Single Trajectory

**The mistake:** Running one simulation and reporting the results as definitive.

**Why it's wrong:** MD trajectories are chaotic — two runs with different initial velocities will diverge within picoseconds and may explore different parts of the energy landscape. A single trajectory is one sample from a vast ensemble.

**Best practice:**
- Run **at least 3 independent replicas** (different random velocity seeds).
- Report **means ± standard errors** across replicas.
- If replicas disagree qualitatively (e.g. different final conformations), you need longer simulations or enhanced sampling.

```
; Different random seeds for replicas
gen_vel     = yes
gen_temp    = 300
gen_seed    = 12345   ; Change this for each replica: 12345, 67890, 24680, ...
```

### 5.5 RMSD Misinterpretation

**Common misinterpretations:**

| Statement | Problem |
|---|---|
| "The RMSD is 0.2 nm, so the protein is stable" | RMSD measures *deviation from a reference*, not stability. A stable but different conformation also has a stable RMSD. |
| "The RMSD increased to 0.5 nm — the protein is unfolding" | Large RMSD can come from a single flexible loop moving. Compute per-residue RMSF (fluctuation) to identify the source. |
| "The RMSD plateaued at 50 ns, so the simulation is converged" | RMSD is a global, rough metric. It can plateau even if local interactions are still rearranging. Use additional metrics (RMSF, contacts, secondary structure, radius of gyration). |

### 5.6 Misusing Hydrogen Bond Cutoffs

**The mistake:** Using different geometric criteria than the original analysis or comparing H-bond counts from studies that used different cutoffs.

**Standard GROMACS defaults for `gmx hbond`:**
- Donor–Acceptor distance: ≤ 0.35 nm
- Donor–H…Acceptor angle: ≤ 30°

Always report your cutoffs. Different tools and publications use different values.

---

# Part B — Mastering GROMACS

---

## 6. Level 1: Competent User

**Goal:** Run a standard protein-in-water simulation end-to-end without errors.

**Skills to develop:**

- [ ] Run `pdb2gmx`, `editconf`, `solvate`, `genion`, `grompp`, `mdrun` in sequence.
- [ ] Understand what every file in the workflow contains (see [Key Files tutorial](Understanding_GROMACS_Key_Files.md)).
- [ ] Write basic `.mdp` files for EM, NVT, NPT, and production MD.
- [ ] Use `gmx energy` to extract and plot temperature, pressure, density, and potential energy.
- [ ] Use `gmx rms`, `gmx rmsf`, `gmx gyrate` for basic structural analysis.
- [ ] Fix PBC artefacts with `gmx trjconv`.
- [ ] Visualise trajectories in VMD, PyMOL, or ChimeraX.
- [ ] Submit and manage GROMACS jobs on an HPC cluster (e.g. Apocrita).
- [ ] Restart a simulation from a checkpoint.

**How to know you've reached this level:** You can reproduce a standard tutorial (e.g. Lemkul's lysozyme tutorial) independently without referring to the instructions.

---

## 7. Level 2: Proficient Practitioner

**Goal:** Customise simulations, handle non-standard systems, and perform rigorous analysis.

**Skills to develop:**

- [ ] Set up systems with **ligands** (parameterise with CGenFF/ACPYPE/GAFF).
- [ ] Set up **membrane protein** simulations using CHARMM-GUI.
- [ ] Use `gmx make_ndx` to create custom index groups for analysis.
- [ ] Perform **free energy perturbation** (alchemical transformations) with GROMACS free energy tools.
- [ ] Use **position restraint release** protocols (gradual relaxation over multiple equilibration stages).
- [ ] Monitor and diagnose LINCS warnings, energy drift, and pressure fluctuations.
- [ ] Run **replica simulations** and compute proper error estimates.
- [ ] Perform **Principal Component Analysis (PCA)** with `gmx covar` and `gmx anaeig`.
- [ ] Use **trjconv** fluently for all export and centering tasks.
- [ ] Write analysis scripts in Python using MDAnalysis or MDTraj.

**How to know you've reached this level:** You can set up a novel system (not from a tutorial) and defend every parameter choice.

---

## 8. Level 3: Expert / Power User

**Goal:** Push the boundaries of what GROMACS can do. Tackle complex systems and methods.

**Skills to develop:**

- [ ] Use **enhanced sampling methods**: replica exchange MD (REMD), metadynamics (via PLUMED), accelerated MD, steered MD.
- [ ] Perform **binding free energy calculations** (alchemical and endpoint methods like MM/PBSA or MM/GBSA).
- [ ] Set up **coarse-grained simulations** with MARTINI and back-map to atomistic resolution.
- [ ] Use **non-equilibrium free energy** methods (Jarzynski equality, Crooks fluctuation theorem).
- [ ] Interface GROMACS with **PLUMED** for collective variable biasing.
- [ ] Perform **constant-pH MD** (λ-dynamics or discrete protonation state methods).
- [ ] Write custom analysis tools that handle GROMACS file formats directly (XTC reader, TPR parser).
- [ ] Debug topology issues by reading `.itp` files line by line.
- [ ] Benchmark and **performance-tune** simulations for specific hardware.
- [ ] Use **GPU offloading** effectively (PME on GPU, bonded on GPU, update on GPU).

**How to know you've reached this level:** You can read a methods section in a paper, reproduce the simulation, and identify whether the choices were well-justified.

---

## 9. Level 4: Developer-Level Mastery

**Goal:** Extend GROMACS, contribute to the codebase, or develop new methods.

**Skills to develop:**

- [ ] Read and navigate the [GROMACS source code](https://gitlab.com/gromacs/gromacs) (C++).
- [ ] Understand the domain decomposition parallelisation scheme.
- [ ] Implement new analysis tools or force field terms.
- [ ] Develop custom integrators or coupling methods.
- [ ] Contribute patches upstream.
- [ ] Understand GPU kernel optimisation (CUDA/SYCL).

---

## 10. Essential Command-Line Mastery

### The Most Important GROMACS Commands

```bash
# ── SYSTEM PREPARATION ──
gmx pdb2gmx    # Convert PDB → topology + processed structure
gmx editconf   # Edit box shape, size, centres
gmx solvate    # Add solvent molecules
gmx genion     # Replace solvent with ions (neutralise + add salt)
gmx grompp     # Pre-process: combine .mdp + .gro + .top → .tpr

# ── SIMULATION ──
gmx mdrun      # Run the simulation

# ── TRAJECTORY PROCESSING ──
gmx trjconv    # Fix PBC, centre, extract frames, convert formats
gmx trjcat     # Concatenate multiple trajectories
gmx eneconv    # Concatenate energy files

# ── ANALYSIS ──
gmx energy     # Extract energy terms (T, P, density, potential, kinetic)
gmx rms        # RMSD vs. reference
gmx rmsf       # Per-residue fluctuation (RMSF)
gmx gyrate     # Radius of gyration
gmx hbond      # Hydrogen bond analysis
gmx sasa       # Solvent-accessible surface area
gmx mindist    # Minimum distance between groups
gmx covar      # Covariance analysis (PCA)
gmx anaeig     # PCA projection and cosine content
gmx cluster    # Conformational clustering
gmx do_dssp    # Secondary structure assignment (requires DSSP)
gmx densmap    # 2D density maps
gmx msd        # Mean squared displacement (diffusion)

# ── UTILITIES ──
gmx make_ndx   # Create custom index files
gmx check      # Verify trajectory/topology files
gmx dump       # Print human-readable contents of binary files (.tpr, .cpt)
gmx view       # Quick trajectory viewer (basic)
```

### Command Patterns

```bash
# Always run grompp before mdrun
gmx grompp -f params.mdp -c structure.gro -p topol.top -o run.tpr
gmx mdrun -deffnm run

# Restart from checkpoint
gmx mdrun -deffnm md -cpi md.cpt

# Extract specific time range
gmx trjconv -f md.xtc -s md.tpr -o subset.xtc -b 50000 -e 100000

# Dump TPR to check settings
gmx dump -s md.tpr | grep -i "nsteps\|dt\|ref-t\|pcoupl\|tcoupl"
```

---

## 11. Building a Robust Simulation Workflow

A disciplined workflow prevents most pitfalls:

```
1. STRUCTURE PREPARATION
   │ ├── Download from PDB
   │ ├── Check for missing residues → model with MODELLER/AlphaFold
   │ ├── Determine protonation states → H++/PropKa
   │ ├── Document all decisions in a README
   │ └── Visualise and inspect
   │
2. TOPOLOGY GENERATION
   │ ├── gmx pdb2gmx (choose FF + water model + protonation states)
   │ ├── Verify disulfide bonds in output
   │ └── Inspect .top and .itp files
   │
3. BOX & SOLVATION
   │ ├── gmx editconf (box shape, size ≥ 1.0 nm padding)
   │ ├── gmx solvate
   │ ├── gmx grompp + gmx genion (neutralise, add salt)
   │ └── Verify [ molecules ] in .top
   │
4. ENERGY MINIMISATION
   │ ├── Run steepest descent
   │ ├── Check: Fmax < 1000 kJ/mol/nm
   │ ├── Visualise minimised structure
   │ └── (Optional) conjugate gradient
   │
5. NVT EQUILIBRATION
   │ ├── Position restraints ON (-DPOSRES)
   │ ├── V-rescale thermostat, gen_vel = yes
   │ ├── 100–500 ps
   │ └── Check: temperature stable at target
   │
6. NPT EQUILIBRATION
   │ ├── Position restraints ON
   │ ├── Berendsen barostat (fast convergence)
   │ ├── 100–500 ps
   │ └── Check: density stable (~1000 kg/m³)
   │
7. PRODUCTION MD
   │ ├── V-rescale thermostat
   │ ├── Parrinello-Rahman barostat
   │ ├── No position restraints
   │ ├── continuation = yes, gen_vel = no
   │ ├── Run 3+ independent replicas (different gen_seed)
   │ └── Use -maxh for cluster time limits
   │
8. ANALYSIS
   ├── Fix PBC (trjconv -pbc whole/mol/nojump)
   ├── Discard equilibration period (-b)
   ├── Compute observables across all replicas
   ├── Report means ± standard errors
   └── Visualise final structures and trajectories
```

---

## 12. Validation Checklist: The "Pre-Flight Check"

Run through this checklist **before** submitting a production simulation:

### Structure
- [ ] All missing residues have been modelled or intentionally omitted (and documented).
- [ ] Protonation states have been determined for the target pH.
- [ ] All disulfide bonds are present and correctly assigned.
- [ ] The structure has been visually inspected in a molecular viewer.

### Topology
- [ ] Force field choice is justified and documented.
- [ ] Water model matches the force field.
- [ ] `[ molecules ]` in `.top` matches the `.gro` file.
- [ ] `grompp` produces **zero fatal errors** and you have checked all warnings.
- [ ] For ligands: parameters are from an appropriate source (CGenFF, GAFF2, manual QM fitting).

### Simulation Box
- [ ] Minimum solute-to-box-edge distance ≥ 1.0 nm (**≥ rcoulomb**).
- [ ] System is electrically neutral (total charge ≈ 0).
- [ ] Physiological salt concentration added if relevant (0.15 M NaCl).

### Energy Minimisation
- [ ] Maximum force < 1,000 kJ/mol/nm.
- [ ] Minimised structure visually inspected.

### Equilibration
- [ ] NVT: temperature stable at target (plot and verify).
- [ ] NPT: density converged to ~1,000 kg/m³ (plot and verify).
- [ ] Position restraints were active during equilibration.

### Production `.mdp` Settings
- [ ] `tcoupl = V-rescale` or `nose-hoover` (NOT Berendsen).
- [ ] `pcoupl = Parrinello-Rahman` (NOT Berendsen).
- [ ] `continuation = yes` and `gen_vel = no`.
- [ ] `constraints = h-bonds` and `dt = 0.002`.
- [ ] `coulombtype = PME`.
- [ ] `cutoff-scheme = Verlet`.
- [ ] `DispCorr = EnerPres`.
- [ ] Output frequencies set appropriately.

### HPC / Logistics
- [ ] `-maxh` set to slightly less than wall time.
- [ ] Checkpoint interval is reasonable (`-cpt 15`).
- [ ] Disk space is sufficient for trajectory output.
- [ ] Multiple replicas planned with different random seeds.

---

## 13. Performance Tuning

### Quick Wins

| Optimisation | How | Impact |
|---|---|---|
| Use GPU | `-nb gpu -pme gpu -bonded gpu` | 2–10× speedup |
| Update on GPU | `-update gpu` (GROMACS 2023+) | Further 10–30% on single GPU |
| Dodecahedron box | `-bt dodecahedron` in `editconf` | ~29% fewer atoms → proportional speedup |
| Constrain H-bonds | `constraints = h-bonds`, `dt = 0.002` | 2× vs. 1 fs time step |
| Auto-tune `nstlist` | `nstlist = 10` with Verlet | GROMACS auto-optimises the buffer |
| Tune MPI/OpenMP split | Test different `-ntmpi` / `-ntomp` ratios | Can be 20–50% difference |

### GPU Performance Checklist

```bash
# Check which GPU GROMACS sees
gmx mdrun -deffnm md -nb gpu -v 2>&1 | head -50

# Single-GPU optimal settings (typical)
gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu -ntomp 8 -pin on

# Multi-GPU (one rank per GPU)
gmx mdrun -deffnm md -nb gpu -pme gpu -ntmpi 2 -ntomp 8 -pin on -gpu_id 01
```

### Benchmarking

Always benchmark before committing to a long production run:

```bash
# Run a short benchmark (5000 steps = 10 ps)
gmx mdrun -deffnm md -nsteps 5000 -resetstep 1000 -v

# Check performance at the end of the log
tail -5 md.log
# Look for: Performance: XX.XX ns/day
```

**Performance expectations (2026 hardware):**

| System Size | Hardware | Typical Performance |
|---|---|---|
| ~50,000 atoms | Single modern GPU (e.g. A100) | 100–300 ns/day |
| ~50,000 atoms | Single consumer GPU (e.g. RTX 4090) | 150–400 ns/day |
| ~200,000 atoms | Single A100 | 30–80 ns/day |
| ~1,000,000 atoms | Multi-GPU node | 5–20 ns/day |

---

## 14. Troubleshooting Decision Tree

```
Simulation crashed
├── Error: "coordinate X Y Z not finite"
│   └── Atom launched to infinity → visualise last frame, check for clash/missing param
├── Error: "LINCS WARNING / ERROR"
│   ├── Occasional warnings → probably OK, but increase lincs_order if frequent
│   └── LINCS ERROR → reduce dt, improve minimisation, check topology
├── Error: "number of coordinates does not match topology"
│   └── [ molecules ] mismatch → check .top, recount after solvate/genion
├── Error: "atom X not found in residue Y in rtp database"
│   └── Non-standard residue → needs custom .rtp entry or separate parameterisation
├── Error: "can not do PME with mixed charges"
│   └── All charges are zero → check topology for missing charges
├── Error: "the cut-off length is X but..."
│   └── Box too small for cutoff → increase box, use editconf -d 1.2
├── Energy drifting upward
│   ├── Check nstlist/Verlet buffer
│   ├── Check DispCorr
│   └── Check for integration errors (conserved-energy from gmx energy)
├── Temperature not stable
│   ├── Check tc-grps cover all atoms
│   ├── Check tau_t is reasonable (0.1 ps for V-rescale)
│   └── Check gen_vel settings
├── Density wrong
│   ├── Check water model
│   ├── Check box geometry
│   └── Check barostat settings
└── Simulation very slow
    ├── Check GPU offloading (`-nb gpu -pme gpu`)
    ├── Check domain decomposition warnings
    ├── Check if writing output too frequently
    └── Benchmark with `-nsteps 5000 -resetstep 1000`
```

---

## 15. Resources for Continued Learning

### Official Documentation
- [GROMACS Manual](https://manual.gromacs.org/) — Definitive reference for every option and algorithm.
- [GROMACS How-To Guides](https://manual.gromacs.org/current/how-to/index.html) — Practical guides for specific tasks.

### Tutorials
- [Justin Lemkul's GROMACS Tutorials](http://www.mdtutorials.com/) — The best step-by-step tutorials for beginners. Covers lysozyme in water, membrane proteins, free energy, and umbrella sampling.
- [CHARMM-GUI](https://www.charmm-gui.org/) — Automated system builder with GROMACS output. Excellent for membrane systems, ligand parameterisation, and complex setups.
- [BioExcel Building Blocks (BioBB)](https://mmb.irbbarcelona.org/biobb/) — Python wrappers for GROMACS workflows.

### Community
- [GROMACS Users Mailing List](https://mailman-1.sys.kth.se/mailman/listinfo/gromacs.org_gmx-users) — Searchable archive of thousands of Q&A threads.
- [BioExcel Forum](https://ask.bioexcel.eu/) — Active Q&A forum for GROMACS questions.
- [GROMACS GitLab](https://gitlab.com/gromacs/gromacs) — Source code, issue tracker.

### Books
- **Frenkel & Smit**, *Understanding Molecular Simulation* — The theory behind every algorithm.
- **Braun et al.** (2019), "Best Practices for Foundations in Molecular Simulations" [*Living J. Comp. Mol. Sci.* 1(1), 5957] — Essential reading for every new MD practitioner.

### Within This Tutorial Series
- [Understanding Protein Files](Understanding_Protein_Files.md) — PDB/mmCIF formats and Python manipulation.
- [Understanding GROMACS Key Files & Force Field Models](Understanding_GROMACS_Key_Files.md) — File formats and force field comparison.
- [Understanding GROMACS Theory](Understanding_GROMACS_Theory.md) — The biophysics behind the algorithms.
- [MD on Apocrita](../MD_On_Apocrita.md) — Practical guide to running GROMACS on QMUL's HPC cluster.

---

*Reference sheet written for the QMUL Molecular Dynamics Teaching Project. Last updated: March 2026.*