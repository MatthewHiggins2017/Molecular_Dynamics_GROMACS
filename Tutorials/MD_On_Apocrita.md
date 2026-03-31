# Molecular Dynamics Simulation of a Small Protein Using GROMACS on Apocrita

This tutorial walks you through running a molecular dynamics (MD) simulation of Factor Xa (PDB: 1FJS) entirely via the command line on the QMUL HPC cluster Apocrita. It is adapted from the official GROMACS introductory tutorial (tutorials.gromacs.org).

**Goal:** Learn step-by-step how to set up and run an MD simulation using GROMACS on a HPC cluster.
**Software:** GROMACS 2025 (2024/2023 also work).
**Visualisation:** All plots and 3D views are handled by the companion Jupyter notebook (`MD_On_Apocrita_Complementary_Notebook.ipynb`) — launch it via [Apocrita OnDemand Jupyter](https://docs.hpc.qmul.ac.uk/ondemand/jupyter/) alongside this tutorial.

## Key References

* Apocrita GROMACS documentation: https://docs.hpc.qmul.ac.uk/apps/chem/gromacs/
* SLURM job scheduling on Apocrita: https://slurm-docs.hpc.qmul.ac.uk/
* GROMACS manual: https://manual.gromacs.org/current/index.html

---

## 1. Getting Started on Apocrita

### 1.1 Start an Interactive Session

SSH into Apocrita and request an interactive session:

```bash
salloc -n 4 --mem-per-cpu=1G -t 4:0:0
```

> **Note:** We request 4 hours here since the full tutorial takes some time. For the longer production simulation you may wish to submit a batch job instead (see Section 9).

### 1.2 Load Modules and Set Up the Environment

```bash
# Load necessary modules
module load nano miniforge gromacs
```

### 1.3 Download and Extract the Tutorial Files

```bash
# Download the tutorial
git clone https://github.com/MatthewHiggins2017/Molecular_Dynamics_GROMACS.git

# Extract and enter directory
unzip md-intro-tutorial-main.zip
cd md-intro-tutorial-main
```

### 1.4 Create and Activate the Conda Environment

Here we are using the version without gmx as have GROMACS installed as a module anyway. 

```bash
conda env create --name GROMACS-md-intro-tutorial --file environment_withoutgmx.yml
conda activate GROMACS-md-intro-tutorial
```

### 1.5 Move into the Data Directory

All work will be done inside the `data` directory, which contains the `input/` folder with prepared parameter files and PDB structures, and a `reference/` folder with pre-computed results.

```bash
cd data
ls
```

---

## 2. Obtaining and Cleaning the Input Structure

### 2.1 About the Protein

We will simulate Factor Xa, a protein critical in blood clot formation and a target for anticoagulant medicines. The 3D crystal structure is available from the RCSB (PDB code: 1FJS) and can be found at `input/1fjs.pdb`.

### 2.2 Clean the PDB File

We only want to simulate the protein itself, so we strip out all non-protein atoms (crystal waters, ligands, glycerol molecules, etc.). These are labelled `HETATM` and `CONECT` in the PDB file:

```bash
grep -v HETATM input/1fjs.pdb > 1fjs_protein_tmp.pdb
grep -v CONECT 1fjs_protein_tmp.pdb > 1fjs_protein.pdb
```

> **Note:** This procedure is not universally appropriate. Removing a tightly bound ligand or functional active-site water molecule can have significant effects on the protein conformation.

### 2.3 Check for Missing Residues/Atoms

Before proceeding with topology generation, it is **critical** to verify the structural completeness of your PDB file. Missing residues or atoms can cause downstream tools to fail or produce invalid results. 

Always check the PDB file for `MISSING` entries, which indicate atoms or residues absent from the crystal structure:

```bash
grep MISSING input/1fjs.pdb
```

**What to look for:**
- **Terminal regions** (N- or C-terminus): Often absent in crystal structures and usually tolerable for MD simulations.
- **Internal loops or domains**: Missing segments within the protein can significantly affect folding dynamics and should be modelled using homology modelling or other structure prediction tools before running MD.
- **Individual atoms**: Missing atoms within an otherwise present residue will cause `gmx pdb2gmx` to fail. These must be added using preparation software (e.g., MODELLER, Swiss-Model).

This is an essential quality control step — incomplete structures lead to failed simulations or unphysical results.


-------

## 3. Generating a Topology

In simple terms, a **topology** is a complete instruction manual for how your protein is built and how all its atoms interact with each other.

Think of it like a LEGO instruction booklet:
- **Atoms** are the individual LEGO bricks
- **Bonds** are the connections between bricks
- **Force field parameters** are the rules (e.g., "this connection should be stiff", "these atoms repel each other")

The topology file tells GROMACS:
1. What atoms exist and their types
2. Which atoms are bonded together (and how strongly)
3. What angles and twists exist between bonded atoms
4. How atoms that aren't bonded should interact (electrostatic attraction/repulsion, van der Waals forces)

Without a topology, the simulation software wouldn't know how to move the atoms or calculate energies. GROMACS uses this information to simulate the protein's natural motion under the laws of physics.

The topology is generated automatically from your cleaned PDB file using `gmx pdb2gmx` — it reads the protein structure and looks up all interaction parameters in a database called the **force field** (in our case, CHARMM27).

### 3.1 Run `gmx pdb2gmx`

The first GROMACS tool we use is `gmx pdb2gmx`, which generates:
- A topology for the molecule (`.itp`)
- A position restraint file (`.itp`)
- A post-processed structure file (`.gro`)

```bash
gmx pdb2gmx -f 1fjs_protein.pdb -o 1fjs_processed.gro -water tip3p -ff "charmm27"
```

Here we chose the **CHARMM27 all-atom force field** and **TIP3P water model**. This is an important choice — always read about each force field and choose the most applicable one.

Commonly used `pdb2gmx` options:

| Option   | Effect |
|----------|--------|
| `-water` | Water model: none, spc, spce, tip3p, tip4p, tip5p, tips3p |
| `-ignh`  | Ignore H atoms in the PDB file (useful for NMR structures) |
| `-ter`   | Interactively assign charge states for N/C-termini |
| `-inter` | Interactively assign charge states for Glu, Asp, Lys, Arg, His; choose disulfide bonds |


### 3.2 Check the Generated Files

```bash
ls
```

You should see new files: `1fjs_processed.gro`, `topol.top`, `topol_Protein_chain_X.itp`, and `posre_Protein_chain_X.itp`.

> **Note:** GROMACS can handle many file formats — `.gro` is simply the default. You can output `.pdb` format by specifying a `.pdb` extension for the output file.

A breakdown of these files can be found [here](./Supporting_Information/Understanding_GROMACS_Key_Files.md)


---

## 4. Understanding the Topology

### 4.1 Inspect the Topology File

```bash
cat topol.top
```

The topology begins with comments, then includes:
1. Force field parameters
2. Protein molecule topologies (`.itp` files)
3. Water topology
4. Ion parameters
5. System-level definitions (`[ system ]` and `[ molecules ]`)

### 4.2 Inspect the Molecule Types

```bash
grep "moleculetype" -A 3 topol_Protein_chain_A.itp
grep "moleculetype" -A 3 topol_Protein_chain_L.itp
```

### 4.3 Inspect Atom Definitions

```bash
grep "atoms" -A 7 topol_Protein_chain_A.itp
```

### 4.4 Inspect Bonded Interactions

```bash
# Bonds
grep "bonds" -A 7 topol_Protein_chain_A.itp

# Angles
grep "angles" -A 7 topol_Protein_chain_A.itp

# Dihedrals
grep "dihedrals" -A 7 topol_Protein_chain_A.itp

# Pairs (1-4 interactions)
grep "pairs" -A 7 topol_Protein_chain_A.itp
```

### 4.5 Position Restraints and Water Topology

```bash
grep "posre" -A 3 -B 3 topol_Protein_chain_A.itp
grep "position_restraints" -A 7 posre_Protein_chain_A.itp
grep "water topology" -A 8 topol.top
```

### 4.6 System-Level Definitions

```bash
grep "system" -A 7 topol.top
```

The `[ molecules ]` directive lists the number and type of each molecule. The order and count must exactly match the coordinate file.

---

## 5. Solvating the System

Solvating a system in GROMACS means adding water molecules (or other solvent) around your protein or molecule of interest. This is a crucial step after generating the topology file.

**Think of it like this:** Imagine studying how a protein behaves in a test tube. In real life, that protein isn't floating in empty space—it's surrounded by water. If you tried to simulate the protein without water, you'd be missing a huge part of the story.

**Key Reasons to Solvate:**

1. **Mimics Reality**
  - Proteins and molecules exist in aqueous (water) environments in cells
  - Without solvent, your simulation won't reflect what actually happens in nature

2. **Stabilizes the System**
  - Water molecules help keep your protein in its proper 3D shape
  - They create a more realistic chemical environment

3. **Enables Important Interactions**
  - Water participates in chemical reactions and bonding
  - It shields charged parts of molecules from each other
  - These interactions affect how your molecule behaves

4. **Prevents Boundary Artifacts**
  - Water fills the space around your molecule
  - Without it, your protein might interact artificially with the simulation box edges
  - This creates unrealistic computational "mirror" effects

Solvation transforms your simulation from studying an isolated molecule in a vacuum to studying it in its natural environment—making your results scientifically meaningful and biologically accurate.


### 5.1 Define the Simulation Box

We use a rhombic dodecahedron (volume ~71% of an equivalent cubic box, saving water molecules):

![rhombic dodecahedron image](./imgs/rhombic_dodecahedron.png)

```bash
gmx editconf -f 1fjs_processed.gro -o 1fjs_newbox.gro -c -d 1.0 -bt dodecahedron
```

- `-c` centres the protein in the box
- `-d 1.0` places it at least 1.0 nm from the box edge
- `-bt dodecahedron` sets the box type

> The distance to the box edge is important: the protein must never interact with its periodic image. With a 1.2 nm non-bonded cutoff and 1.0 nm box margin, we get at least 2.0 nm between periodic images.

### 5.2 Fill the Box with Water

```bash
gmx solvate -cp 1fjs_newbox.gro -cs spc216.gro -o 1fjs_solv.gro -p topol.top
```

Check that the topology was updated:

```bash
tail topol.top
```

A `SOL` entry should now appear under `[ molecules ]` with the number of added water molecules.

> **Note:** `gmx solvate` *adds* water entries to the topology — it does not overwrite. Running it multiple times will cause a mismatch.

---

## 6. Adding Ions

### 6.1 Determine the System Charge - [COME BACK AND CHECK THIS]

Check the total charge of each chain by looking at the last `qtot` value before the `[ bonds ]` section:

```bash
grep "bonds" -B 20 topol_Protein_chain_A.itp | tail -25
grep "bonds" -B 20 topol_Protein_chain_L.itp | tail -25
```

The protein has a net charge of -2e (qtot_A = 1, qtot_L = -3).

### 6.2 Prepare the Input for `gmx genion`

Create an empty `.mdp` file and run `gmx grompp` to generate a `.tpr`:

```bash
touch ions.mdp
gmx grompp -f ions.mdp -c 1fjs_solv.gro -p topol.top -o ions.tpr
```

### 6.3 Add Ions

```bash
printf "SOL\n" | gmx genion -s ions.tpr -o 1fjs_solv_ions.gro -conc 0.15 \
  -p topol.top -pname NA -nname CL -neutral
```

We select the `SOL` group for ion placement (we never want to replace protein atoms with ions). The `-neutral` flag ensures the system is charge-neutral, and `-conc 0.15` sets the NaCl concentration to 0.15 M.

> **Warning:** Only run `gmx genion` once. It edits the topology in-place and does not check for existing ions.

Verify the topology:

```bash
tail -6 topol.top
```

You should see two more NA than CL ions (to compensate the -2e charge).

---

## 7. Energy Minimisation

### 7.1 Inspect the EM Parameter File

```bash
cat input/emin-charmm.mdp
```

### 7.2 Prepare and Run Energy Minimisation

```bash
gmx grompp -f input/emin-charmm.mdp -c 1fjs_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 4
```

> **Note on `-v` flag:** The `-v` (verbose) flag prints progress at every step. This is fine for interactive sessions but should **not** be used in batch job scripts — it produces excessive output.

### 7.3 Check the Results

Two key indicators of successful energy minimisation:
1. **Potential energy (Epot):** Should be negative, on the order of 10^5 kJ/mol
2. **Maximum force (Fmax):** Should be below the target set in `emin-charmm.mdp` (typically `emtol = 1000.0` kJ/(mol nm))

### 7.4 Analyse the Potential Energy

```bash
printf "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg
```

**Visualise:** Open the companion notebook and run the **Section 7 — Potential Energy** cell.

The plot should show steady convergence of the potential energy.

---

## 8. Equilibration

### 8.1 NVT Equilibration (Temperature)

The first equilibration phase uses the NVT ensemble (constant Number of particles, Volume, Temperature) for 100 ps with position restraints on the protein.

Inspect the parameter file:

```bash
cat input/nvt-charmm.mdp
```

Key parameters:
- `gen_vel = yes` — generates initial velocities
- `tcoupl = V-rescale` — velocity rescaling thermostat
- `pcoupl = no` — no pressure coupling
- `define = -DPOSRES` — enables position restraints

Run the NVT equilibration:

```bash
gmx grompp -f input/nvt-charmm.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -ntmpi 1 -ntomp 4 -v -deffnm nvt
```

Analyse the temperature:

```bash
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -b 20
```

**Visualise:** Open the companion notebook and run the **Section 8.1 — Temperature** cell.

The temperature should quickly reach 300 K and remain stable.


### 8.2 NPT Equilibration (Pressure)

The second equilibration phase uses the NPT ensemble (constant Number of particles, Pressure, Temperature) for 100 ps, also with position restraints.

Inspect the parameter file:

```bash
cat input/npt-charmm.mdp
```

Key changes from NVT:
- `continuation = yes` — continuing from NVT
- `gen_vel = no` — velocities are read from the checkpoint

Run the NPT equilibration:

```bash
gmx grompp -f input/npt-charmm.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -ntmpi 1 -ntomp 4 -v -deffnm npt
```

Analyse pressure and density:

```bash
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
echo "Density" | gmx energy -f npt.edr -o density.xvg
```

**Visualise:** Open the companion notebook and run the **Section 8.2 — Pressure** and **Density** cells.

Pressure will fluctuate widely (RMSE ~100 bar is normal) — what matters is that the *average* is close to the target (1 bar). The density should be close to ~1000 kg/m^3.



---

## 9. Production MD Run

The system is now equilibrated. We release position restraints and run a 1 ns production simulation.

### 9.1 Inspect the Production Parameter File

```bash
cat input/md-charmm.mdp
```

### 9.2 Prepare the Run

```bash
gmx grompp -f input/md-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
```

### 9.3 Option A: Run Interactively

For a quick test (if your interactive allocation has enough time remaining):

```bash
gmx mdrun -ntmpi 1 -ntomp 4 -v -deffnm md
```

### 9.4 Option B: Submit as a SLURM Batch Job (Recommended)

For the production run, it is better to submit a batch job. Create a script called `run_md.sh`:

```bash
nano run_md.sh
```

Paste the following content:

```bash
#!/bin/bash
#SBATCH --job-name=gromacs_md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=md_%j.out
#SBATCH --error=md_%j.err

module load miniforge gromacs
conda activate GROMACS-md-intro-tutorial

gmx mdrun -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} -deffnm md
```

Submit the job:

```bash
sbatch run_md.sh
```

Monitor the job:

```bash
squeue -u $USER
```

---

## 10. Analysis

### 10.1 Post-Process the Trajectory

The protein may diffuse across periodic boundaries during the simulation. We fix this in three steps:

```bash
# Make all molecules whole
printf "1" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole

# Remove jumps across periodic boundaries
printf "1" | gmx trjconv -s md.tpr -f md_whole.xtc -o md_nojump.xtc -pbc nojump

# Centre the protein in the box
printf "1\n1" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol
```

In each step, we select group `1` (Protein).

### 10.2 Visualise the Trajectory

**Visualise:** Open the companion notebook and run the **Section 10.2 — Trajectory** cell. This uses NGLView + MDAnalysis to give you an interactive 3D viewer with playback controls — no VMD required.

### 10.3 Check the Minimum Image Convention

Verify the protein never interacts with its own periodic image:

```bash
printf "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg
```

**Visualise:** Open the companion notebook and run the **Section 10.3 — Minimum Image Distance** cell.

The minimum distance should always exceed the non-bonded cutoff (1.2 nm).

### 10.4 RMSD (Structural Stability)

Calculate the backbone RMSD relative to the crystal structure:

```bash
printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns
```

**Visualise:** Open the companion notebook and run the **Section 10.4 — RMSD** cell.

Here, group `4` (Backbone) is used for both the least-squares fit and RMSD calculation. The RMSD should level off around 0.15 nm (1.5 A), indicating a stable structure.

> **Note:** The reference structure must have the same number of atoms as the trajectory.

### 10.5 Radius of Gyration (Compactness)

```bash
echo "1" | gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg
```

**Visualise:** Open the companion notebook and run the **Section 10.5 — Radius of Gyration** cell.

A stable Rg value indicates the protein remains in its compact, folded form.


### 10.7 Report Methods

Generate a summary of the simulation settings (useful for publications):

```bash
gmx report-methods -s md.tpr
```

---

## 11. Improving Performance with Longer Timesteps (Optional)

Two methods can increase simulation speed at the cost of slightly altered dynamics.

### 11.1 Hydrogen Mass Repartitioning

Increases the mass of hydrogen atoms (transferring mass from bonded heavy atoms) to slow down the fastest motions, allowing a 4 fs timestep:

```bash
cat input/md-mass-charmm.mdp
gmx grompp -f input/md-mass-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md-mass.tpr
gmx mdrun -v -deffnm md-mass -ntmpi 1 -ntomp 4
```

### 11.2 Multiple Time-Stepping (MTS)

Computes some force types less frequently to reduce computational cost:

```bash
cat input/md-mts-charmm.mdp
gmx grompp -f input/md-mts-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md-mts.tpr
gmx mdrun -v -deffnm md-mts -ntmpi 1 -ntomp 4
```

Typical performance improvements (on 8 CPU cores):

| Method | Performance (ns/day) |
|--------|--------------------:|
| Regular MD | 26.9 |
| Mass repartitioning (4 fs) | 52.3 |
| MTS (`longrange-nonbonded`) | 28.3 |
| MTS (`nonbonded longrange-nonbonded`) | 50.8 |

---

## 12. Summary of the Full Workflow

For quick reference, here are all the key commands in order:

```bash
# --- Setup ---
salloc -n 4 --mem-per-cpu=1G -t 4:0:0
module load nano miniforge gromacs
mkdir -p Learning_GROMACS && cd Learning_GROMACS
wget https://gitlab.com/gromacs/online-tutorials/md-intro-tutorial/-/archive/main/md-intro-tutorial-main.zip
unzip md-intro-tutorial-main.zip && cd md-intro-tutorial-main
conda env create --name GROMACS-md-intro-tutorial --file environment.yml
conda activate GROMACS-md-intro-tutorial
cd data

# --- Clean PDB ---
grep -v HETATM input/1fjs.pdb > 1fjs_protein_tmp.pdb
grep -v CONECT 1fjs_protein_tmp.pdb > 1fjs_protein.pdb

# --- Generate Topology ---
gmx pdb2gmx -f 1fjs_protein.pdb -o 1fjs_processed.gro -water tip3p -ff "charmm27"

# --- Solvate ---
gmx editconf -f 1fjs_processed.gro -o 1fjs_newbox.gro -c -d 1.0 -bt dodecahedron
gmx solvate -cp 1fjs_newbox.gro -cs spc216.gro -o 1fjs_solv.gro -p topol.top

# --- Add Ions ---
touch ions.mdp
gmx grompp -f ions.mdp -c 1fjs_solv.gro -p topol.top -o ions.tpr
printf "SOL\n" | gmx genion -s ions.tpr -o 1fjs_solv_ions.gro -conc 0.15 \
  -p topol.top -pname NA -nname CL -neutral

# --- Energy Minimisation ---
gmx grompp -f input/emin-charmm.mdp -c 1fjs_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 4

# --- NVT Equilibration ---
gmx grompp -f input/nvt-charmm.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -ntmpi 1 -ntomp 4 -v -deffnm nvt

# --- NPT Equilibration ---
gmx grompp -f input/npt-charmm.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -ntmpi 1 -ntomp 4 -v -deffnm npt

# --- Production MD ---
gmx grompp -f input/md-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -ntmpi 1 -ntomp 4 -v -deffnm md

# --- Analysis ---
printf "1" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole
printf "1" | gmx trjconv -s md.tpr -f md_whole.xtc -o md_nojump.xtc -pbc nojump
printf "1\n1" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol
printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns
echo "1" | gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg
```

---

