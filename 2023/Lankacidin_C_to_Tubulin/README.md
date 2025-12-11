# Protein–Ligand MD/MMPBSA Pipeline

This folder contains the scripts and parameters used to prepare a receptor–GDP–ligand complex, run GROMACS simulations on a Slurm cluster, and compute binding free energies with g_mmpbsa.

## Layout
- `bash_scripts/` – all workflow scripts (prep, solvation, minimization, equilibration, production, continuation, concatenation, analysis, MMPBSA, Slurm job templates).
- `gromacs_mdp/` – GROMACS parameter files (`em.mdp`, `ions.mdp`, `nvt.mdp`, `npt.mdp`, `md.mdp`, `pbsa.mdp`).
- `ligand_paramatrization/` – LEaP inputs and force-field pieces for GDP/ligand (`gdpleap.in`, `LIGleap.in`, `gdp.prep`, `gdp_GMX.itp`, `frcmod.phos`, `leaprc.ff99SBildn`, `acpype.py`).
- `binding_free_energy/` – `g_mmpbsa.tar` source (build to obtain `bin/g_mmpbsa`).

## Requirements
- Slurm (`sbatch`).
- GROMACS (scripts reference 2019-fosscuda-2018b; MMPBSA scripts reference 2016.3).
- AmberTools (Antechamber/Parmchk2/TLeap; scripts used AmberTools/17).
- OpenBabel (used for ligand preprocessing).
- ACPYPE (Python 2.7).
- APBS or compatible PB solver as configured in `gromacs_mdp/pbsa.mdp`.

Edit the hard-coded paths in the scripts (e.g., `neededFiles`, `ACPYPE`, `pbsa`) to match this folder. A quick option:
```bash
export BASE_DIR=/Users/mohamed/Documents/Research/Projects/GitHub_repo/Papers/2023/neededFiles
export PARAM_DIR=$BASE_DIR/ligand_paramatrization
export MDP_DIR=$BASE_DIR/gromacs_mdp
export BFE_DIR=$BASE_DIR/binding_free_energy
```

## Workflow (run from `2023/neededFiles`)
1) **Ligand & GDP parameters** – `bash_scripts/1_ligandPreparation.sh complex.pdb`  
   Splits the input PDB into receptor (`rec.pdb`), GDP (`gdp.pdb`), and ligand (`LIG.pdb`), then runs Antechamber/ACPYPE using files in `ligand_paramatrization/`.

2) **Build topology** – `bash_scripts/2_prepareComplex.sh`  
   Runs `pdb2gmx`, merges receptor + ligand + GDP coordinates into `complex.gro`, and inserts ligand/GDP includes into `topol.top` and `LIG_GMX.itp`.

3) **Solvate & ionize** – `bash_scripts/3_solvIons.sh`  
   Uses `gromacs_mdp/ions.mdp` to solvate (`solv.gro`) and add 0.15 M ions (`solv_ions.gro`, updated `topol.top`).

4) **Energy minimization** – `bash_scripts/4_em.sh`  
   Copies `gromacs_mdp/em.mdp`, builds `em.tpr`, and submits `em.sh` via Slurm (template `bash_scripts/gpujob.sh`).

5) **Equilibration (NVT → NPT)** – `bash_scripts/5_equil.sh`  
   Generates ligand/GDP position restraints, builds `index.ndx`, then runs `nvt.mdp` followed by `npt.mdp` (both from `gromacs_mdp/`) through Slurm.

6) **Production MD**  
   - First segment: `bash_scripts/6_md1.sh 50` (for 50 ns; adjusts `md.mdp` from `gromacs_mdp/`).
   - Continuation: `bash_scripts/7_md2.sh md_00_50.xtc 50` (extends another 50 ns from previous outputs).
   - Concatenate segments: `bash_scripts/8_cat.sh` → `md_00_<end>.xtc`.

7) **Binding free energy (MMPBSA)**  
   - Quick path: `bash_scripts/9_mmpbsa.sh md_00_50.xtc` then `cd mmpbsa && sbatch mmpbsajob.sh` (uses `gromacs_mdp/pbsa.mdp`).  
   - Legacy split run: `bash_scripts/1-mmpbsa.sh` then `bash_scripts/2-mmpbsa.sh`.  
   Update `pbsa` path inside the scripts to point to the built `binding_free_energy/g_mmpbsa*/bin/g_mmpbsa`.

8) **Trajectory analysis**  
   - RMSD/Gyration/RMSF: `bash_scripts/analysis.sh md_00_50_center.xtc md_00_50.gro em.tpr` (outputs to `analysis/`).  
   - Variant with custom index: `bash_scripts/analysis-RMSD-script.sh ...`.  
   - Frame extraction: `bash_scripts/PDB-extraction.sh traj.xtc <dt_ps>` (uses `index2.ndx`).  
   - Energy breakdown: `bash_scripts/frameVsComponents.sh` post-processes g_mmpbsa outputs.

## Notes
- Slurm job templates (`gpujob.sh`, `job.sh`, headers in MMPBSA scripts) contain email/account fields—update them for your cluster.
- The scripts still default to historic module names; change the `module load` lines to your available GROMACS/AmberTools/OpenBabel/ACPYPE modules or paths.
- Run all scripts from this folder so relative paths resolve correctly after the re-organization.

## Provenance
These workflows were used for the tubulin–lankacidin studies:
- Ayoub, A. T., Elrefaiy, M. A., et al. (2022) “Bioinspired computational design of lankacidin derivatives for improvement in antitumor activity.” Future Medicinal Chemistry 14, 1349–1360. DOI:10.4155/fmc-2022-0134.
- Ayoub, A. T., Elrefaiy, M. A., & Arakawa, K. (2019) “Computational prediction of the mode of binding of antitumor lankacidin C to tubulin.” ACS Omega 4(2), 4461–4471. DOI:10.1021/acsomega.8b03470.
