# Gas-phase MD scripts to accompany chargePlacer

A set of scripts for performing consistent gas-phase MD simulations of protein structures with GROMACS, particularly those with proton assignments predicted by **`chargePlacer`**.

## Usage

Perform **`chargePlacer`** calculations on an initial PDB file (`protein_init.pdb`) and generate the corresponding `proton_sites.txt` file for a given charge state. Copy the contents of this folder into the same directory as `protein_init.pdb` and `proton_sites.txt`. Then run:
```bash ./gmx_command_list.sh
```

This will iterate through the following steps automatically:
1. Apply the charges determined by **`chargePlacer`** to `protein_init.pdb` (creating `protein.gro`).
2. Expand the boundary box around the protein to a 900x900x900nm cube (creating `protein_box.gro`).
3. Minimise the structure, according to `min.mdp` (creating `min.gro`, and plotting the energy curve `min_energy.png`).
4. Create 4 new directores for the room temperature (RT) and thermal gradient (A, B and C) simulations.
5. Perform an equilibration at 298 K for each of these directories, according to `eq.mdp` (creating `eq.gro`)
6. Run production simulations at 298 K (RT) and 298–798 K (A, B and C), according to `prod_298K.mdp` and `prod_798K.mdp`, respectively.
7. Run `gmx` analysis functions on the production simulations and plot them.

## Notes
* The automated charge assignment currently does not process chain identity, therefore multi-chain structures should be renumbered sequentially such that each residue has a unique number.
* The scripts are written in a combination of `bash`, `expect` and `gnuplot`, these must all be installed for the script to function.
