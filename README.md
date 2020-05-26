# chargePlacer: Python scripts for charge distribution in gas-phase molecular dynamics

This is a command line script to determine a reasonablly energy
minimised proton sequence for an input PDB file (**INPUT**) for a given
charge state (**CHARGE**). A search algorithm is used to sample proton
permutations across chargeable side-chains and termini represented as
point charges. This algorithm produces a reproducible output proton
sequence in far fewer steps than required for sampling all permutations.

The optimised proton sequence is saved to file (**OUTPUT**) with the format:
`<Residue Name> <Residue Number> <Chain Identifier>`

A choice of energies are calculated and used for determination. By
default, `E_tot` is used. This is the Coulomb energy minus the proton
binding energy (i.e. the summed proton affinities of protonated
residues). `Coulomb-only` is the alternative mode, where only the
Coulomb energy is taken into account.

This software also provides the option to perform in silico alanine
scanning, where each chargeable side-chain is removed and the minimised
proton sequence determined for each 'mutant'.
