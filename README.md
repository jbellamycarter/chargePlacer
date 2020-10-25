# chargePlacer: Python scripts for charge distribution in Gas-Phase molecular dynamics

This is a command line script to determine a reasonably energy minimised proton sequence for an input PDB file (**input**) for a given charge state (**charge**). See [Basic Usage](#basic-usage).
A search algorithm is used to sample proton permutations across chargeable side-chains and termini represented as point charges. The algorithm is inspired by the method described by refs [1–3](#references), it produces a reproducible output proton sequence in far fewer steps than required for sampling all permutations.

### Citation

If you use chargePlacer in your research please cite the following article:

> Bellamy-Carter J, O'Grady L, Passmore M, Jenner M and Oldham NJ. Decoding Protein Gas‐Phase Stability with Alanine Scanning and Collision Induced Unfolding Ion Mobility Mass Spectrometry. *Anal. Sens.* (2020). Accepted Author Manuscript. [doi:10.1002/anse.202000019](https://doi.org/10.1002/anse.202000019)

## Other Info

A choice of energies are calculated and used for determination:
* By default, `E_tot` is used. This is the Coulomb energy minus the proton binding energy (i.e. the summed proton affinities of protonated residues). The values used herein are derived from simplified versions of each amino acid, per ref [4](#references). These are PA<sub>NT</sub>=886.6, PA<sub>ASP-</sub>=1453.5, PA<sub>GLU-</sub>=1448.5, PA<sub>HIS</sub>=958, PA<sub>LYS</sub>=918, PA<sub>ARG</sub>=1002, PA<sub>CT-</sub>=1430 kJ/mol.
* `Coulomb-only` is the alternative mode (activated by `-c`), where only the Coulomb energy is taken into account.

This software also provides the option to perform in silico alanine scanning (activated by `-a`), where each chargeable side-chain is removed and the minimised proton sequence determined for each 'mutant'.

Both default and `alanine_scan` modes output two tab-separated files `proton_sites.txt` and `charges.txt`.
* `proton_sites.txt` : A text file listing the side-chains and termini that are protonated in the energy minimised sequence. Each row has the form `<RESN> <CHAIN> <RESI>`. The energies calculated are included in the file header. This file only ever contains the data for the 'wild-type' protein sequence.
* `charges.txt` : A text file containing the charges of each residue and terminus in the energy minimised sequence. For the `alanine_scan` method, each subsequent row contains each mutant variant with the mutated residue charge indicated by `nan`.
An additional file `energies.txt` is generated in `alanine_scan` mode, this contains the calculated energies for each charge sequence in `charges.txt`.

## Basic Usage

In the simplest case, a PDB file **`input`** is provided along with a target charge state **`charge`**.
`chargePlacer` will then import the atom coordinates for pre-defined point charge atoms for all chargeable residues and termini. 

```shell
python chargePlacer.py input charge
```

If you move `chargePlacer` to a directory in your system `$PATH` then you may omit the `python` command:

```shell
chargePlacer.py input charge
```

To minimise for Coulomb energy only, append `-c`:

```shell
python chargePlacer.py input charge -c
```

To perform an `alanine_scan`, append `-a`:

```shell
python chargePlacer.py input charge -a
```

If you would like to capture the full output of the program you can use piping:

```shell
python chargePlacer.py input charge > log.txt
```

this may be especially powerful when combined with the `verbose` option (`-v`), which will print the energies calculated at each step.

```shell
python chargePlacer.py input charge -v > log.txt
```

See [Command Line Options](#command-line-options) below for further details. 

### Command Line Options

#### Required Arguments
| Parameter | Description                                   |
|-----------|-----------------------------------------------|
| **input**, e.g. `input.pdb`     | input PDB file for which to determine charges |
| **charge**, e.g. `7`    | target charge state                           |

#### Optional arguments
| Option                    |  Description                    |
|---------------------------|---------------------------------|
| `-h`, `--help`            | show this help message and exit |
| `-v`, `--verbose`         | verbose output |
| `-c`, `--coulomb_only`    | minimise for Coulomb repulsion only, ignores proton affinity |
| `-r` , `--relative_permittivity` | relative permittivity to use (default: 1) |
| `-a`, `--alanine_scan`    | perform in silico alanine scanning for all chargeable residues |
| `-p`, `--protect`         | protect listed chains from alanine scanning, e.g. `-p ABC` |
| `-o`, `--output` OUTPUT | prefix for output files (default: ""). Gives `*proton_sites.txt` and `*charges.txt` |


## References
1. V. Popa, D. A. Trecroce, R. G. McAllister and L. Konermann, J. Phys. Chem. B, 2016, 120, 5114–5124.
2. M. Bakhtiari and L. Konermann, J. Phys. Chem. B, 2019, 123, 1784–1796.
3. S. K. Fegan and M. Thachuk, J. Chem. Theory Comput., 2013, 9, 2531–2539.
4. A. Moser, K. Range and D. M. York, J. Phys. Chem. B, 2010, 114, 13911–13921.
