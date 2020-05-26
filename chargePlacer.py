#!/usr/bin/env python3
#-*- coding:utf-8 -*-
# Copyright 2020 Jedd Bellamy-Carter
# MIT License
"""chargePlacer: Python implementation of a charge positioning algorithm

This is a command line script to determine a reasonablly energy
minimised proton sequence for an input PDB file (INPUT) for a given
charge state (CHARGE). A search algorithm is used to sample proton
permutations across chargeable side-chains and termini represented as
point charges. This algorithm produces a reproducible output proton
sequence in far fewer steps than required for sampling all permutations.

The optimised proton sequence is saved to file (*proton_sites.txt) with 
the format:
<Residue Name> <Residue Number> <Chain Identifier>

The charges matching to this proton sequences are saved to file 
(*charges.txt).

A choice of energies are calculated and used for determination. By
default, `E_tot` is used. This is the Coulomb energy minus the proton
binding energy (i.e. the summed proton affinities of protonated
residues). `Coulomb-only` is the alternative mode, where only the
Coulomb energy is taken into account.

This software also provides the option to perform in silico alanine
scanning, where each chargeable side-chain is removed and the minimised
proton sequence determined for each 'mutant'. The charges for each 
mutant proton sequence are appended to `*charges.txt`.
"""

from __future__ import division, print_function, absolute_import

# Standard Python Modules
import argparse
import os
import sys
import time

# Additional Modules
import numpy as np

#%% GLOBAL VARIABLES

_E_CONST = 2.307E-18  # J.Å -- (elementary charge)^2 / (4 x PI x vacuum permittivity)
relative_permittivity = 1 # Relative permittivity of vacuum
AVOGADRO = 6.022E+23
E_CONST = _E_CONST * 0.001 * AVOGADRO / relative_permittivity # kJ.Å/mol

#%% LOOKUP DICTIONARIES 

# Atoms to get coordinates for point-charge assignment
point_charge_dict = {
    'ASP': 'OD2',
    'GLU': 'OE2',
    'LYS': 'NZ', 
    'ARG': 'NH2', 
    'HIS': 'CB', 
    'NT': 'N', 
    'CT': 'C'
}
# Charge of the deprotonated residue
deprot_charge_dict = {
    'ASP': -1, 
    'GLU': -1, 
    'LYS': 0, 
    'ARG': 0, 
    'HIS': 0, 
    'NT': 0, 
    'CT': -1
}
# Proton affinities in kJ/mol
proton_affinity_dict = {
    'ASP': 1453.5, 
    'GLU': 1448.5, 
    'LYS': 918., 
    'ARG': 1002.,
    'HIS': 958., 
    'NT': 886.6, 
    'CT': 1430.
}

#%% CLASSES

class PDB():
    """
    Class for parsing PDB files that follow the accepted format
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

    Extracts ATOM and HETATM records. 
    HETATM with HOH resname are removed.

    Dictionary format:

    Structure (list)>
        Model (dict)>
            Chain (dict)>
                Residue (list)>
                    Atom (dict)
    """

    def __init__(self, filename=None):
        if not os.path.splitext(filename)[1] in ['.pdb', '.pdbqt']:
            raise ValueError('Incorrect file extension, must be .pdb or .pdbqt')
        else:
            self.filename = filename

        self.AA_3to1    = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
                           'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G',
                           'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
                           'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
                           'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        self.ATOM_STRING = "{}{:5d} {:4}{:.1}{:.3} {:.1}{:>4d}{:.1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:.2}  \n"

        self.structure = [None]
        self.num_models = 0
        self.num_chains = 0
        self.chains = {}
        self._parse()

    def _parse(self):
        """Parse the PDB file as a series of record entries into a dictionary"""
        
        with open(self.filename, 'r') as entries:
            model = 0
            open_model = False
            chain = None
            resnum = None
            for entry in entries:
                record_type = entry[0:6]

                if record_type == 'ATOM  ' or record_type == 'HETATM' and not entry[17:20] == 'HOH':
                    if not open_model:
                        model += 1
                        open_model = True
                        self.structure.append(dict())
                    if not entry[21] == chain:
                        chain = entry[21]
                        if chain == ' ':
                            chain = 'A'
                        if not chain in self.structure[model]:
                            self.structure[model][chain] = dict()
                    if not int(entry[22:26]) == resnum:
                        resnum = int(entry[22:26])
                        self.structure[model][chain][resnum] = []
                    self.structure[model][chain][resnum].append({'type': record_type,
                                                                 'serial': int(entry[6:11]),
                                                                 'name': entry[12:16].strip(),
                                                                 'altLoc': entry[16],
                                                                 'resname': entry[17:20].strip(),
                                                                 'icode': entry[26],
                                                                 'x': float(entry[30:38]),
                                                                 'y': float(entry[38:46]),
                                                                 'z': float(entry[46:54]),
                                                                 'occupancy': float(entry[54:60]),
                                                                 'bfactor': -1.00,
                                                                 'element': entry[76:78].strip()
                                                                 })
                
                elif record_type == 'MODEL ':
                    open_model = True
                    _model = int(entry.split()[1])
                    if not _model == model+1:
                        print('MODEL records should be sequentially numbered beginning with 1.')
                        print('MODEL {} renumbered to {}'.format(_model, model+1))
                    model += 1
                    self.structure.append(dict())
    
                elif record_type == 'ENDMDL':
                    open_model = False
                    chain = None
                    resnum = None
                
                elif record_type[:3] == 'TER':
                    chain = None
                    resnum = None

            self.chains = {chain:self._get_sequence(self.structure[1][chain]) for chain in self.structure[1]}
            self.num_chains = len(self.chains)
            self.num_models = len(self.structure)

    def _get_sequence(self, chain):
        """Parse single character amino acid code from chain"""
        _sequence = ''

        for residue in chain:
            resname = chain[residue][0]['resname']
            if chain[residue][0]['type'] == 'ATOM':
                try:
                    _sequence += self.AA_3to1[resname]
                except KeyError:
                    print('The residue name {} is not standard amino acid. Not added to FASTA sequence.'.format(resname))

        return _sequence
    
    def get_coords(self, model, chain, resnum, atom_name):
        """Return the xyz coordinates of specified atom"""
        for atom in self.structure[model][chain][resnum]:
            if atom['name'] == atom_name:
                return atom['x'], atom['y'], atom['z']
        print('No atom with name {} found for this residue'.format(atom_name))
        

#%% FUNCTIONS

def joules_to_kj_per_mol(value):
    """Converts values in J units to kJ/mol"""
    return value*0.001*AVOGADRO
    
def kj_per_mol_to_joules(value):
    """Converts values in kJ/mol units to J"""
    return value*1000/AVOGADRO

def distance_matrix(a, b):
    """
    Calculate the Euclidean distance between groups of atoms.
    Uses NumPy broadcasting logic to efficiently determine vectors
    before returning the Euclidean distance.
    Functionally equivalent to `scipy.spatial.distance_matrix` for p=2.
    
    Parameters
    ----------
    a : list or array of xyz coordinates (Nx3)
    b : list or array of xyz coordinates (Mx3)
    
    Returns
    -------
    distances : array of Euclidean distances (NxM)
    """
    
    _a = np.asarray(a)[:,np.newaxis,:]
    _b = np.asarray(b)[np.newaxis,:,:]
    return np.sum((_a - _b)**2, axis = -1)**0.5

def symmetric_matrix(vector):
    """
    Create NxN symmetric matrix from 1xN vector
    
    Parameters
    ----------
    vector : array (N)
    
    Returns
    -------
    matrix : symmetric matrix (NxN)
    """
    
    return vector*vector[:, np.newaxis]

def moveable_protons(deprotonated_charges, target_charge):
    """
    Calculate the number of moveable protons required to attain the `target_charge`.
    Performing `np.sum(proton_vector + deprotonated_charges)` gives `target_charge`.
    
    Parameters
    ----------
    deprotonated_charges : np.ndarray of deprotonated charges (N)
    target_charge : the final charge of the protein
    
    Returns
    -------
    proton_vector : a randomised guess for the `proton_sequence` vector (N)
    """
    num_protons = target_charge - deprotonated_charges.sum()
    proton_vector = np.zeros_like(deprotonated_charges)
    proton_vector[:num_protons] = 1
    np.random.shuffle(proton_vector)
    return proton_vector
    
def parse_coordinates(pdb_file, point_charges=point_charge_dict,
                      deprot_charges=deprot_charge_dict, 
                      proton_affinities=proton_affinity_dict):
    """Parses coordinates for chargeable residues from `pdb_file`.
    
    Parameters
    ----------
    
    Returns
    -------
    resn : list
    resi : list
    deprot_charges : ndarray
    xyz : ndarray
    affinities : ndarray
    """
    structure = PDB(pdb_file)
    
    # Initialise lists
    resn = []
    resi = []
    deprot_charges = []
    affinities = []
    xyz = []
    
    # Preselect model and chain in case of multiple of either
    # TODO: Add multi-chain and multi-model behaviour 
    model = 1
    chain = 'A'
    
    residues = sorted(structure.structure[model][chain])

    # Add N-terminus to lists
    resn.append('NT')
    resi.append(residues[0])
    deprot_charges.append(deprot_charge_dict['NT'])
    affinities.append(proton_affinity_dict['NT'])
    xyz.append(structure.get_coords(model, chain, residues[0], point_charge_dict['NT']))

    for residue in residues:
        resname = structure.structure[model][chain][residue][0]['resname']
        if resname in ['ASP', 'GLU', 'LYS', 'ARG', 'HIS']:
            resn.append(resname)
            resi.append(residue)
            deprot_charges.append(deprot_charge_dict[resname])
            affinities.append(proton_affinity_dict[resname])
            xyz.append(structure.get_coords(model, chain, residue, point_charge_dict[resname]))

    # Add C-terminus to lists
    resn.append('CT')
    resi.append(residues[-1])
    deprot_charges.append(deprot_charge_dict['CT'])
    affinities.append(proton_affinity_dict['CT'])
    xyz.append(structure.get_coords(model, chain, residues[-1], point_charge_dict['CT']))

    return resn, resi, np.array(deprot_charges), np.array(xyz), np.array(affinities)

def coulomb_energy(charge_seq, distances, mask):
    """Calculate Coulomb energy for charges by distance
    
    Parameters
    ----------
    charge_seq : ndarray
    distances : ndarray
    mask : ndarray
    
    Returns
    -------
    coulomb energy : float
       Coulomb energy (electrostatic energy) in Joules 
    """
    charge_mat = symmetric_matrix(charge_seq)
    return E_CONST * np.sum(charge_mat[mask] / distances)
    
def binding_energy(proton_seq, affinities):
    """Calculate binding energy from proton affinities
    
    Parameters
    ----------
    proton_seq : ndarray
        protonation state of residues
    affinities : ndarray
        proton affinities for residues
    
    Returns
    -------
    binding energy : float
        Binding energy in Joules
    """
    return affinities[proton_seq.nonzero()[0]].sum()

def minimise_energy(deprot_charges, affinities, xyz, charge, coulomb_only=False, verbose=True):
    """Minimising function. Iterates through a proton sequence to find the
       combination that provides the miminal energy.
    
    Parameters
    ----------
    deprot_charges : charge of residue when deprotonated (1xN array) 
    affinities : proton affinities of each residue (1xN array)
    xyz : coordinates for point charges (Nx3)
    charge : target charge state (int)
    coulomb_only : whether to only calculate Coulomb energy (boolean)
    verbose : whether to print results (boolean)
    
    Returns
    -------
    proton_seq : current best proton sequence after minimisation
    e_total : total energy of `proton_seq` after minimisation 
                (Only if `coulomb_only`=False)
    e_coulomb : Coulomb energy of `proton_seq` after minimisation
    e_proton : binding energy of `proton_seq` after minimisation
    """
    
    # Initialise local variables
    proton_seq = moveable_protons(deprot_charges, charge)
    mask = np.mask_indices(len(proton_seq), np.triu, 1)
    distances = distance_matrix(xyz, xyz)[mask]
    
    if coulomb_only:
        get_energy = lambda: coulomb_energy(proton_seq+deprot_charges, distances, mask)
    else:
        get_energy = lambda: (coulomb_energy(proton_seq+deprot_charges, distances, mask)
                              - binding_energy(proton_seq, affinities))
    
    # Initial energies
    current_min = get_energy()
    shunt_min = current_min
    counters = [time.process_time(), 0, 0]

    while (shunt_min <= current_min):
        counters[1] += 1
        if verbose:
            print('Shunt={}'.format(counters[1]))
        shunt_min = get_energy()
        best_shunt = [0, 0]
        deprot_sequence = np.where(proton_seq == 0)[0]
        
        for p in proton_seq.nonzero()[0]:
            proton_seq[p] = 0
            # For all protonatable sites
            for d in deprot_sequence:
                counters[2] += 1
                proton_seq[d] = 1
                e_tot = get_energy()
                if verbose:
                    print('Step {}, {:10.2f} kJ/mol'.format(counters[2], e_tot))
                if e_tot <= shunt_min:
                    shunt_min = e_tot
                    best_shunt = [p, d]
                proton_seq[d] = 0
            proton_seq[p] = 1
        if verbose:
            print('Shunt {} minimum energy {:.2f} kJ/mol'.format(counters[1], shunt_min))
        # Update `proton_seq` to best values
        if (shunt_min >= current_min):
            e_coulomb = coulomb_energy(proton_seq+deprot_charges, distances, mask)
            e_proton = binding_energy(proton_seq, affinities)
            counters[0] = time.process_time() - counters[0]
            print("Best Sequence\n-------------\n{}".format(proton_seq))
            print("Coulomb energy = {:.2f} kJ/mol".format(e_coulomb))
            if not coulomb_only:
                print("Binding energy = {:.2f} kJ/mol".format(e_proton))
                print("Total energy = {:.2f} kJ/mol".format(current_min))
            print("Optimisation completed in {:.2f} seconds after {} shunts in a total of {} steps.".format(*counters))
            break
        # Reset `proton_seq` to best sequence to reseed
        proton_seq[best_shunt[0]] = 0
        proton_seq[best_shunt[1]] = 1
        current_min = shunt_min
    
    return proton_seq, current_min, e_coulomb, e_proton
    
def alanine_scan(resn, resi, deprot_charges, affinities, xyz, charge, coulomb_only=False, verbose=True):
    """In silico alanine scanning of chargeable side-chains
    
    Parameters
    ----------

    
    Returns
    -------

    """
    # TODO: fill this function and refactor
    ignore_mask = np.ones_like(deprot_charges, dtype=bool)
    for r, res in enumerate(resn):
        if res in ['NT', 'CT']:
            continue
        print('\n{} {} -> ALA...'.format(res, resi[r]))
        ignore_mask[r] = False
        mutant_min_energy = minimise_energy(deprot_charges[ignore_mask], 
                                            affinities[ignore_mask], 
                                            xyz[ignore_mask], 
                                            charge, 
                                            coulomb_only, 
                                            verbose)
        ignore_mask[r] = True
    
def save_proton_sequence(filename, proton_sequence, e_coulomb, e_proton, resn, resi, pdb_file=None):
    """Save minimised proton sequence to file in the form:
    
    Parameters
    ----------
    filename : string
    proton_sequence : ndarray
        Proton sequence that gives a (reasonably) minimised energy.
    resn : list
    resi : list
    pdb_file : string
        Name of input PDB file from which proton_sequence was generated
    
    Returns
    -------
    outfile : file
        Text file containing all protonated residues and termini by row.
        <RESN>\t<RESI>\t<CHAIN>
    """
    with open(filename+'proton_sites.txt', 'w') as outfile:
        outfile.write('# This file was generated by chargePlacer.py, {}\n'.format(time.strftime("%d %b %Y %H:%M:%S")))
        if pdb_file:
            outfile.write('# From {}\n'.format(pdb_file))
        outfile.write('# Coulomb Energy = {:.2f} kJ/mol, Proton Binding Energy = {:.2f} kJ/mol\n'.format(e_coulomb, e_proton))
        for r, res in enumerate(resi):
            if proton_sequence[r]:
                outfile.write('{}\t{}\t{}\n'.format(resn[r], res, 'A'))
    print('Proton sequence successfully saved to {}'.format(filename+'proton_sites.txt'))



if __name__ == '__main__':
    #print(__doc__)
    argparser = argparse.ArgumentParser(description=__doc__, 
                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    argparser.add_argument('input', metavar='INPUT', help='input PDB file to determine charges for')
    argparser.add_argument('charge', metavar='CHARGE', help='target charge state', type=int)
    argparser.add_argument('-v', '--verbose', 
                           help='verbose output', action='store_true')
    argparser.add_argument('-c', '--coulomb_only', 
                           help='minimise for Coulomb repulsion only, ignores proton affinity', 
                           action='store_true')
    argparser.add_argument('-r', '--relative_permittivity', metavar='',
                           help='relative permittivity to use (default: 1)', 
                           default=1, type=float)
    argparser.add_argument('-a', '--alanine_scan', 
                           help='perform in silico alanine scanning for all chargeable residues', 
                           action='store_true')
    argparser.add_argument('-o', '--output', 
                           help='prefix for output files (default: ""). Gives *proton_sites.txt and *charges.txt', 
                           default='')
    args = argparser.parse_args()

    if args.relative_permittivity != 1:
        E_CONST = _E_CONST / args.relative_permittivity
    
    print('Opening {} and parsing coordinates...'.format(args.input))
    resn, resi , deprot_charges, xyz, affinities = parse_coordinates(args.input)
  
    print('Minimising energy of proton sequence...')
    min_energy = minimise_energy(deprot_charges, 
                                 affinities, 
                                 xyz, 
                                 args.charge,
                                 coulomb_only = args.coulomb_only, 
                                 verbose = args.verbose)
    
    if args.alanine_scan:
        print('Beginning in silico alanine scanning...')
        alanine_scan(resn,
                     resi,
                     deprot_charges, 
                     affinities, 
                     xyz, 
                     args.charge,
                     coulomb_only = args.coulomb_only, 
                     verbose = args.verbose)
    
    save_proton_sequence(args.output, 
                         min_energy[0], 
                         min_energy[2], 
                         min_energy[3], 
                         resn, 
                         resi, 
                         args.input)
