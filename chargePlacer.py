#!/usr/bin/env python3
#-*- coding:utf-8 -*-

"""
A python implementation of a charge positioning algorithm for gas-phase
molecular dynamics.
"""

# Copyright 2020 Jedd Bellamy-Carter
# MIT License
from __future__ import division, print_function, absolute_import

import argparse
import os
import sys
import time
import numpy as np

# GLOBAL VARIABLES

_E_CONST = 2.307E-18  # J.Å
relative_permittivity = 1
E_CONST = _E_CONST / relative_permittivity  #Relative permittivity of water
AVOGADRO = 6.022E+23

# LOOKUP DICTIONARIES

# Atoms to get coordinates for point-charge assignment
point_charge_dict = {'ASP': 'OD2',
                     'GLU': 'OE2',
                     'LYS': 'NZ', 
                     'ARG': 'NH2', 
                     'HIS': 'CB', 
                     'NT': 'N', 
                     'CT': 'C'}
# Charge of the deprotonated residue
deprot_charge_dict = {'ASP': -1, 
                      'GLU': -1, 
                      'LYS': 0, 
                      'ARG': 0, 
                      'HIS': 0, 
                      'NT': 0, 
                      'CT': -1}
# Proton affinities in Joules
proton_affinity_dict = {'ASP': 2.414e-18,
                        'GLU': 2.405e-18,
                        'LYS': 1.524e-18, 
                        'ARG': 1.664e-18, 
                        'HIS': 1.591e-18, 
                        'NT': 1.472e-18, 
                        'CT': 2.375e-18}

# CLASSES

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

        self.AA_3to1    = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E',
                           'GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K',
                           'MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
                           'TYR':'Y','VAL':'V'}

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
        

# FUNCTIONS

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

def sym_mat(vector):
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


def print_energy(charge_product):
    """Prints formatted energy returned from calcualtion, in kJ/mol"""
    print("{:G} kJ/mol".format(charge_product*0.001*AVOGADRO))
    return
    
def parse_coordinates(pdb_file, point_charges=point_charge_dict,
                      deprot_charges=deprot_charge_dict, 
                      proton_affinities=proton_affinity_dict):
    """Parses coordinates for chargeable residues from `pdb_file`.
    
    Parameters
    ----------
    
    Returns
    -------
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

    return np.array(deprot_charges), np.array(xyz), np.array(affinities)

def matrix_impl(charge_seq, dist_m, mask):
    charge_mat = sym_mat(charge_seq)
    return E_CONST*np.sum(charge_mat[mask]/dist_m)
    
def minimise_energy(proton_sequence, deprot_charges, affinities, dist_m, mask, correct_affinity=True, verbose=True, shuffle=False):
    """Minimising function. Iterates through a proton sequence to find the
       combination that provides the miminal energy.
    
    Parameters
    ------
    proton_sequence : initial proton sequence (1xN array)
    deprot_charges : charge of residue when deprotonated (1xN array) 
    affinities : proton affinities of each residue (1xN array)
    verbose : whether to print results (boolean)
    shuffle : if True, the input `proton_sequence` will be shuffled before starting
    
    Returns
    -------
    current_seq : current best proton sequence after minimisation
    current_min : total energy of `current_seq` after minimisation
    e_coulomb : Coulomb energy of `current_seq` after minimisation
    e_proton : binding energy of `current_seq` after minimisation
    """
    
    if not correct_affinity:
        affinities = np.zeros_like(proton_sequence)
    if shuffle:
        current_seq = np.random.permutation(proton_sequence)
    else:
        current_seq = proton_sequence.copy()
    current_min = matrix_impl(current_seq+deprot_charges, dist_m, mask) - affinities[current_seq.nonzero()[0]].sum()
    best_seqs = [[0],[current_min*0.001*AVOGADRO],[current_seq.copy()]] #initialise with starting sequence
    if verbose:
        print("Starting sequence\n------------------------------\n{}".format(current_seq))
        print("Coulomb energy = {:.0f} kJ/mol".format(matrix_impl(current_seq+deprot_charges, dist_m, mask)*0.001*AVOGADRO))
        print("Binding energy = {:.0f} kJ/mol".format(affinities[current_seq.nonzero()[0]].sum()*0.001*AVOGADRO))
        print("Total energy = {:.0f} kJ/mol".format(current_min*0.001*AVOGADRO))
    shunt_min = current_min
    counters = [time.process_time(), 0, 0]

    while (shunt_min <= current_min):
        counters[1] += 1
        shunt_min = matrix_impl(current_seq+deprot_charges, dist_m, mask) - affinities[current_seq.nonzero()[0]].sum()
        best_shunt = [0, 0]
        deprot_sequence = np.where(current_seq == 0)[0]
        for p in current_seq.nonzero()[0]:
            current_seq[p] = 0
            # For all protonatable sites
            for d in deprot_sequence:
                counters[2] += 1
                current_seq[d] = 1
                e_coulomb = matrix_impl(current_seq+deprot_charges, dist_m, mask)
                e_proton = affinities[current_seq.nonzero()[0]].sum()
                e_tot = e_coulomb - e_proton
                if e_tot <= shunt_min:
                    #print(current_seq)
                    best_seqs[0].append(counters[2])
                    best_seqs[1].append(e_tot*0.001*AVOGADRO)
                    best_seqs[2].append(current_seq.copy())
                    shunt_min = e_tot
                    best_shunt = [p, d]
                current_seq[d] = 0
            current_seq[p] = 1

        # Update `current_seq` to best values
        if (shunt_min >= current_min):
            e_coulomb = matrix_impl(current_seq+deprot_charges, dist_m, mask)*0.001*AVOGADRO
            e_proton = affinities[current_seq.nonzero()[0]].sum()*0.001*AVOGADRO
            if verbose:
                counters[0] = time.process_time() - counters[0]
                print("\nBest Sequence\n------------------------------\n{}".format(current_seq))
                print("Coulomb energy = {:.0f} kJ/mol".format(e_coulomb))
                print("Binding energy = {:.0f} kJ/mol".format(e_proton))
                print("Total energy = {:.0f} kJ/mol".format(current_min*0.001*AVOGADRO))
                print("\nOptimisation completed in {:.2f} seconds after {} shunts in a total of {} steps.".format(*counters))
            break
        current_seq[best_shunt[0]] = 0
        current_seq[best_shunt[1]] = 1
        #best_seqs.append(current_seq.copy())
        #best_seqs.append(shunt_min)
        current_min = shunt_min
    
    return current_seq, current_min*0.001*AVOGADRO, e_coulomb, e_proton, best_seqs

# Importing Files
# TODO: Add command line functionality
pdb_file = '/home/jedd/Documents/PostDoc/MD/protein_init.pdb'

deprot_charges, xyz, affinities = parse_coordinates(pdb_file)

# Generate for 8+
target_charge = 8
proton_seq = moveable_protons(deprot_charges, target_charge)
# Mask for arrays (speeds up calculation)
mask = np.mask_indices(len(proton_seq), np.triu, 1)
distance_mat = distance_matrix(xyz, xyz)  # Create distance matrix. r_ij
dist_m = distance_mat[mask]

min_e1 = minimise_energy(proton_seq, deprot_charges, affinities, dist_m, mask, shuffle=True)

