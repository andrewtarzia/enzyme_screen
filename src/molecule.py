#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import glob
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from numpy import average

from rdkit_functions import calc_molecule_diameter


class Molecule:
    """
    Class that defines the molecules extracted from a database.

    """

    def __init__(self, name, role, DB, DB_ID, params, structure_file):
        self.name = name
        if role == 'Substrate':
            self.role = 'reactant'
        else:
            self.role = role.lower()
        self.DB = DB
        self.DB_ID = DB_ID
        self.params = params
        self.structure_file = structure_file
        self.SMILES = None
        self.logP = None
        self.logS = None
        self.Synth_score = None
        self.DB_list = [DB]
        self.min_diam = None
        self.mid_diam = None
        self.max_diam = None
        self.rat_1 = None
        self.rat_2 = None

    def write_structure(self, mol):
        Chem.MolToMolFile(mol, self.structure_file)

    def read_structure_to_mol(self):
        return Chem.MolFromMolFile(self.structure_file)

    def read_structure_to_smiles(self):
        return Chem.MolToSmiles(
            Chem.MolFromMolFile(self.structure_file)
        )

    def get_properties(self, check=True):
        """
        Calculate some general molecule properties from SMILES

        From RDKit:
            - synthetic accesibility:
                https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3225829/
                (1 = easy to make, 10 = harder)
            - logP (hydrophobicity):
                https://pubs.acs.org/doi/10.1021/ci990307l
                (smaller = more hydrophilic)
            - logS (aqueous solubility):
                https://github.com/PatWalters/solubility
                (smaller = less water soluble)

        """

        print('collect molecular properties using RDKit.')
        # logP and SA from RDKIT with SMILES:
        rdkitmol = Chem.MolFromSmiles(self.SMILES)
        rdkitmol = Chem.AddHs(rdkitmol)
        rdkitmol.Compute2DCoords()
        self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
        self.logS = get_logSw(rdkitmol)
        self.Synth_score = get_SynthA_score(rdkitmol)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name}, role={self.role}, DB={self.DB}, '
            f'ID={self.DB_ID}, SMILEs:\n{self.SMILES})'
        )

    def __repr__(self):
        return str(self)


def get_logSw(mol):
    """
    Get water solubility using RDKit as described here:
    https://github.com/PatWalters/solubility.
    Using the newly paramterized function.

    """
    from solubility.esol import ESOLCalculator
    esol_calculator = ESOLCalculator()
    logS = esol_calculator.calc_esol(mol)
    return logS


def get_SynthA_score(mol):
    """
    Get synthetic accesibility score from RDKIT contrib (SA_score).

    """
    from SA_score import sa_scores
    s = sa_scores.calculateScore(mol)
    return s


def yield_molecules(directory, file=False):
    """
    Read molecule list from all pkl files in directory or from a
    text file.

    """
    if file is False:
        files = sorted(glob.glob(directory+'ATRS_*.gpkl'))
    else:
        files = []
        with open(file, 'r') as f:
            for line in f.readlines():
                files.append(line.rstrip())
    for f in files:
        # print('----')
        # print('doing', count, 'of', len(files))
        # print(f)
        # load in molecule
        mol = Molecule.load_molecule(f, verbose=False)
        yield mol


def populate_all_molecules(directory, vdwScale, boxMargin, spacing,
                           N_conformers, MW_thresh, mol_file=False):
    """
    Populate all molecules in pickle files in directory.

    """
    count = 0
    for mol in yield_molecules(directory=directory, file=mol_file):
        # K_count = 0
        # for R in mol.rs_pkls:
        #     if 'KEGG' in R:
        #         K_count += 1
        # if K_count == 0:
        #     continue
        print('--------------------------------------------------')
        print('populating:', mol.name, mol.pkl)
        print('----', count)
        count += 1
        print('--------------------------------------------------')
        # properties to get:
        # iupac name
        # if mol.iupac_name is None and mol.cirpy_done is False:
        #     print('getting IUPAC name from CIRPY...')
        #     mol.cirpy_to_iupac()
        #     if mol.iupac_name is None:
        #         mol.iupac_name
        # logP
        if mol.logP is None:
            print('getting logP from RDKit...')
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            if rdkitmol is None:
                mol.logP = 'not found'
                continue
            rdkitmol.Compute2DCoords()
            mol.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
        # water solubility
        if mol.logS is None:
            print('getting logS from RDKit...')
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            if rdkitmol is None:
                mol.logS = 'not found'
                continue
            rdkitmol.Compute2DCoords()
            mol.logS = get_logSw(rdkitmol)
        # synthetic accessibility
        if mol.Synth_score is None:
            print('getting synthetic accessibility from RDKit...')
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            if rdkitmol is None:
                mol.Synth_score = 'not found'
                continue
            rdkitmol.Compute2DCoords()
            mol.Synth_score = get_SynthA_score(rdkitmol)
        # # check if compound is charged -
        # if so, set logS, logP and XlogP
        # if Chem.GetFormalCharge(Chem.MolFromSmiles(mol.SMILES)) != 0:
        #     print('setting charged compounds')
        #     mol.logP = 'charged'
        #     mol.logS = 'charged'
        #     mol.XlogP = 'charged'
        # diameters and ratios
        # assume if one is None then all are None!
        if mol.mid_diam is None:
            print('doing size calculation...')
            # name: smiles
            res = calc_molecule_diameter(mol.name, mol.SMILES,
                                         out_dir=directory,
                                         vdwScale=vdwScale,
                                         boxMargin=boxMargin,
                                         spacing=spacing,
                                         N_conformers=N_conformers,
                                         MW_thresh=MW_thresh)
            if res is None or len(res) == 0:
                # get the min values of all diameters of all conformers
                mol.min_diam = 0
                mol.mid_diam = 0
                mol.max_diam = 0
                # get avg values of all ratios of all conformers
                mol.rat_1 = 0
                mol.rat_2 = 0
            else:
                # get the min values of all diameters of all conformers
                mol.min_diam = round(min(res['diam1']), 3)
                mol.mid_diam = round(min(res['diam2']), 3)
                mol.max_diam = round(min(res['diam3']), 3)
                # get avg values of all ratios of all conformers
                mol.rat_1 = round(average(res['ratio_1']), 3)
                mol.rat_2 = round(average(res['ratio_2']), 3)
        print('--------------------------------------------------')
        # save object
        mol.save_object(mol.pkl)
