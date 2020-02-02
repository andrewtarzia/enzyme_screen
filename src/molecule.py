#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os.path import exists, join
import json
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors

import rdkit_functions as rdkf


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

    def read_prop_file(self):
        """
        Read properties file.

        """
        name = self.name
        prop_file = join(
            self.params['molec_dir'],
            name+'_prop.json'
        )
        if not exists(prop_file):
            raise FileNotFoundError(f'{prop_file} not found')

        with open(prop_file, 'r') as f:
            prop_dict = json.load(f)

        return prop_dict

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

        print('>>> collect molecular properties using RDKit.')
        # logP and SA from RDKIT with SMILES:
        rdkitmol = Chem.MolFromSmiles(self.SMILES)
        rdkitmol = Chem.AddHs(rdkitmol)
        rdkitmol.Compute2DCoords()
        self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
        self.logS = rdkf.get_logSw(rdkitmol)
        self.Synth_score = rdkf.get_SynthA_score(rdkitmol)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name}, role={self.role}, DB={self.DB}, '
            f'ID={self.DB_ID})'
        )

    def __repr__(self):
        return str(self)
