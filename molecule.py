#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import cirpy
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
import PUBCHEM_IO


class molecule:
    """Class that defines the molecules extracted from a database.

    """

    def __init__(self, name, role, DB, DB_ID):
        self.name = name
        if role == 'Substrate':
            self.role = 'reactant'
        else:
            self.role = role.lower()
        self.DB = DB
        self.DB_ID = DB_ID
        self.InChi = None
        self.iupac_name = None
        self.mid_diam = None
        self.SMILES = None
        self.logP = None
        self.Synth_score = None

    def PUBCHEM_last_shot(self):
        """Use PUBCHEM search for last chance at getting structure information.

        """
        # check for pubchem entry based on name
        print('search PUBCHEM')
        smiles = PUBCHEM_IO.get_SMILES_from_name(self.name)
        if smiles is not None:
            self.SMILES = smiles
            # assigned a PUBCHEM SMILES and IUPAC name
            rdkitmol = Chem.MolFromSmiles(self.SMILES)
            rdkitmol.Compute2DCoords()
            # remove molecules with generalised atoms
            if '*' in self.SMILES:
                self.mol = None
            else:
                self.mol = rdkitmol

    def get_compound(self):
        """Get reaction system from SABIO reaction ID (rID).

        """
        if self.DB == 'SABIO':
            from SABIO_IO import get_cmpd_information
            # set DB specific properties
            self.cID = self.DB_ID
            get_cmpd_information(self)
            if self.SMILES is None and self.mol is None:
                self.PUBCHEM_last_shot()
        elif self.DB == 'KEGG':
            from CHEBI_IO import get_cmpd_information
            # set DB specific properties
            self.chebiID = self.DB_ID
            self.change_name = True  # name is set to KEGG C-ID at this point
            get_cmpd_information(self)
            if self.SMILES is None and self.mol is None:
                self.PUBCHEM_last_shot()
        elif self.DB == 'BKMS':
            if self.SMILES is not None:
                # assigned a PUBCHEM SMILES and IUPAC name
                rdkitmol = Chem.MolFromSmiles(self.SMILES)
                rdkitmol.Compute2DCoords()
                # remove molecules with generalised atoms
                if '*' in self.SMILES:
                    self.mol = None
                else:
                    self.mol = rdkitmol
            else:
                from CHEBI_IO import get_cmpd_information
                # set DB specific properties
                self.chebiID = self.DB_ID
                get_cmpd_information(self)
                if self.SMILES is None and self.mol is None:
                    self.PUBCHEM_last_shot()
        elif self.DB == 'BRENDA':
            if self.SMILES is not None:
                # assigned a PUBCHEM SMILES and IUPAC name
                rdkitmol = Chem.MolFromSmiles(self.SMILES)
                rdkitmol.Compute2DCoords()
                # remove molecules with generalised atoms
                if '*' in self.SMILES:
                    self.mol = None
                else:
                    self.mol = rdkitmol
            else:
                from CHEBI_IO import get_cmpd_information
                # set DB specific properties
                self.chebiID = self.DB_ID
                get_cmpd_information(self)
                if self.SMILES is None and self.mol is None:
                    self.PUBCHEM_last_shot()

    def get_properties(self):
        """Calculate some general molecule properties from SMILES

        From RDKit:
            - synthetic accesibility:
                https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3225829/
                (1 = easy to make, 10 = harder)
            - logP (hydrophobicity):
                https://pubs.acs.org/doi/10.1021/ci990307l
                (smaller = more hydrophilic)
        """
        if self.SMILES is not None:
            rdkitmol = Chem.MolFromSmiles(self.SMILES)
            rdkitmol.Compute2DCoords()
            self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
            self.Synth_score = molecule.get_SynthA_score(rdkitmol)

    def cirpy_to_iupac(self):
        """Attempt to resolve IUPAC molecule name using CIRPY.

        Returns None if not possible.

        """
        self.iupac_name = cirpy.resolve(self.name, 'iupac_name')


def get_SynthA_score(mol):
    """Get synthetic accesibility score from RDKIT contrib (SA_score).

    """
    from SA_score import sa_scores
    s = sa_scores.calculateScore(mol)
    return s
