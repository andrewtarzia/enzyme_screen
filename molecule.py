#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
# ensure cpickle usage
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import cirpy
import glob
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
import PUBCHEM_IO
from rdkit_functions import calc_molecule_diameter
from numpy import average


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
        self.complexity = None
        self.XlogP = None
        self.pkl = self.get_pkl()
        self.DB_list = [DB]
        self.min_diam = None
        self.mid_diam = None
        self.max_diam = None
        self.rat_1 = None
        self.rat_2 = None

    def molecule_db_dir(self):
        return '/home/atarzia/psp/molecule_DBs/atarzia/'

    def molecule_db_prefix(self):
        return 'ATRS_'

    def determine_ID(self):
        """Determine an ID based on what is in my molecule_db_dir.

        """
        dir = self.molecule_db_dir()
        pre = self.molecule_db_prefix()

        existing_pkls = glob.glob(dir+pre+'*.pkl')
        existing_ids = [int(i.replace(dir+pre, '').replace('.pkl', ''))
                        for i in existing_pkls]
        if len(existing_ids) > 0:
            max_id = max(existing_ids)
        else:
            max_id = 0
        return str(max_id+1)

    def get_pkl(self):
        pkl = self.molecule_db_dir()+self.molecule_db_prefix()
        pkl += self.determine_ID()+'.pkl'
        return pkl

    def save_object(self, filename):
        """Pickle molecule object to file.

        """
        # Overwrites any existing file.
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(self, filename, verbose=True):
        """unPickle molecule object from file.

        """
        if verbose:
            print('loading:', filename)
        with open(filename, 'rb') as input:
            self = pickle.load(input)
            return self

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

    def get_properties(self, check=True):
        """Calculate some general molecule properties from SMILES

        From RDKit:
            - synthetic accesibility:
                https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3225829/
                (1 = easy to make, 10 = harder)
            - logP (hydrophobicity):
                https://pubs.acs.org/doi/10.1021/ci990307l
                (smaller = more hydrophilic)
        From PUBCHEM:
            - molecule complexity:
                https://pubchemdocs.ncbi.nlm.nih.gov/glossary$Complexity
                (0 to inf)
            - XlogP:
                https://pubchemdocs.ncbi.nlm.nih.gov/glossary$XLogP
                (smaller = more hydrophilic)
        """
        if check:
            try:
                if self.SMILES is not None:
                    rdkitmol = Chem.MolFromSmiles(self.SMILES)
                    rdkitmol.Compute2DCoords()
                    if self.logP is None:
                        self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
                    if self.Synth_score is None:
                        self.Synth_score = get_SynthA_score(rdkitmol)
                if self.XlogP is None:
                    if self.iupac_name is not None:
                        self.XlogP = PUBCHEM_IO.get_logP_from_name(self.iupac_name)
                    else:
                        self.XlogP = PUBCHEM_IO.get_logP_from_name(self.name)
                if self.complexity is None:
                    if self.iupac_name is not None:
                        self.complexity = PUBCHEM_IO.get_complexity_from_name(self.iupac_name)
                    else:
                        self.complexity = PUBCHEM_IO.get_complexity_from_name(self.name)
            except AttributeError:
                print('remove this!')
                self.XlogP = None
                self.complexity = None
        else:
            if self.SMILES is not None:
                rdkitmol = Chem.MolFromSmiles(self.SMILES)
                rdkitmol.Compute2DCoords()
                self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
                self.Synth_score = get_SynthA_score(rdkitmol)
            if self.iupac_name is not None:
                self.XlogP = PUBCHEM_IO.get_logP_from_name(self.iupac_name)
                self.complexity = PUBCHEM_IO.get_complexity_from_name(self.iupac_name)
            else:
                self.XlogP = PUBCHEM_IO.get_logP_from_name(self.name)
                self.complexity = PUBCHEM_IO.get_complexity_from_name(self.name)

    def cirpy_to_iupac(self):
        """Attempt to resolve IUPAC molecule name using CIRPY.

        Returns None if not possible.

        """
        self.iupac_name = cirpy.resolve(self.name, 'iupac_name')


def load_molecule(filename, verbose=True):
    """unPickle molecule object from file.

    """
    if verbose:
        print('loading:', filename)
    with open(filename, 'rb') as input:
        mol = pickle.load(input)
        return mol


def check_molecule_unique(molec, molecules):
    """Check if a molecule is unique compared to those in molecule DB.

    Returns the pkl file of the original molecule also.

    """


    return unique, original_pkl


def get_all_molecules_from_rxn_systems(rxns):
    """From list of reactions, collect all molecules into molecule DB.

    This is a one off function to update the molecule database because it was
    written after reaction system collection.

    """
    for rs in rxns:
        for m in rs.components:
            new_mol = molecule(name=m.name, role=m.role,
                               DB=m.DB, DB_ID=m.DB_ID)
            # check if unique
            molecules = glob.glob('/home/atarzia/psp/molecule_DBs/atarzia/ATRS_*.pkl')
            unique, old_pkl = check_molecule_unique(m, molecules)
            if unique is True:
                # copy old object properties to new but only overwrite None or NaN
                for key, val in m.__dict__.items():
                    if key not in new_mol.__dict__:
                        new_mol.__dict__[key] = val
                    elif new_mol.__dict__[key] is None and val is not None:
                        new_mol.__dict__[key] = val
                new_mol.save_object(new_mol.pkl)
            else:
                # check if different database

                # add database to list of DBs

                # save object
                new_mol.save_object(new_mol.pkl)
            break
        break


def get_SynthA_score(mol):
    """Get synthetic accesibility score from RDKIT contrib (SA_score).

    """
    from SA_score import sa_scores
    s = sa_scores.calculateScore(mol)
    return s


if __name__ == "__main__":
    print('testing and debugging')
    rs_dir = '/home/atarzia/psp/screening_results/biomin_search/'
    import rxn_syst
    get_all_molecules_from_rxn_systems(rxn_syst.yield_rxn_syst(rs_dir))

    test_mol = '/home/atarzia/psp/molecule_DBs/atarzia/ATRS_13.pkl'
    mol = load_molecule(test_mol)
    print(mol.mol)
