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
            self.SMILES2MOL()

    def get_compound(self):
        """Get compound SMILES from available identifiers.

        """
        # check if molecule exists in molecule database already
        old_pkl = search_molecule_by_ident(self)
        if old_pkl is not None:
            DB = self.DB
            # load old pkl
            self = load_molecule(old_pkl, verbose=True)
            print('collected from personal DB!')
            if DB not in self.DB_list:
                self.DB_list.append(DB)
        elif self.DB == 'SABIO':
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
            # name is set to KEGG C-ID at this point
            self.change_name = True
            get_cmpd_information(self)
            if self.SMILES is None and self.mol is None:
                self.PUBCHEM_last_shot()
        elif self.DB == 'BKMS':
            if self.SMILES is not None:
                self.SMILES2MOL()
            else:
                from CHEBI_IO import get_cmpd_information
                # set DB specific properties
                self.chebiID = self.DB_ID
                get_cmpd_information(self)
                if self.SMILES is None and self.mol is None:
                    self.PUBCHEM_last_shot()
        elif self.DB == 'BRENDA':
            if self.SMILES is not None:
                self.SMILES2MOL()
            else:
                from CHEBI_IO import get_cmpd_information
                # set DB specific properties
                self.chebiID = self.DB_ID
                get_cmpd_information(self)
                if self.SMILES is None and self.mol is None:
                    self.PUBCHEM_last_shot()

    def SMILES2MOL(self):
        # assigned a PUBCHEM SMILES and IUPAC name
        rdkitmol = Chem.MolFromSmiles(self.SMILES)
        rdkitmol.Compute2DCoords()
        # remove molecules with generalised atoms
        if '*' in self.SMILES:
            self.mol = None
        else:
            self.mol = rdkitmol

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


def get_SynthA_score(mol):
    """Get synthetic accesibility score from RDKIT contrib (SA_score).

    """
    from SA_score import sa_scores
    s = sa_scores.calculateScore(mol)
    return s


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

    1 - check for same SMIlEs


    Returns the pkl file of the original molecule also.

    """
    if molec.SMILES is not None:
        for original_pkl in molecules:
            o_mol = load_molecule(original_pkl, verbose=False)
            if o_mol.SMILES == molec.SMILES:
                unique = False
                return unique, original_pkl
    return True, None


def search_molecule_by_ident(molec, molecules):
    """Search for molecule in molecule database.

    1 - check for same SMIlEs if SMILEs is not None
    2 - check for same IUPAC name if not None
    3 - check for same name
    4 - check for same DB and DB_ID
    5 - check for same KEGG_ID

    Returns the pkl file of the original molecule also.

    """
    for o_pkl in molecules:
        o_mol = load_molecule(o_pkl, verbose=False)
        if molec.SMILES is not None:
            if o_mol.SMILES == molec.SMILES:
                return o_pkl
        elif molec.iupac_name is not None:
            if o_mol.iupac_name == molec.iupac_name:
                return o_pkl
        elif o_mol.name == molec.name:
            return o_pkl
        elif o_mol.iupac_name == molec.name:
            return o_pkl
        elif o_mol.DB == molec.DB:
            if o_mol.DB_ID == molec.DB_ID:
                return o_pkl
        else:
            try:
                if o_mol.KEGG_ID == molec.KEGG_ID:
                    return o_pkl
            except AttributeError:
                pass
    return None


def get_all_molecules_from_rxn_systems(rxns):
    """From list of reactions, collect all molecules into molecule DB.

    This is a one off function to update the molecule database because it was
    written after reaction system collection.

    """
    count = 0
    for rs in rxns:
        print(rs.pkl, '-----', count)
        count += 1
        if rs.components is None:
            continue
        for m in rs.components:
            if m.SMILES is None:
                # we do not want a pkl file for all molecules without SMILES
                continue
            print(m.name)
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
                # we do not change the new molecule, but we update the old mol
                old_mol = load_molecule(old_pkl, verbose=False)
                print(old_pkl)
                # check if different database
                # add database to list of DBs
                for i in new_mol.DB_list:
                    if i not in old_mol.DB_list:
                        old_mol.DB_list.append(i)
                # add any attributes to the old mol that the new mol has
                for key, val in new_mol.__dict__.items():
                    if key not in old_mol.__dict__:
                        old_mol.__dict__[key] = val
                        print('added:', key)
                    elif old_mol.__dict__[key] is None and val is not None:
                        old_mol.__dict__[key] = val
                        print('added:', key)
                # save object
                old_mol.save_object(old_mol.pkl)


def yield_molecules(directory):
    files = glob.glob(directory+'ATRS_*.pkl')
    for f in files:
        print(f)
        # load in molecule
        mol = load_molecule(f, verbose=False)
        yield mol


def populate_all_molecules(directory, vdwScale, boxMargin, spacing,
                           N_conformers, MW_thresh):
    """Populate all molecules in pickle files in directory.

    """
    for mol in yield_molecules(directory=directory):
        # properties to get:
        # iupac name
        if mol.iupac_name is None:
            mol.cirpy_to_iupac()
        # logP
        if mol.logP is None:
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            rdkitmol.Compute2DCoords()
            mol.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
        # XlogP
        if mol.XlogP is None:
            if mol.iupac_name is not None:
                mol.XlogP = PUBCHEM_IO.get_logP_from_name(mol.iupac_name)
            else:
                mol.XlogP = PUBCHEM_IO.get_logP_from_name(mol.name)
        # synthetic accessibility
        if mol.Synth_score is None:
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            rdkitmol.Compute2DCoords()
            mol.Synth_score = get_SynthA_score(rdkitmol)
        # complexity
        if mol.complexity is None:
            if mol.iupac_name is not None:
                mol.complexity = PUBCHEM_IO.get_complexity_from_name(mol.iupac_name)
            else:
                mol.complexity = PUBCHEM_IO.get_complexity_from_name(mol.name)
        # diameters and ratios
        if mol.mid_diam is None:  # assume if one is None then all are None!
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

        # save object
        mol.save_object(mol.pkl)


if __name__ == "__main__":
    import rxn_syst
    import os
    import sys

    if (not len(sys.argv) == 4):
        print('Usage: molecule.py get_mol pop_mol update_KEGG\n')
        print("""    get_mol: T for overwrite and collection of molecules from RS
        in current dir (F for skip)
        ---- this function is useful if you update the base attributes of the molecule class.
        """)
        print('    pop_mol: T to run population of molecule properties (does not overwrite).')
        print('    update_KEGG: T to update KEGG translator.')
        sys.exit()
    else:
        get_mol = sys.argv[1]
        pop_mol = sys.argv[2]
        update_KEGG = sys.argv[3]

    if get_mol == 'T':
        print('extract all molecules from reaction syetms in current dir...')
        curr_dir = os.getcwd()
        print(curr_dir)
        get_all_molecules_from_rxn_systems(rxn_syst.yield_rxn_syst(curr_dir+'/'))

    if pop_mol == 'T':
        vdwScale = 0.8
        boxMargin = 4.0
        spacing = 0.6
        N_conformers = 50
        MW_thresh = 2000
        print('settings:')
        print('    VDW scale:', vdwScale)
        print('    Box Margin:', boxMargin, 'Angstrom')
        print('    Grid spacing:', spacing, 'Angstrom')
        print('    No Conformers:', N_conformers)
        print('    MW threshold:', MW_thresh, 'g/mol')
        inp = input('happy with these? (T/F)')
        if inp == 'F':
            sys.exit('change them in the source code')
        elif inp != 'T':
            sys.exit('I dont understand, T or F?')
        print('populate the properties attributes for all molecules in DB...')
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        populate_all_molecules(directory=directory,
                               vdwScale=vdwScale,
                               boxMargin=boxMargin,
                               spacing=spacing,
                               N_conformers=N_conformers,
                               MW_thresh=MW_thresh)

    if update_KEGG == 'T':
        translator = '/home/atarzia/psp/molecule_DBs/KEGG/translator.txt'
        # iterate over all molecules in DB and if they have a KEGG ID then
        # write translation to CHEBI ID
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        for mol in yield_molecules(directory=directory):
            if 'KEGG' in mol.DB_list:
                KID = mol.KEGG_ID
                pkl = mol.pkl
                with open(translator, 'a') as f:
                    f.write(KID+'__'+pkl+'\n')

    # for i in yield_molecules(directory=directory):
    #     print(i.DB_list)
    #     if 'KEGG' in i.DB_list:
    #         print(i.KEGG_ID)

    # a = '/home/atarzia/psp/molecule_DBs/atarzia/ATRS_9.pkl'
    # b = load_molecule(a)
    # b.name
    # print(b.SMILES)
    # b.KEGG_ID
