#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import pickle
import gzip
import cirpy
import pandas as pd
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

        existing_pkls = glob.glob(dir+pre+'*.gpkl')
        existing_ids = [int(i.replace(dir+pre, '').replace('.gpkl', ''))
                        for i in existing_pkls]
        if len(existing_ids) > 0:
            max_id = max(existing_ids)
        else:
            max_id = 0
        return str(max_id+1)

    def get_pkl(self):
        pkl = self.molecule_db_dir()+self.molecule_db_prefix()
        pkl += self.determine_ID()+'.gpkl'
        return pkl

    def save_object(self, filename):
        """Pickle molecule object to file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        # Overwrites any existing file.
        with gzip.GzipFile(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(self, filename, verbose=True):
        """unPickle molecule object from file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        if verbose:
            print('loading:', filename)
        with gzip.GzipFile(filename, 'rb') as input:
            self = pickle.load(input)
            return self

    def PUBCHEM_last_shot(self):
        """Use PUBCHEM search for last chance at getting structure information.

        """
        # check for pubchem entry based on name
        print('search PUBCHEM...')
        # smiles = PUBCHEM_IO.get_SMILES_from_name(self.name)
        smiles = PUBCHEM_IO.hier_name_search(self, 'CanonicalSMILES')
        if smiles is not None:
            self.SMILES = smiles
            self.SMILES2MOL()

    def get_compound(self, dataset, search_mol=True):
        """Get compound SMILES from available identifiers.

        """
        print('collecting information for:', self.name)
        # check if molecule exists in molecule database already
        old_pkl = None
        if search_mol:
            old_pkl = search_molecule_by_ident(self, dataset)
        if old_pkl is not None:
            DB = self.DB
            # load old pkl
            print('collected', self.name, 'from personal DB:')
            old_mol = load_molecule(old_pkl, verbose=True)
            # copy old_mol DB object properties to self
            # only overwrite None or NaN
            for key, val in old_mol.__dict__.items():
                if key not in self.__dict__:
                    self.__dict__[key] = val
                elif self.__dict__[key] is None and val is not None:
                    self.__dict__[key] = val
            if DB not in self.DB_list:
                self.DB_list.append(DB)
            return self
        elif self.DB == 'SABIO':
            # from SABIO_IO import get_cmpd_information
            # # set DB specific properties
            # self.cID = self.DB_ID
            # get_cmpd_information(self)
            # if self.SMILES is None and self.mol is None:
            #     self.PUBCHEM_last_shot()
            # return self
            from CHEBI_IO import get_cmpd_information
            get_cmpd_information(self)
            if self.SMILES is None and self.mol is None:
                self.PUBCHEM_last_shot()
            return self
        elif self.DB == 'KEGG':
            from CHEBI_IO import get_cmpd_information
            # set DB specific properties
            self.chebiID = self.DB_ID
            # name is set to KEGG C-ID at this point
            self.change_name = True
            get_cmpd_information(self)
            if self.SMILES is None and self.mol is None:
                self.PUBCHEM_last_shot()
            return self
        elif self.DB == 'BKMS':
            print(self.SMILES)
            if self.SMILES is not None:
                print('have SMILES already')
                self.SMILES2MOL()
                return self
            else:
                from CHEBI_IO import get_cmpd_information
                # set DB specific properties
                self.chebiID = self.DB_ID
                get_cmpd_information(self)
                if self.SMILES is None and self.mol is None:
                    self.PUBCHEM_last_shot()
                return self
        elif self.DB == 'BRENDA':
            if self.SMILES is not None:
                self.SMILES2MOL()
                return self
            else:
                from CHEBI_IO import get_cmpd_information
                # set DB specific properties
                self.chebiID = self.DB_ID
                get_cmpd_information(self)
                if self.SMILES is None and self.mol is None:
                    self.PUBCHEM_last_shot()
                return self

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
        print(self.SMILES)
        print(self.__dict__)
        if self.SMILES is not None:
            # logP and SA from RDKIT with SMILES:
            rdkitmol = Chem.MolFromSmiles(self.SMILES)
            rdkitmol.Compute2DCoords()
            self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
            self.Synth_score = get_SynthA_score(rdkitmol)
            # XlogP and complexity from PUBCHEM with SMILES:
            # check if already calculated
            # set NAME based on available DB IDs
            print(self.DB)
            print(self.DB_ID)
            if check:
                if self.XlogP is None:
                    self.XlogP = PUBCHEM_IO.hier_name_search(self, 'XLogP')
                if self.complexity is None:
                    self.complexity = PUBCHEM_IO.hier_name_search(self, 'complexity')
            else:
                self.XlogP = PUBCHEM_IO.hier_name_search(self, 'XLogP')
                self.complexity = PUBCHEM_IO.hier_name_search(self, 'complexity')
            print('check', check)
            print('xlogp', self.XlogP)
            print('comp', self.complexity)

    def cirpy_to_iupac(self):
        """Attempt to resolve IUPAC molecule name using CIRPY.

        Returns None if not possible.

        """
        cirpy_res = cirpy.resolve(self.name, 'iupac_name')
        if type(cirpy_res) is list:
            self.iupac_name = cirpy_res[0]
        else:
            self.iupac_name = cirpy_res
        self.cirpy_done = True


def get_SynthA_score(mol):
    """Get synthetic accesibility score from RDKIT contrib (SA_score).

    """
    from SA_score import sa_scores
    s = sa_scores.calculateScore(mol)
    return s


def load_molecule(filename, verbose=True):
    """unPickle molecule object from file.

    """
    filename = filename.replace('.pkl', '.gpkl')
    filename = filename.replace('.bpkl', '.gpkl')
    if verbose:
        print('loading:', filename)
    with gzip.GzipFile(filename, 'rb') as input:
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


def search_molecule_by_ident(molec, dataset):
    """Search for molecule in molecule database using lookup file.

    1 - check for same SMILEs if SMILEs is not None
    2 - check for same IUPAC name if not None
    3 - check for same name
    4 - check for same DB and DB_ID
    5 - check for same KEGG_ID
    6 - check for same InChiKey
    7 - check for same chebiID

    Returns the pkl file of the original molecule also.

    """
    for idx, row in dataset.iterrows():
        if molec.SMILES is not None:
            if row['SMILES'] == molec.SMILES:
                return row['pkl']
        elif molec.iupac_name is not None:
            if row['iupac'] == molec.iupac_name:
                return row['pkl']
        elif row['name'] == molec.name:
            return row['pkl']
        elif row['iupac'] == molec.name:
            return row['pkl']
        elif row['DB'] == molec.DB:
            if row['DB_ID'] == molec.DB_ID:
                return row['pkl']
        else:
            try:
                if row['KEGG_ID'] == molec.KEGG_ID:
                    return row['pkl']
            except AttributeError:
                pass
            try:
                if row['InChiKey'] == molec.InChiKey:
                    return row['pkl']
            except AttributeError:
                pass
            try:
                if row['CHEBI_ID'] == molec.chebiID:
                    return row['pkl']
            except AttributeError:
                pass
    return None


def get_all_molecules_from_rxn_systems(rxns, done_file, from_scratch='F'):
    """From list of reactions, collect all molecules into molecule DB.

    This is a one off function to update the molecule database because it was
    written after reaction system collection.

    """
    if scratch == 'T':
        with open(done_file, 'w') as f:
            f.write('pkls\n')
    count = 0
    for rs in rxns:
        already_done = False
        if scratch == 'F':
            with open(done_file, 'r') as f:
                for line in f.readlines():
                    if rs.pkl in line.rstrip():
                        already_done = True
        if already_done:
            continue
        print(rs.pkl, '-----', count)
        count += 1
        if rs.components is None:
            continue
        for m in rs.components:
            if m.SMILES is None:
                # we do not want a pkl file for all molecules without SMILES
                continue
            new_mol = molecule(name=m.name, role=m.role,
                               DB=m.DB, DB_ID=m.DB_ID)
            # check if unique
            molecules = glob.glob('/home/atarzia/psp/molecule_DBs/atarzia/ATRS_*.gpkl')
            unique, old_pkl = check_molecule_unique(m, molecules)
            print(m.name, 'u', unique)
            if unique is True:
                # copy old object properties to new but only overwrite None or NaN
                for key, val in m.__dict__.items():
                    if key not in new_mol.__dict__:
                        new_mol.__dict__[key] = val
                    elif new_mol.__dict__[key] is None and val is not None:
                        new_mol.__dict__[key] = val
                # add rxn syst pkl name
                try:
                    if rs.pkl not in new_mol.rs_pkls:
                        new_mol.rs_pkls.append(rs.pkl)
                except AttributeError:
                    new_mol.rs_pkls = []
                    new_mol.rs_pkls.append(rs.pkl)
                new_mol.save_object(new_mol.pkl)
            else:
                # we do not change the new molecule, but we update the old mol
                old_mol = load_molecule(old_pkl, verbose=False)
                print(old_pkl)
                # check if different database
                # add database to list of DBs
                try:
                    for i in m.DB_list:
                        if i not in old_mol.DB_list:
                            old_mol.DB_list.append(i)
                except AttributeError:
                    m.DB_list = [m.DB]
                    rs.save_object(rs.pkl)
                    for i in m.DB_list:
                        if i not in old_mol.DB_list:
                            old_mol.DB_list.append(i)
                # add any attributes to the old mol that the new mol has
                for key, val in m.__dict__.items():
                    if key not in old_mol.__dict__:
                        old_mol.__dict__[key] = val
                        print('added:', key)
                    elif old_mol.__dict__[key] is None and val is not None:
                        old_mol.__dict__[key] = val
                        print('added:', key)
                # add rxn syst pkl name
                try:
                    if rs.pkl not in old_mol.rs_pkls:
                        old_mol.rs_pkls.append(rs.pkl)
                except AttributeError:
                    old_mol.rs_pkls = []
                    old_mol.rs_pkls.append(rs.pkl)
                # save object
                old_mol.save_object(old_mol.pkl)
        # add rs.pkl to done_file
        with open(done_file, 'a') as f:
            f.write(rs.pkl+'\n')


def yield_molecules(directory, file=False):
    """Read molecule list from all pkl files in directory or from a text file.

    """
    if file is False:
        files = glob.glob(directory+'ATRS_*.gpkl')
    else:
        files = []
        with open(file, 'r') as f:
            for line in f.readlines():
                files.append(line.rstrip())
    count = 1
    for f in files:
        # print('----')
        # print('doing', count, 'of', len(files))
        # print(f)
        # load in molecule
        mol = load_molecule(f, verbose=False)
        count += 1
        yield mol


def populate_all_molecules(directory, vdwScale, boxMargin, spacing,
                           N_conformers, MW_thresh, mol_file=False):
    """Populate all molecules in pickle files in directory.

    """
    for mol in yield_molecules(directory=directory, file=mol_file):
        # properties to get:
        # iupac name
        if mol.iupac_name is None and mol.cirpy_done is False:
            print('getting IUPAC name from CIRPY...')
            mol.cirpy_to_iupac()
            if mol.iupac_name is None:
                mol.iupac_name
        # logP
        if mol.logP is None:
            print('getting logP from RDKit...')
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            if rdkitmol is None:
                mol.logP = 'not found'
                continue
            rdkitmol.Compute2DCoords()
            mol.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
        # XlogP
        if mol.XlogP is None:
            print('getting XlogP from PUBCHEM...')
            mol.XlogP = PUBCHEM_IO.hier_name_search(mol, 'XLogP')
            # if mol.iupac_name is not None:
            #     mol.XlogP = PUBCHEM_IO.get_logP_from_name(mol.iupac_name)
            # else:
            #     mol.XlogP = PUBCHEM_IO.get_logP_from_name(mol.name)
        # synthetic accessibility
        if mol.Synth_score is None:
            print('getting synthetic accessibility from RDKit...')
            rdkitmol = Chem.MolFromSmiles(mol.SMILES)
            if rdkitmol is None:
                mol.Synth_score = 'not found'
                continue
            rdkitmol.Compute2DCoords()
            mol.Synth_score = get_SynthA_score(rdkitmol)
        # complexity
        if mol.complexity is None:
            print('getting complexity from PUBCHEM...')
            mol.complexity = PUBCHEM_IO.hier_name_search(mol, 'complexity')
            # if mol.iupac_name is not None:
            #     mol.complexity = PUBCHEM_IO.get_complexity_from_name(mol.iupac_name)
            # else:
            #     mol.complexity = PUBCHEM_IO.get_complexity_from_name(mol.name)
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


def check_arbitrary_names(comp):
    """Check a database entry for an arbitrarily difficult name.

    i.e. swap 'H2O' with 'water'

    """
    list_of_changes = {'H2O': 'water',
                       'CO2': 'carbon dioxide',
                       'H+': 'hydron',
                       'NH3': 'ammonia'}
    if comp[0] in list(list_of_changes.keys()):
        new_comp = (list_of_changes[comp[0]], comp[1])
        return new_comp
    return comp


def charge_except():
    """Return a list of molecule SMILES that we allow even though they are
    charged.

    Hydroxide and Hydrogen ion so far.

    """
    li = ['[H+]', '[OH-]']
    return li


def check_mol_diam_per_pkl(filename):
    """This function is for debugging.

    """
    a = load_molecule(filename)
    print('name:', a.name)
    print('diam:', a.mid_diam)
    if input('do calc?') == 'T':
        res = calc_molecule_diameter(a.name, a.SMILES,
                                     out_dir=directory,
                                     vdwScale=0.8,
                                     boxMargin=4.0,
                                     spacing=0.6,
                                     N_conformers=50,
                                     MW_thresh=500)
        print('result:', res)
        if input('set to zero?') == 'T':
            a.min_diam = 0
            a.mid_diam = 0
            a.max_diam = 0
            # get avg values of all ratios of all conformers
            a.rat_1 = 0
            a.rat_2 = 0
            a.save_object(a.pkl)
        else:
            print('I did nothing.')


def change_all_pkl_suffixes(directory):
    """Change the suffixes of pkl file names in all molecules in directory.

    For Debugging

    """
    for i in yield_molecules(directory=directory):
        i.pkl = i.pkl.replace('.pkl', '.gpkl')
        i.pkl = i.pkl.replace('.bpkl', '.gpkl')
        new_rs_pkls = []
        try:
            for j in i.rs_pkls:
                new_rs_pkls.append(j.replace('.pkl', '.gpkl').replace('.bpkl', '.gpkl'))
        except AttributeError:
            pass
        i.rs_pkls = new_rs_pkls
        i.save_object(i.pkl)


if __name__ == "__main__":
    import rxn_syst
    import os
    import sys

    if (not len(sys.argv) == 7):
        print('Usage: molecule.py get_mol pop_mol mol_file update_KEGG update_lookup\n')
        print("""    get_mol: T for overwrite and collection of molecules from RS
        in current dir (F for skip)
        ---- this function is useful if you update the base attributes of the molecule class.
        """)
        print('    scratch: T for restart collection of molecules from RS.')
        print('    pop_mol: T to run population of molecule properties (does not overwrite).')
        print('    mol_file: file name of list of molecules, F if not specified.')
        print('    update_KEGG: T to update KEGG translator.')
        print('    update_lookup: T to update lookup file.')
        sys.exit()
    else:
        get_mol = sys.argv[1]
        scratch = sys.argv[2]
        pop_mol = sys.argv[3]
        mol_file = sys.argv[4]
        update_KEGG = sys.argv[5]
        update_lookup = sys.argv[6]

    if get_mol == 'T':
        print('extract all molecules from reaction syetms in current dir...')
        curr_dir = os.getcwd()
        done_file = curr_dir+'/done_RS.txt'
        print(curr_dir)
        get_all_molecules_from_rxn_systems(rxn_syst.yield_rxn_syst(curr_dir+'/'),
                                           from_scratch=scratch,
                                           done_file=done_file)

    if pop_mol == 'T':
        vdwScale = 0.8
        boxMargin = 4.0
        spacing = 0.6
        N_conformers = 50
        MW_thresh = 500
        print('settings:')
        print('    VDW scale:', vdwScale)
        print('    Box Margin:', boxMargin, 'Angstrom')
        print('    Grid spacing:', spacing, 'Angstrom')
        print('    No Conformers:', N_conformers)
        print('    MW threshold:', MW_thresh, 'g/mol')
        print('    Molecule file:', mol_file)
        inp = input('happy with these? (T/F)')
        if inp == 'F':
            sys.exit('change them in the source code')
        elif inp != 'T':
            sys.exit('I dont understand, T or F?')
        print('populate the properties attributes for all molecules in DB...')
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        if mol_file == 'F':
            mol_file = False
        populate_all_molecules(directory=directory,
                               vdwScale=vdwScale,
                               boxMargin=boxMargin,
                               spacing=spacing,
                               N_conformers=N_conformers,
                               MW_thresh=MW_thresh,
                               mol_file=mol_file)

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

    if update_lookup == 'T':
        lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
        with open(lookup_file, 'w') as f:
            f.write('SMILES___iupac___name___DB___DB_ID___KEGG_ID___InChiKey___CHEBI_ID___pkl\n')
        # iterate over all molecules in DB and if they have a KEGG ID then
        # write translation to CHEBI ID
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        for mol in yield_molecules(directory=directory):
            smiles = '-'
            iupac = '-'
            name = '-'
            DB = '-'
            DB_ID = '-'
            KEGG_ID = '-'
            CHEBI_ID = '-'
            IKEY = '-'
            pkl = '-'
            if mol.SMILES is not None:
                smiles = mol.SMILES
            if mol.iupac_name is not None:
                if type(mol.iupac_name) != list:
                    iupac = mol.iupac_name
            if mol.name is not None:
                name = mol.name
            if mol.DB is not None:
                DB = mol.DB
            if mol.DB_ID is not None:
                DB_ID = mol.DB_ID
            if 'KEGG' in mol.DB_list:
                KEGG_ID = mol.KEGG_ID
            if mol.chebiID is not None:
                CHEBI_ID = mol.chebiID
            if mol.InChIKey is not None:
                IKEY = mol.InChIKey
            if mol.pkl is not None:
                pkl = mol.pkl
            print(smiles, iupac, name, DB, DB_ID, KEGG_ID, CHEBI_ID, pkl)
            with open(lookup_file, 'a') as f:
                f.write(smiles+'___'+iupac+'___'+name+'___')
                f.write(DB+'___'+DB_ID+'___'+KEGG_ID+'___'+IKEY+'___')
                f.write(CHEBI_ID+'___'+pkl+'\n')
    # if do_plots = 'T':
    #     #######
    #     # molecule distributions
    #     #######
    #     # implement plotting of molecule properties without consideration of RS?
    #     # # categorize all molecules in mol output file
    #     # mol_categorical(mol_output_file=search_mol_output_file,
    #     #                                threshold=size_thresh,
    #     #                                output_dir=search_output_dir)
    #     # plot a distribution of all molecule complexity
    #     # mol_dist_complexity(output_dir=search_output_dir,
    #                         # generator=yield_rxn_syst(search_output_dir))
    sys.exit()

    from SABIO_IO import get_cmpd_information
    directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
    # change_all_pkl_suffixes(directory)
    for i in yield_molecules(directory=directory):
        if i.DB == 'SABIO':
            print(i.InChi)
            print(i.name, i.PubChemID)
            i.InChi = None
            get_cmpd_information(i)
            print(i.name, i.PubChemID)
            break

    i.rs_pkls
    i.name
    i.pkl
    i.SMILES
    print(i.min_diam)
    i.mid_diam
    i.max_diam
    i.rat_1
    mol_file = '/home/atarzia/psp/molecule_DBs/atarzia/ATRS_4782.gpkl'
    a = load_molecule(mol_file)
    a.__dict__
    a.SMILES = None
    a.save_object(a.pkl)
    # print(a.__dict__)
    check_mol_diam_per_pkl(mol_file)
