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
from os import getcwd
import glob
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from numpy import average
from molvs import standardize_smiles

import IO
from PUBCHEM_IO import hier_name_search_pcp
from rdkit_functions import calc_molecule_diameter


class Molecule:
    """
    Class that defines the molecules extracted from a database.

    """

    def __init__(self, name, role, DB, DB_ID, params):
        self.name = name
        if role == 'Substrate':
            self.role = 'reactant'
        else:
            self.role = role.lower()
        self.DB = DB
        self.DB_ID = DB_ID
        self.params = params
        self.InChi = None
        self.iupac_name = None
        self.mid_diam = None
        self.SMILES = None
        self.mol = None
        self.logP = None
        self.logS = None
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
        self.cirpy_done = False

    def molecule_db_prefix(self):
        return 'ATRS_'

    def determine_ID(self):
        """
        Determine an ID based on what is in my molecule_db_dir.

        """
        dir = self.params['molec_dir']
        pre = self.molecule_db_prefix()

        existing_pkls = glob.glob(dir+pre+'*.gpkl')
        existing_ids = [
            int(i.replace(dir+pre, '').replace('.gpkl', ''))
            for i in existing_pkls
        ]
        if len(existing_ids) > 0:
            max_id = max(existing_ids)
        else:
            max_id = 0
        return str(max_id+1)

    def get_pkl(self):
        pkl = self.params['molec_dir']+self.molecule_db_prefix()
        pkl += self.determine_ID()+'.gpkl'
        return pkl

    def save_molecule(self, filename):
        """
        Pickle molecule object to file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        # Overwrites any existing file.
        with gzip.GzipFile(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_molecule(self, filename, verbose=True):
        """
        unPickle molecule object from file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        if verbose:
            print('loading:', filename)
        with gzip.GzipFile(filename, 'rb') as input:
            self = pickle.load(input)
            return self

    def write_structure(self, filename):
        Chem.MolToMolFile(self.mol, filename)

    def PUBCHEM_last_shot(self):
        """
        Use PUBCHEM search for last chance at getting structure
        information.

        """
        # check for pubchem entry based on name
        print('search PUBCHEM...')
        # smiles = PUBCHEM_IO.get_SMILES_from_name(self.name)
        smiles_search = hier_name_search_pcp(self, 'CanonicalSMILES')
        if smiles_search is not None:
            if len(smiles_search) == 2:
                smiles = smiles_search[0]
                option = smiles_search[1]
            elif len(smiles_search) > 0:
                smiles = smiles_search
                option = 0
            else:
                smiles = None
        else:
            smiles = None
        if smiles is not None:
            self.SMILES = smiles
            print('>> SMILES:', self.SMILES)
            self.SMILES2MOL()
            self.InChiKey = hier_name_search_pcp(self, 'InChiKey',
                                                 option=option)
            print('>> IKEY:', self.InChiKey)
            self.iupac_name = hier_name_search_pcp(self, 'IUPACName',
                                                   option=option)
            print('>> iupac_name', self.iupac_name)

    def get_compound(self, dataset, search_mol=True):
        """Get compound SMILES from available identifiers.

        In this function we need to pick one CHEBI ID if available that
        describes the compound structure.

        """
        print('collecting information for:', self.name)
        # check if molecule exists in molecule database already
        old_pkl = None
        if search_mol:
            lookup_file = (
                '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
            )
            dataset = IO.read_molecule_lookup_file(
                lookup_file=lookup_file
            )
            old_pkl = search_molecule_by_ident(self, dataset)
        if old_pkl is not None:
            DB = self.DB
            # load old pkl
            print('collecting', self.name, 'from personal DB...')
            print('old pkl:', old_pkl)
            old_mol = Molecule.load_molecule(old_pkl, verbose=True)
            # copy old_mol DB object properties to self
            # only overwrite None or NaN
            for key in old_mol.__dict__:
                val = old_mol.__dict__[key]
                if key not in self.__dict__:
                    self.__dict__[key] = val
                elif self.__dict__[key] is None and val is not None:
                    self.__dict__[key] = val
            if DB not in self.DB_list:
                self.DB_list.append(DB)
            # change pkl of current molecule
            self.pkl = old_pkl
            return self
        elif self.SMILES is not None:
            print('have SMILES already...')
            self.SMILES2MOL()
        elif self.chebiID is not None:
            print('get compound using ChebiID...')
            if self.DB == 'KEGG':
                # name is set to KEGG C-ID at this point
                self.change_name = True
            from ercollect.CHEBI_IO import get_cmpd_information
            get_cmpd_information(self)
        if self.SMILES is None or self.chebiID is None:
            if self.DB == 'SABIO':
                print(
                    'No CHEBI structure, '
                    'try get_compound using SABIO DB...'
                )
                from ercollect.SABIO_IO import get_cmpd_information
                # set DB specific properties
                self.cID = self.DB_ID
                get_cmpd_information(self)
        if isinstance(self.chebiID, list):
            print(self.name, self.chebiID, self.SMILES, self.pkl)
            self.chebiID = ' '.join(self.chebiID)
        if isinstance(self.DB_ID, list):
            self.DB_ID = ' '.join(self.DB_ID)
            # input('why? - should only occur if no structure found')
        # if self.SMILES is None and self.mol is None:
        #     print('get compound using PUBCHEM as last shot...')
        #     self.PUBCHEM_last_shot()

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
        From PUBCHEM:
            - molecule complexity:
                https://pubchemdocs.ncbi.nlm.nih.gov/glossary
                $Complexity
                (0 to inf)
            - XlogP:
                https://pubchemdocs.ncbi.nlm.nih.gov/glossary$XLogP
                (smaller = more hydrophilic)
        """
        print('collect molecular properties using RDKit and PUBCHEM.')
        # logP and SA from RDKIT with SMILES:
        rdkitmol = Chem.MolFromSmiles(self.SMILES)
        rdkitmol.Compute2DCoords()
        self.logP = Descriptors.MolLogP(rdkitmol, includeHs=True)
        self.logS = get_logSw(rdkitmol)
        self.Synth_score = get_SynthA_score(rdkitmol)
        # XlogP and complexity from PUBCHEM with SMILES:
        # check if already calculated
        # set NAME based on available DB IDs
        if check:
            if self.XlogP is None or self.complexity is None:
                result = hier_name_search_pcp(molecule=self,
                                              property=['XLogP',
                                                        'complexity'])
                if result is not None:
                    self.XlogP, self.complexity = result
        else:
            result = hier_name_search_pcp(molecule=self,
                                          property=['XLogP',
                                                    'complexity'])
            if result is not None:
                self.XlogP, self.complexity = result

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

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(name={self.name}, role={self.role}, DB={self.DB}, '
            f'ID={self.DB_ID}, SMILEs={self.SMILES}, pkl={self.pkl})'
        )

    def __repr__(self):
        return str(self)


def get_logSw(mol):
    """Get water solubility using RDKit as described here:
    https://github.com/PatWalters/solubility.
    Using the newly paramterized function.

    """
    from ercollect.solubility.esol import ESOLCalculator
    esol_calculator = ESOLCalculator()
    logS = esol_calculator.calc_esol(mol)
    return logS


def get_SynthA_score(mol):
    """Get synthetic accesibility score from RDKIT contrib (SA_score).

    """
    from ercollect.SA_score import sa_scores
    s = sa_scores.calculateScore(mol)
    return s


def check_molecule_unique(molec, molecules):
    """Check if a molecule is unique compared to those in molecule DB.

    1 - check for same SMIlEs


    Returns the pkl file of the original molecule also.

    """
    if molec.SMILES is not None:
        for original_pkl in molecules:
            o_mol = Molecule.load_molecule(original_pkl, verbose=False)
            if o_mol.SMILES == molec.SMILES:
                unique = False
                return unique, original_pkl
    return True, None


def define_ident_list(molec):
    """Define ident search list.

    """
    ident_list = []
    if molec.SMILES is not None:
        ident_list.append(('SMILES', molec.SMILES))
    if molec.iupac_name is not None:
        ident_list.append(('iupac', molec.iupac_name))
    ident_list.append(('name', molec.name))
    ident_list.append(('name', molec.name.lower()))
    ident_list.append(('iupac', molec.name))
    try:
        if molec.KEGG_ID is not None:
            ident_list.append(('KEGG_ID', molec.KEGG_ID))
    except AttributeError:
        pass
    try:
        if molec.InChiKey is not None:
            ident_list.append(('InChiKey', molec.InChiKey))
    except AttributeError:
        pass
    try:
        if molec.CHEBI_ID is not None and isinstance(
            molec.CHEBI_ID, list
        ) is False:
            ident_list.append(('CHEBI_ID', molec.CHEBI_ID))
    except AttributeError:
        pass
    return ident_list


def search_molecule_by_ident(molec, dataset):
    """Search for molecule in molecule database using lookup file.

    1 - check for same SMILEs if SMILEs is not None
    2 - check for same IUPAC name if not None
    3 - check for same name
    4 - check for same DB and DB_ID -- removed on 05/12/18
    5 - check for same KEGG_ID
    6 - check for same InChiKey
    7 - check for same chebiID

    Returns the pkl file of the original molecule also.

    """
    ident_list = define_ident_list(molec)
    for search in ident_list:
        match = dataset[dataset[search[0]] == search[1]]
        if len(match) > 0:
            break
    if len(match) == 0:
        return None
    if len(match) > 1:
        print(match.T)
        print(match.pkl)
        import sys
        sys.exit('problem here with multiple hits in molecule DB')
    # the following code could be improved with the use of the above
    # code its on the todo list
    for idx, row in match.iterrows():
        if molec.SMILES is not None:
            if row['SMILES'] == molec.SMILES:
                print('>> found match with SMILES')
                return row['pkl']
        if molec.iupac_name is not None:
            if row['iupac'] == molec.iupac_name:
                print('>> found match with IUPAC name')
                return row['pkl']
        if row['name'] == molec.name:
            print('>> found match with name')
            return row['pkl']
        if row['name'] == molec.name.lower():
            print('>> found match with name')
            return row['pkl']
        if row['iupac'] == molec.name:
            print('>> found match with IUPAC name')
            return row['pkl']
        try:
            if row['KEGG_ID'] == molec.KEGG_ID:
                print('>> found match with KEGG ID')
                return row['pkl']
        except AttributeError:
            pass
        try:
            if row['InChiKey'] == molec.InChiKey:
                print('>> found match with InChiKey')
                return row['pkl']
        except AttributeError:
            pass
        try:
            if row['CHEBI_ID'] == molec.chebiID:
                print('>> found match with Chebi ID')
                return row['pkl']
        except AttributeError:
            pass
    return None


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


def iterate_rs_components(rs, molecule_dataset):
    """Iterate over all components in a reaction system and collect the
    component structures/properties.

    All changes to the reaction system should be done in-place.

    Arguments:
        rs (rxn_syst.reaction) - reaction system being tested
        molecule_dataset (Pandas DataFrame) - look up for known
        molecules

    """
    for m in rs.components:
        # # need to make sure that the role of this molecule matches
        # this RS
        # SET_role = m.role
        # translation only applies to molecules with KEGG IDs
        # which means we were able to collect all properties already.
        try:
            if m.translated is True:
                continue
        except AttributeError:
            m.translated = False
        m.get_compound(dataset=molecule_dataset,
                       search_mol=True)
        if m.SMILES is None:
            print('One SMILES not found in get_compound - skip.')
            rs.skip_rxn = True
            rs.skip_reason = 'one component has no SMILES'
            IO.fail_list_write(
                new_name=m.name,
                directory='/home/atarzia/psp/molecule_DBs/atarzia/',
                file_name='failures.txt')
            break
        else:
            # standardize SMILES
            print("smiles:", m.SMILES)
            try:
                m.SMILES = standardize_smiles(m.SMILES)
                print("standardized smiles:", m.SMILES)
            except ValueError:
                print('standardization failed - therefore assume')
                print('SMILES were invalid - skip')
                m.SMILES = None
                rs.skip_rxn = True
                rs.skip_reason = 'one component has invalid SMILES'
                IO.fail_list_write(
                    new_name=m.name,
                    directory=(
                        '/home/atarzia/psp/molecule_DBs/atarzia/'
                    ),
                    file_name='failures.txt'
                )
                # import sys
                # sys.exit()
                break
            # check for charge in SMILES
            # if check_charge_on_SMILES(m.SMILES):
            #     if charge_except(m.SMILES):
            #         # charged SMILES is in excepted cases
            #         pass
            #     else:
            #         # skip rxn
            #         print('One SMILES is charged - skip.')
            #         print('>>>>', m.name)
            #         print('>>>>', m.pkl)
            #         print('>>>>', m.SMILES)
            #         rs.skip_rxn = True
            #         rs.skip_reason = 'one component has charged
            # SMILES'
            #         fail_list_write(
            #             new_name=m.name,
            #             directory='/home/atarzia/psp/molecule_DBs/
            # atarzia/',
            #             file_name='failures.txt')
            #         break
            # check for wildcard in SMILES
            if '*' in m.SMILES:
                # skip rxn
                print('One SMILES contains wildcard - skip.')
                rs.skip_rxn = True
                rs.skip_reason = 'one component has wildcard SMILES'
                IO.fail_list_write(
                    new_name=m.name,
                    directory=(
                        '/home/atarzia/psp/molecule_DBs/atarzia/'
                    ),
                    file_name='failures.txt'
                )
                break
        m.get_properties()
        # add rxn syst pkl name
        # try:
        #     if rs.pkl not in m.rs_pkls:
        #         m.rs_pkls.append(rs.pkl)
        # except AttributeError:
        #     m.rs_pkls = []
        #     m.rs_pkls.append(rs.pkl)
        # # need to make sure that the role of this molecule matches
        #  this RS
        # m.role = SET_role
        # # save to molecule database
        # if os.path.isfile(m.pkl) is False:
        #     m.save_object(m.pkl)
        #     # update kegg translation and molecule look up file
        #     update_KEGG_translator()
        #     # update molecule look up file
    # once al components have been collected and skip_rxn is False
    # update the molecule DB and reread lookup_file
    if rs.skip_rxn is not True:
        print('--- updating molecule DB ---')
        done_file = getcwd()+'/done_RS.txt'
        # reload molecule data set
        lookup_file = (
            '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
        )
        molecule_dataset = IO.read_molecule_lookup_file(
            lookup_file=lookup_file
        )
        IO.update_molecule_DB(rxns=[rs], done_file=done_file,
                           dataset=molecule_dataset, from_scratch='T')
        # done_list_write(
        #     new_name=m.name,
        #     pkl=m.pkl,
        #     directory='/home/atarzia/psp/molecule_DBs/atarzia/')


def populate_all_molecules(directory, vdwScale, boxMargin, spacing,
                           N_conformers, MW_thresh, mol_file=False):
    """Populate all molecules in pickle files in directory.

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
        # XlogP + complexity
        if mol.XlogP is None or mol.complexity is None:
            print('getting XlogP + complexity from PUBCHEM...')
            result = hier_name_search_pcp(molecule=mol,
                                          property=['XLogP',
                                                    'complexity'])
            if result is not None:
                mol.XlogP, mol.complexity = result
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


def check_arbitrary_names(comp):
    """Check a database entry for an arbitrarily difficult name.

    i.e. swap 'H2O' with 'water'

    """
    list_of_changes = {'H2O': 'water',
                       'CO2': 'carbon dioxide',
                       'H+': 'hydron',
                       'NH3': 'ammonia',
                       'O2': 'dioxygen',
                       'H2O2': 'hydrogen peroxide',
                       'NO': 'nitric oxide',
                       'CN-': 'cyanide'}
    if type(comp) is tuple:
        if comp[0] in list(list_of_changes.keys()):
            new_comp = (list_of_changes[comp[0]], comp[1])
            return new_comp
        else:
            return comp
    elif type(comp) is str:
        if comp in list(list_of_changes.keys()):
            return list_of_changes[comp]
        else:
            return comp


def check_charge_on_SMILES(SMILES):
    """Return True if a SMILES is formally charged.

    This function is duplicated in PUBCHEM_IO.py.

    """
    mol = Chem.MolFromSmiles(SMILES)
    Chem.AddHs(mol)
    FC = Chem.GetFormalCharge(mol)
    if FC == 0:
        return False
    else:
        return True


def charge_except(SMILES):
    """
    Return True if SMILES matches one of the exempt charge species or
    if charge is balanced.

    Hydroxide and Hydrogen ion so far.

    """
    li = ['[H+]', '[OH-]']
    if SMILES in li:
        return True

    # certain SMILES strings that we allow
    # nitro groups
    li = ['[N+]([O-])=O', 'O=[N+]([O-])', '[N+](=O)[O-]']
    for l in li:
        # return True if SMILES contains target strings and SMILES
        # without target strings is uncharged
        if l in SMILES:
            if check_charge_on_SMILES(SMILES.replace(l, '')) is False:
                return True
    return False


def check_mol_diam_per_pkl(filename):
    """
    This function is for debugging.

    """
    a = Molecule.load_molecule(filename)
    print('name:', a.name)
    print('diam:', a.mid_diam)
    if input('do calc?') == 'T':
        res = calc_molecule_diameter(
            a.name,
            a.SMILES,
            out_dir=directory,
            vdwScale=0.8,
            boxMargin=4.0,
            spacing=0.6,
            N_conformers=50,
            MW_thresh=500
        )
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
