#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of BKMS DB.

Using downloaded flat files. Offline.

Author: Andrew Tarzia

Date Created: 30 Aug 2018
"""

import DB_functions
import CHEBI_IO
from rdkit.Chem import AllChem as Chem
import pandas as pd


def init_BKMS():
    """Output some hardcoded things for BKMS usage.

    """
    bkms_dir = '/home/atarzia/psp/BKMS_react/'

    # NaNs in Remark column replaced with 999 - all others remain
    bkms_data = pd.read_table(bkms_dir+'Reactions_BKMS.tab', delimiter='\t',
                              names=['ID', 'EC', 'Rec. Name', 'rxn',
                                     'RID_Brenda', 'RID_KEGG', 'RID_metacyc',
                                     'RID_SABIO', 'Brenda_pathway_name',
                                     'KEGG_pathway_ID', 'KEGG_pathway_name',
                                     'metacyc_pathway_id',
                                     'metacyc_pathway_name',
                                     'stoich_check', 'missing substrate',
                                     'missing product', 'KEGG_comments',
                                     'metacyc_comments', 'remark']).fillna(value={'remark': 999})

    print("BKMS Stats:")
    print("The table contains actual data of BRENDA (release 2018.2, only",
          "reactions with naturally occuring substrates), MetaCyc (version",
          "21.5), SABIO-RK (02/05/2018) and KEGG data, downloaded 23/04/2012.",
          "(Downloading more recent KEGG data cannot be offered because a",
          "KEGG license agreement would be necessary.)", sep='\n')
    print('----------------')
    print("independant EC No.:", len(list(set(bkms_data['EC']))))
    print("independant Rxns:", len(bkms_data['EC']))
    print("independant non-generic Rxns:", len([i for i in list(bkms_data['remark'][bkms_data['remark'] != 999]) if 'generic' not in i]))
    return bkms_dir, bkms_data


def check_for_radicals(mol_list):
    """Check if list of molecules contains a radical.

    i.e. lignin implies polymeric species.

    """

    to_remove = []
    for name in mol_list:
        if 'radical' in name:
            # new_name = name.replace('radical', '').lstrip().rstrip()
            # if new_name in mol_list:
            print('removing', name, 'as it is a radical')
            to_remove.append(name)

    return to_remove


def skip_names(mol_list):
    """Check if a list of molecules contains nomenclature that suggests it can
    be skipped.

    i.e. lignin implies polymeric species.

    """

    # lignin - a polymeric species that can form under oxidative polymerisation
    #          many forms exist. HRP produces some.

    # add white space to avoid getting the word in the middle of other terms
    skippable_strings = [' lignin']

    for name in mol_list:
        for i in skippable_strings:
            if i in name:
                print('skipping this molecule because it contains:', i)
                return True
    return False


def get_SMILES_for_molecule_list(mol_list, DBs='any'):
    """Convert list of molecule names to Canonical SMILEs by searching DBs.

    Keywords:
        mol_list (list) - list of molecule names
        DBs (str) - DB to use, defaults to a hierachical search through all
            available

    DBs available (online.offline) (in 'any' order):
        - CHEBI (offline) - not yet
        - CHEMBL (online) - not yet
        - SABIO (online) - not yet
        - KEGG (online) - not yet

    Returns:
        mol_dict (dict) - {name: (SMILEs, DB, DB_ID, iupac_name)}
            DB_ID is the molecule ID within the given DB.

    """
    mol_dict = {}

    # get properties of DBs
    DB_prop = DB_functions.get_DB_prop(DB=DBs)

    db_dir = DB_prop['CHEBI'][0]

    # chebi code - temp
    compounds_file = db_dir+'compounds.tsv'
    names_file = db_dir+'names.tsv'
    structures_file = db_dir+'structures.csv'

    for mol in mol_list:
        print(mol)
        # search for name in compound file
        res = CHEBI_IO.search_for_compound_by_name(compounds_file, mol)
        if res is None:
            # search for formula in names file
            res = CHEBI_IO.search_for_name_by_name(names_file, mol)
            if res is None:
                print('no match in DB')
                continue
            else:
                ID, name = res
                parent_id = None
        else:
            ID, parent_id, name, star = res

        # make sure is parent compound
        if parent_id != 'null':
            res = CHEBI_IO.convert_nameID_to_parent(compounds_file, nameID=ID)
            if res is None:
                print("this should not happen - error with cross reference")
                print('check this!')
                import sys
                sys.exit()
            ID, parent_id, name, star = res

        # get structure using CHEBI ID
        # structures.csv - read in, get COMPOUND ID match then extract the
        # structure
        # read into RDKIT straight up OR save the SMILE or save a pickle
        # mol_dict should refeence the SMILES or the file its saved to.
        structure, s_type = CHEBI_IO.get_structure(structures_file, ID)
        if structure is not None:
            # is structure a MolBlock or Smiles
            if s_type == 'mol':
                # convert structure to SMILEs
                rdkitmol = Chem.MolFromMolBlock(structure)
                rdkitmol.Compute2DCoords()
                smile = Chem.MolToSmiles(rdkitmol)
            elif s_type == 'SMILES':
                smile = structure
            elif s_type == 'InChIKey':
                rdkitmol = Chem.MolFromInchi(structure)
                rdkitmol.Compute2DCoords()
                smile = Chem.MolToSmiles(rdkitmol)
        else:
            smile = '-'
            print('molecule does not have recorded structure in DB')
            print('need to search another DB')

        DB = DBs  # temp
        DB_ID = ID
        iupac_name = name
        mol_dict[mol] = (smile, DB, DB_ID, iupac_name)

    return mol_dict
