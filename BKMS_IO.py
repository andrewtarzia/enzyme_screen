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
import rxn_syst
import os
import molecule


def init_BKMS(bkms_dir):
    """Output some hardcoded things for BKMS usage.

    """

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
    return bkms_data


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


def get_chebiID_from_BKMS(mol_name):
    """Convert molecule name to chebiID using CHEBI DB files.

    Offline.

    Keywords:
        mol_name (str) - molecule name

    Returns:
        ID (str) - chebiID
    """

    DB_prop = DB_functions.get_DB_prop('CHEBI')
    compounds_file = DB_prop[0]+DB_prop[1]['cmpds_file']
    names_file = DB_prop[0]+DB_prop[1]['names_file']
    structures_file = DB_prop[0]+DB_prop[1]['strct_file']

    # search for name in compound file
    res = CHEBI_IO.search_for_compound_by_name(compounds_file, mol_name)
    if res is None:
        # search for formula in names file
        res = CHEBI_IO.search_for_name_by_name(names_file, mol_name)
        if res is None:
            print('no match in DB')
            return None
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

    return ID


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


def get_rxn_systems(EC, output_dir):
    """Get reaction systems from BKMS entries in one EC and output to Pickle.

    """
    DB_prop = DB_functions.get_DB_prop('BKMS')

    # read in BKMS data
    bkms_data = init_BKMS(DB_prop[0])

    # get EC specific entries
    EC_data = bkms_data[bkms_data['EC'] == EC]

    # iterate over reactions
    count = 0
    for idx, row in EC_data.iterrows():
        # get BKMS ID
        bkms_id = row['ID']
        print('DB: BKMS - EC:', EC, '-',
              'DB ID:', bkms_id, '-', count, 'of', len(EC_data))
        # ignore those with 'generic' in remark
        if 'generic' in row['remark']:
            count += 1
            continue
        # ignore those with missing substrates
        if pd.isna(row['missing substrate']) is False:
            count += 1
            continue

        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'BKMS', bkms_id)
        if os.path.isfile(output_dir+rs.pkl) is True:
            print('-----------------------------------')
            count += 1
            continue

        # there are no KEGG specific properties (for now)
        # could append pathways and orthologies

        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.DB_ID, row)
        if rs.skip_rxn is False:
            # append compound information - again DB specific
            for m in rs.components:
                m.get_compound()

        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
        print('-----------------------------------')
        count += 1


def get_rxn_system(rs, ID, row):
    """Get reaction system from BKMS ID.

    Offline. Use PANDAS DataFrame row.

    Keywords:
        rs (class) - reaction system object
        ID (str) - DB reaction ID
        row (Pandas Series) - row associated with BKMS ID.

    """
    # define reactants and products
    if '<=>' in row['rxn']:
        # implies it is reversible
        reactants, products = row['rxn'].split("<=>")
        print("reaction is reversible")
    else:
        reactants, products = row['rxn'].split("=")
    reactants = reactants.split(" + ")
    products = products.split(" + ")
    # remove white space on the left and right ends
    reactants = [i.lstrip().rstrip() for i in reactants]
    products = [i.lstrip().rstrip() for i in products]
    # remove stoichiometry
    new_reactants = []
    for r in reactants:
        # have stoich?
        if r.split(' ')[0].isnumeric() is True:
            # yes
            new_reactants.append(' '.join(r.split(' ')[1:]))
        else:
            # no
            new_reactants.append(r)
    new_products = []
    for r in products:
        # have stoich?
        if r.split(' ')[0].isnumeric() is True:
            # yes
            new_products.append(' '.join(r.split(' ')[1:]))
        else:
            # no
            new_products.append(r)

    # check if any of the reactants or products should be skipped
    # because of their names.
    if skip_names(new_reactants+new_products) is True:
        # skip whole reaction if one component has skipped name
        rs.skip_rxn = True

    # check if the reactants or products contain the term radical
    # here we will assume that structurally the non radical can
    # represent the radical - and the radical component is ignored.
    rad_to_remove = check_for_radicals(new_reactants+new_products)
    if len(rad_to_remove) > 0:
        new_reactants = [i for i in new_reactants if i not in rad_to_remove]
        new_products = [i for i in new_products if i not in rad_to_remove]

    # deinfe component list
    comp_list = []
    for i in new_reactants:
        comp_list.append((i, 'reactant'))
    for i in new_products:
        comp_list.append((i, 'product'))

    rs.components = []
    for comp in comp_list:
        chebiID = get_chebiID_from_BKMS(comp[0])

        if chebiID is None:
            rs.skip_rxn = True
            continue
        new_mol = molecule.molecule(comp[0], comp[1], 'BKMS', chebiID)
        # add new_mol to reaction system class
        rs.components.append(new_mol)

    return rs
