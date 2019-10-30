#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of BKMS DB.

Using downloaded flat files. Offline.

Author: Andrew Tarzia

Date Created: 30 Aug 2018
"""
from ercollect import DB_functions
from ercollect import CHEBI_IO
import pandas as pd
from ercollect import rxn_syst
import os
from ercollect.molecule import (
    molecule, iterate_rs_components, check_arbitrary_names
)
from ercollect.molecule import fail_list_read, fail_list_write
from ercollect import PUBCHEM_IO


def init_BKMS(bkms_dir, verbose=False):
    """Output some hardcoded things for BKMS usage.

    """

    # NaNs in Remark column replaced with 999 - all others remain
    bkms_data = pd.read_table(
        bkms_dir+'Reactions_BKMS.tab',
        delimiter='\t',
        names=[
            'ID', 'EC', 'Rec. Name', 'rxn',
            'RID_Brenda', 'RID_KEGG', 'RID_metacyc',
            'RID_SABIO', 'Brenda_pathway_name',
            'KEGG_pathway_ID', 'KEGG_pathway_name',
            'metacyc_pathway_id',
            'metacyc_pathway_name',
            'stoich_check', 'missing substrate',
            'missing product', 'KEGG_comments',
            'metacyc_comments', 'remark'
        ]
    ).fillna(value={'remark': 999})

    if verbose:
        print("BKMS Stats:")
        print(
            "The table contains actual data of BRENDA (release 2018.2,"
            " only"
            "reactions with naturally occuring substrates), MetaCyc "
            "(version",
            "21.5), SABIO-RK (02/05/2018) and KEGG data, downloaded "
            "23/04/2012.",
            "(Downloading more recent KEGG data cannot be offered "
            "because a",
            "KEGG license agreement would be necessary.)",
            sep='\n'
        )
        print('----------------')
        print("independant EC No.:", len(list(set(bkms_data['EC']))))
        print("independant Rxns:", len(bkms_data['EC']))
        print("independant non-generic Rxns:", len([
            i for i in list(
                bkms_data['remark'][bkms_data['remark'] != 999]
            ) if 'generic' not in i
        ]))
    return bkms_data


def check_for_radicals(mol_list):
    """Check if list of molecules contains a radical.

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
    """
    Check if a list of molecules contains nomenclature that suggests
    it can be skipped.

    i.e. lignin implies polymeric species.

    """

    # lignin - a polymeric species that can form under oxidative
    # polymerisation
    #          many forms exist. HRP produces some.

    # add white space to avoid getting the word in the middle of other
    # terms
    skippable_strings = [' lignin']

    for name in mol_list:
        for i in skippable_strings:
            if i in name:
                print('skipping this molecule because it contains:', i)
                return True
    return False


def get_rxn_systems(EC, output_dir, molecule_dataset,
                    clean_system=False, verbose=False):
    """
    Get reaction systems from BKMS entries in one EC and output to
    Pickle.

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
        # # check if SABIO ID, BRENDA ID or KEGG ID for reaction has
        # already been
        # # collected
        # BRENDA_ID = str(row['RID_Brenda'])
        # KEGG_ID = str(row['RID_KEGG'])
        # SABIO_ID = str(row['RID_SABIO'])
        # # if BRENDA_ID != 'nan':
        # #     # check for BRENDA files of same ID
        # #     # skip
        # if KEGG_ID != 'nan':
        #     # check for KEGG files of same ID
        #     # skip
        #     kegg_files = glob.glob(output_dir+'*KEGG*')
        #     kegg_files = [i.replace(output_dir+'sRS-', '') for i in
        #     kegg_files]
        #     kegg_files = [i.replace(EC.replace('.' , '_'), '') for
        #     i in
        #     kegg_files]
        #     kegg_files = [i.replace('-KEGG-', '') for i in
        #     kegg_files]
        #     RIDs = [i.replace('.pkl', '') for i in kegg_files]
        #     if KEGG_ID in RIDs:
        #         continue
        # if SABIO_ID != 'nan':
        #     # check for SABIO files of same ID
        #     # skip
        #     sabio_files = glob.glob(output_dir+'*SABIO*')
        #     sabio_files = [i.replace(output_dir+'sRS-', '')
        #     for i in sabio_files]
        #     sabio_files = [i.replace(EC.replace('.' , '_'), '')
        #     for i in sabio_files]
        #     sabio_files = [i.replace('-SABIO-', '')
        #     for i in sabio_files]
        #     RIDs = [i.replace('.pkl', '') for i in sabio_files]
        #     if SABIO_ID in RIDs:
        #         continue

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
            if clean_system is False:
                count += 1
                continue
        if verbose:
            print('DB: BKMS - EC:', EC, '-',
                  'DB ID:', bkms_id, '-', count, 'of', len(EC_data))
        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.DB_ID, row)
        if rs.skip_rxn is False:
            # append compound information
            iterate_rs_components(
                rs,
                molecule_dataset=molecule_dataset
            )
        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
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
        reactants, products = row['rxn'].split(" <=> ")
        print("reaction is reversible")
    else:
        reactants, products = row['rxn'].split(" = ")
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
        rs.skip_reason = 'a component is in skip_names'
        return rs
    # check if the reactants or products contain the term radical
    # here we will assume that structurally the non radical can
    # represent the radical - and the radical component is ignored.
    # don't actually do this for now 06/12/18 #######################
    # rad_to_remove = check_for_radicals(new_reactants+new_products)
    # if len(rad_to_remove) > 0:
    #     new_reactants = [
    #       i for i in new_reactants if i not in rad_to_remove
    #   ]
    #     new_products = [
    #       i for i in new_products if i not in rad_to_remove
    #   ]

    # define component list
    comp_list = []
    for i in new_reactants:
        comp_list.append((i, 'reactant'))
    for i in new_products:
        comp_list.append((i, 'product'))

    rs.components = []
    fail_list = fail_list_read(
        directory='/home/atarzia/psp/molecule_DBs/atarzia/',
        file_name='failures.txt'
    )
    for comp in comp_list:
        # check if component name should be changed to a common name
        comp = check_arbitrary_names(comp)
        if comp[0] in fail_list:
            rs.skip_rxn = True
            rs.skip_reason = 'one component failed resolution'
            print('one molecule in fail list - skipping...')
            print('>>> ', comp[0], 'failed')
            break
        chebiID = CHEBI_IO.get_chebiID(comp[0])
        new_mol = molecule(comp[0], comp[1], 'BKMS', chebiID)
        if chebiID is None:
            result = PUBCHEM_IO.pubchem_synonym(new_mol)
            if result is not None:
                chebiID = result
                new_mol.DB_ID = chebiID
                new_mol.chebiID = chebiID
        if chebiID is None:
            new_mol, result = PUBCHEM_IO.pubchem_check_smiles(new_mol)
            if result is None:
                rs.skip_rxn = True
                rs.skip_reason = 'one component failed resolution'
                print('all failed - add to fail list + skipping...')
                fail_list_write(
                    new_name=comp[0],
                    directory=(
                        '/home/atarzia/psp/molecule_DBs/atarzia/'
                    ),
                    file_name='failures.txt')
                break
        # add new_mol to reaction system class
        else:
            new_mol.chebiID = chebiID
        rs.components.append(new_mol)

    return rs
