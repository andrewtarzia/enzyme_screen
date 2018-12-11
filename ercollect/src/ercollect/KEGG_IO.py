#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of KEGG DB.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import json
import pickle
import gzip
import requests
from ercollect import DB_functions
from ercollect import rxn_syst
import os
from ercollect.molecule import molecule, iterate_rs_components, load_molecule


def check_translator(ID):
    """Check for KEGG ID of molecule in KEGG translation file that links to
    molecules already collected from online.

    """
    translator = '/home/atarzia/psp/molecule_DBs/KEGG/translator.txt'
    # read translator
    with gzip.GzipFile(translator, 'rb') as output:
        translation = pickle.load(output)
    try:
        return translation[ID]
    except KeyError:
        return None


def get_EC_rxns_from_JSON(JSON_DB, EC):
    """Get reactions associated with EC number in KEGG.

    """
    try:
        rxn_list = JSON_DB[EC]
        return rxn_list
    except KeyError:
        return None


def get_rxn_systems(EC, output_dir, molecule_dataset,
                    clean_system=False, verbose=False):
    """Get reaction systems from KEGG entries in one EC and output to Pickle.

    """
    DB_prop = DB_functions.get_DB_prop('KEGG')
    # read in JSON file of whole DB
    rxn_DB_file = DB_prop[0]+DB_prop[1]['JSON_EC_file']
    with open(rxn_DB_file, 'r') as data_file:
        rxn_DB = json.load(data_file)

    # get EC specific entries
    EC_rxns = get_EC_rxns_from_JSON(rxn_DB, EC)
    if EC_rxns is None:
        return None
    # iterate over reactions
    count = 0
    for rxn in EC_rxns:
        # get KEGG rxn id
        string = rxn['name']
        K_Rid = string.split(' ')[0].rstrip()
        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'KEGG', K_Rid)
        if os.path.isfile(output_dir+rs.pkl) is True and clean_system is False:
            count += 1
            continue
        if verbose:
            print('======================================================')
            print('DB: KEGG - EC:', EC, '-',
                  'DB ID:', K_Rid, '-', count, 'of', len(EC_rxns))
        # there are no KEGG specific properties (for now)
        # could append pathways and orthologies
        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.DB_ID)
        if rs.skip_rxn is False:
            # append compound information
            iterate_rs_components(rs, molecule_dataset=molecule_dataset)
        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
        count += 1


def get_rxn_system(rs, ID):
    """Get reaction system from KEGG reaction ID (rID).

    Use KEGG API - online.

    Keywords:
        rs (class) - reaction system object
        ID (str) - DB reaction ID

    """
    # get Reaction information from KEGG API
    URL = 'http://rest.kegg.jp/get/reaction:'+ID
    request = requests.post(URL)
    if request.text != '':
        request.raise_for_status()
    else:
        rs.skip_rxn = True
        rs.skip_reason = 'No result for KEGG URL search - likely outdated'
        return rs
    # collate request output
    rs.components = []
    # because of the formatting of KEGG text - this is trivial
    equations_string = request.text.split('EQUATION    ')[1].split('\n')[0].rstrip()
    if '<=>' in equations_string:
        # implies it is reversible
        reactants, products = equations_string.split("<=>")
        print("reaction is reversible")
    else:
        reactants, products = equations_string.split("=")

    reactants = [i.lstrip().rstrip() for i in reactants.split("+")]
    products = [i.lstrip().rstrip() for i in products.split("+")]
    # collect KEGG Compound/Glycan ID
    # check if reactant or product are compound or glycan
    # remove stoichiometry for now
    comp_list = []
    for r in reactants:
        if 'G' in r:
            # is glycan
            KID = 'G'+r.split('G')[1].rstrip()
        elif 'C' in r:
            # is compound
            KID = 'C'+r.split('C')[1].rstrip()
        elif 'D' in r:
            # is drug
            KID = 'D'+r.split('D')[1].rstrip()
        comp_list.append((KID, 'reactant'))
    for r in products:
        if 'G' in r:
            # is glycan
            KID = 'G'+r.split('G')[1].rstrip()
        elif 'C' in r:
            # is compound
            KID = 'C'+r.split('C')[1].rstrip()
        elif 'D' in r:
            # is drug
            KID = 'D'+r.split('D')[1].rstrip()
        comp_list.append((KID, 'product'))

    for comp in comp_list:
        # check for CHEBI ID in KEGG translations file
        translated = check_translator(comp[0])
        if translated is not None:
            pkl = translated
            print('collecting KEGG molecule using translator:', comp[0])
            new_mol = load_molecule(pkl, verbose=True)
            new_mol.KEGG_ID = comp[0]
            new_mol.translated = True
            # need to make sure tranlated molecule has the correct role
            new_mol.role = comp[1]
            rs.components.append(new_mol)
        else:
            chebiIDs = KEGGID_to_CHEBIID(KEGG_ID=comp[0])
            if chebiIDs is None:
                print('CHEBI ID not available - skipping whole reaction.')
                rs.skip_rxn = True
                rs.skip_reason = 'CHEBI ID not available for one component'
                return rs
            elif chebiIDs is not None:
                new_mol = molecule(comp[0], comp[1], 'KEGG', chebiIDs)
                # add new_mol to reaction system class
                new_mol.KEGG_ID = comp[0]
                new_mol.translated = False
                new_mol.chebiID = chebiIDs
                rs.components.append(new_mol)
            else:
                new_mol = molecule(comp[0], comp[1], 'KEGG', comp[0])
                # add new_mol to reaction system class
                new_mol.KEGG_ID = comp[0]
                new_mol.translated = False
                new_mol.chebiID = None
                rs.components.append(new_mol)
    return rs
