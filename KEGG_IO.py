#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of SABIO DB.

Modified code from:
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/searchPython.gsp

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import json
import requests
import DB_functions
import rxn_syst
import os
import molecule


def get_EC_rxns_from_JSON(JSON_DB, EC):
    """Get reactions associated with EC number in KEGG.

    """
    EC_heir = [int(i) for i in EC.split('.')]
    for i in JSON_DB['children'][EC_heir[0]-1]['children']:
        for j in i['children']:
            for k in j['children']:
                EC_sub = k['name']
                if EC == EC_sub:
                    EC_rxn = k['children']
                    return EC_rxn
    return None


def get_rxn_systems(EC, output_dir, clean_system=False):
    """Get reaction systems from KEGG entries in one EC and output to Pickle.

    """
    DB_prop = DB_functions.get_DB_prop('KEGG')
    # read in JSON file of whole DB
    rxn_DB_file = DB_prop[0]+DB_prop[1]['JSON_file']
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
        print('DB: KEGG - EC:', EC, '-',
              'DB ID:', K_Rid, '-', count, 'of', len(EC_rxns))

        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'KEGG', K_Rid)
        if os.path.isfile(output_dir+rs.pkl) is True and clean_system is False:
            print('-----------------------------------')
            count += 1
            continue

        # there are no KEGG specific properties (for now)
        # could append pathways and orthologies

        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.DB_ID)
        if rs.skip_rxn is False:
            # append compound information - again DB specific
            for m in rs.components:
                m.get_compound()

        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
        print('-----------------------------------')
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
    request.raise_for_status()
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
        comp_list.append((KID, 'reactant'))
    for r in products:
        if 'G' in r:
            # is glycan
            KID = 'G'+r.split('G')[1].rstrip()
        elif 'C' in r:
            # is compound
            KID = 'C'+r.split('C')[1].rstrip()
        comp_list.append((KID, 'product'))

    for comp in comp_list:
        # get compound information from KEGG API
        # just convert to CHEBI ID and use CHEBI functions
        if 'C' in comp[0]:
            URL = 'http://rest.kegg.jp/conv/chebi/compound:'+comp[0]
        elif 'G' in comp[0]:
            URL = 'http://rest.kegg.jp/conv/chebi/glycan:'+comp[0]
        request = requests.post(URL)
        request.raise_for_status()
        # get CHEBI ID
        # because of the formatting of KEGG text - this is trivial
        if 'chebi' in request.text:
            chebiID = request.text.split('chebi:')[1].split('\n')[0].rstrip()
        else:
            print('CHEBI ID not available - skipping whole reaction.')
            rs.skip_rxn = True
            return rs

        new_mol = molecule.molecule(comp[0], comp[1], 'KEGG', chebiID)
        # add new_mol to reaction system class
        new_mol.KEGG_ID = comp[0]
        rs.components.append(new_mol)
    return rs
