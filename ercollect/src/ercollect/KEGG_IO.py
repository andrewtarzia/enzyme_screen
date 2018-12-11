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
from sys import exit
from rdkit.Chem import AllChem as Chem
import requests
from ercollect import DB_functions
from ercollect import rxn_syst
from os import getcwd
from os.path import isfile
from ercollect.molecule import molecule, load_molecule, fail_list_read, \
                               fail_list_write, read_molecule_lookup_file, \
                               update_molecule_DB


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
        if isfile(output_dir+rs.pkl) is True and clean_system is False:
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
            iterate_rs_components_KEGG(rs, molecule_dataset=molecule_dataset)
        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
        count += 1


def get_components(equations_string):
    """Get reactants and products from a KEGG equation string.

    Returns bool (reversible), list of components.
    """
    if '<=>' in equations_string:
        # implies it is reversible
        reactants, products = equations_string.split("<=>")
        reversible = True
    else:
        reactants, products = equations_string.split("=")
        reversible = False
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
    return reversible, comp_list


def modify_MOLBlock(string):
    """Modify Mol block string from KEGG API to have a file source type.
    This means setting line 2 to be '     RDKit          2D'. Without this, no
    chiral information would be collected.

    """
    string = string.split('\n')
    string[1] = '     RDKit          2D'
    string = '\n'.join(string)
    return string


def convert_MOL_to_SMILES(string):
    """Convert MOL (in text) from KEGG website into rdkit Molecule object and
    SMILES.

    Returns RDKIT molecule, SMILES
    """
    string = modify_MOLBlock(string)
    mol = Chem.MolFromMolBlock(string)
    if mol is not None:
        smiles = Chem.MolToSmiles(mol)
        return mol, smiles
    else:
        # handle cases where conversion cannot occur
        # will add more cases as I discover them
        if ' R ' in string or ' * ' in string:
            # assume this means generic structure
            return 'generic'


def KEGGID_to_MOL(KEGG_ID):
    """Use KEGG API to collect MOL file using KEGG API.

    Examples: https://www.kegg.jp/kegg/rest/keggapi.html#get

    """
    # get compound information from KEGG API
    URL = 'http://rest.kegg.jp/get/'+KEGG_ID+'/mol'
    request = requests.post(URL)
    # KEGG API will give a 404 if MOL file cannot be collected
    if request.status_code == 404:
        print(KEGG_ID, 'has no mol file associated with it')
        return None
    elif request.status_code == 200:
        output = request.text
        print('convert', KEGG_ID, 'to RDKit MOL...')
        result = convert_MOL_to_SMILES(output)
        return result
    else:
        print("haven't come across this yet, figure out how to handle.")
        print('KEGG ID:', KEGG_ID, 'gives this error')
        exit('exitting....')


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
    fail_list = fail_list_read(
                    directory='/home/atarzia/psp/molecule_DBs/atarzia/',
                    file_name='failures.txt')
    # because of the formatting of KEGG text - this is trivial
    equations_string = request.text.split('EQUATION    ')[1].split('\n')[0].rstrip()
    rs.components = []
    rs.reversible, comp_list = get_components(equations_string)
    for comp in comp_list:
        if comp[0] in fail_list:
            rs.skip_rxn = True
            rs.skip_reason = 'KEGG ID could not be converted to MOL'
            print('one molecule in fail list - skipping...')
            print('>>> ', comp[0], 'failed')
            break
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
            result = KEGGID_to_MOL(KEGG_ID=comp[0])
            if result is None:
                print('KEGG ID could not be converted to MOL')
                print(' - skipping whole reaction.')
                print('>>>>>>', comp)
                rs.skip_rxn = True
                rs.skip_reason = 'KEGG ID could not be converted to MOL'
                # write KEGG ID to fail list
                fail_list_write(
                    new_name=comp[0],
                    directory='/home/atarzia/psp/molecule_DBs/atarzia/',
                    file_name='failures.txt')
                return rs
            elif result == 'generic':
                print('KEGG ID gave generic structure')
                print(' - skipping whole reaction.')
                rs.skip_rxn = True
                rs.skip_reason = 'KEGG ID gave generic structure'
                # write KEGG ID to fail list
                fail_list_write(
                    new_name=comp[0],
                    directory='/home/atarzia/psp/molecule_DBs/atarzia/',
                    file_name='failures.txt')
                return rs
            else:
                new_mol = molecule(comp[0], comp[1], 'KEGG', comp[0])
                # add new_mol to reaction system class
                new_mol.mol = result[0]
                new_mol.SMILES = result[1]
                new_mol.KEGG_ID = comp[0]
                new_mol.translated = False
                rs.components.append(new_mol)
    return rs


def iterate_rs_components_KEGG(rs, molecule_dataset):
    """Iterate over all components in a reaction system and collect the
    component structures/properties.

    > This is a version of the same function in molecule.py
    (iterate_rs_components) that was designed to handle
    KEGG collected structures implemented 12/12/18.

    All changes to the reaction system should be done in-place.

    Arguments:
        rs (rxn_syst.reaction) - reaction system being tested
        molecule_dataset (Pandas DataFrame) - look up for known molecules

    """
    for m in rs.components:
        # # need to make sure that the role of this molecule matches this RS
        # SET_role = m.role
        # translation only applies to molecules with KEGG IDs
        # which means we were able to collect all properties already.
        if m.translated is True:
            continue
        # all molecules should have a mol and SMILES attribute at this point
        # so remove m.get_compound() and other checks applied in the
        # molecule.py version
        # check for wildcard in SMILES
        if '*' in m.SMILES:
            # skip rxn
            print('One SMILES contains wildcard - skip.')
            rs.skip_rxn = True
            rs.skip_reason = 'one component has wildcard SMILES'
            fail_list_write(
                new_name=m.name,
                directory='/home/atarzia/psp/molecule_DBs/atarzia/',
                file_name='failures.txt')
            break
        m.get_properties()
    # once all components have been collected and skip_rxn is False
    # update the molecule DB and reread lookup_file
    if rs.skip_rxn is not True:
        print('--- updating molecule DB ---')
        done_file = getcwd()+'/done_RS.txt'
        # reload molecule data set
        lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
        molecule_dataset = read_molecule_lookup_file(lookup_file=lookup_file)
        update_molecule_DB(rxns=[rs], done_file=done_file,
                           dataset=molecule_dataset, from_scratch='T')
