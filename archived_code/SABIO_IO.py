#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of SABIO DB.

Modified code from:
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/searchPython.gsp

Author: Andrew Tarzia

Date Created: 15 Aug 2018

"""
import requests
from rdkit.Chem import AllChem as Chem
import os
from ercollect import rxn_syst
from ercollect.molecule import (
    molecule,
    iterate_rs_components,
    check_arbitrary_names,
    fail_list_read,
    load_molecule
)
from ercollect import CHEBI_IO
from ercollect.KEGG_IO import check_translator, KEGGID_to_CHEBIID
from molvs import standardize_smiles


def get_cmpd_information(molec):
    """Get information from SABIO Database of a compound with ID cID.

    """
    QUERY_URL = (
        'http://sabiork.h-its.org/sabioRestWebServices/'
        'searchCompoundDetails'
    )

    # input: SabioCompoundID
    # valid output fields: "fields[]":["Name","ChebiID",
    #                           "PubChemID","InChI",
    #                        "SabioCompoundID","KeggCompoundID"]
    params = {"SabioCompoundID": molec.cID,
              "fields[]": ["Name", "ChebiID", "PubChemID", "InChI"]}
    if molec.InChi is None:
        request = requests.post(QUERY_URL, params=params)
        request.raise_for_status()
        if request.text == 'No results found for query':
            molec.mol = None
        else:
            # results
            txt = request.text.split('\n')[1].split('\t')
            _, _, _, molec.InChi = txt
    if molec.InChi != 'null':
        print('collect SMILES from SABIO InChi')
        molec.mol = get_rdkit_mol_from_InChi(molec.InChi)
        if molec.mol is not None:
            smiles = Chem.MolToSmiles(Chem.RemoveHs(molec.mol))
            molec.SMILES = smiles
            try:
                molec.SMILES = standardize_smiles(molec.SMILES)
            except ValueError:
                print('standardization failed - therefore assume')
                print('SMILES were invalid - skip')
                molec.SMILES = None
                molec.mol = None
                # import sys
                # sys.exit()
        else:
            molec.SMILES = None
    else:
        molec.mol = None
        molec.SMILES = None


def get_rdkit_mol_from_InChi(InChi, AddHs=True):
    """Converts InChi into RDKIT molecule object.

    Arguments:
        InChi (str) - InChi code for molecule
        AddHs (bool) - default True

    Returns:
        mol (RDKIT molecule object)
    """
    try:
        mol = Chem.inchi.MolFromInchi(InChi, treatWarningAsError=True)
    except Chem.inchi.InchiReadWriteError:
        print(
            '>> RDKIT warning when converting InChi (SABIO) '
            'to molecule.'
        )
        return None
    if AddHs is True:
        mol = Chem.AddHs(mol)
    return mol


def get_entries_per_EC(EC):
    """Collect SABIO entry IDs associated with an EC no.

    Arguments:
        EC (str) - format: X.; -- X.X.X.X

    Returns:
        entries (list) - list of entries

    """
    ENTRYID_QUERY_URL = (
        'http://sabiork.h-its.org/sabioRestWebServices/'
        'searchKineticLaws/entryIDs'
    )

    # ask SABIO-RK for all EntryIDs matching a query
    query_dict = {"ECNumber": EC}
    query_string = ' AND '.join([
        '%s:%s' % (k, query_dict[k]) for k in query_dict
    ])
    query = {'format': 'txt', 'q': query_string}
    # make GET request
    request = requests.get(ENTRYID_QUERY_URL, params=query)
    request.raise_for_status()  # raise if 404 error
    if request.text == 'No results found for query':
        print(EC, ':', request.text)
        return []
    # each entry is reported on a new line
    entries = [int(x) for x in request.text.strip().split('\n')]
    print(EC, ': %d matching entries found.' % len(entries))

    return entries


def get_rxnID_from_eID(eID):
    """Collect SABIO reaction ID and properties from entry ID.

    """
    PARAM_QUERY_URL = (
        'http://sabiork.h-its.org/entry/exportToExcelCustomizable'
    )
    # encode next request, for parameter data given entry IDs
    data_field = {'entryIDs[]': eID}
    query = {'format': 'tsv',
             'fields[]': ['EntryID',
                          'Organism',
                          'ECNumber',
                          'SabioReactionID',
                          'UniprotID',
                          'ReactionEquation',
                          'EnzymeType']}
    # make POST request
    request = requests.post(
        PARAM_QUERY_URL,
        params=query,
        data=data_field
    )
    request.raise_for_status()
    # print(request.text)
    if request.text == 'No results found for query':
        print(request.text)
        return None, None, None
    # results
    RES = request.text.split('\n')[1].split('\t')
    _, organism, _, rxn_id, UniprotID, RE, enzymetype = RES
    return organism, rxn_id, UniprotID, RE, enzymetype


def get_rxn_systems(EC, output_dir, molecule_dataset,
                    clean_system=False, verbose=False):
    """Get reaction systems from SABIO entries in one EC and output to
    Pickle.

    """
    # get all SABIO entries
    entries = get_entries_per_EC(EC)
    # iterate over entries
    count = 0
    for eID in entries:
        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'SABIO', eID)
        if os.path.exists(output_dir+rs.pkl) is True:
            if clean_system is False:
                count += 1
                continue
        if verbose:
            print('===============================================')
            print('DB: SABIO - EC:', EC, '-',
                  'DB ID:', eID, '-', count, 'of', len(entries))
        # get reaction ID
        # SABIO specific properties
        RES = get_rxnID_from_eID(eID)
        rs.organism, rs.rID, rs.UniprotID, _, rs.etype = RES
        # etype = wildtype OR mutant
        # skip all mutants
        if 'wildtype' not in rs.etype:
            rs.skip_rxn = True
            rs.skip_reason = 'SABIO E-ID is for mutant'
            print('SABIO E-ID is mutant - skipping...')
        if rs.rID is None:
            rs.skip_rxn = True
            rs.skip_reason = 'SABIO R-ID not found'
            print('SABIO R-ID not found - skipping...')
        # get reaction system using DB specific function
        if rs.skip_rxn is False:
            rs = get_rxn_system(rs, rs.rID)
        if rs.skip_rxn is False:
            # append compound information
            iterate_rs_components(
                rs,
                molecule_dataset=molecule_dataset
            )
        # pickle reaction system object to file
        # prefix sRS + EC + DB + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
        count += 1


def get_rxn_system(rs, ID):
    """Get reaction system from SABIO reaction ID (rID).

    Uses SABIO API - online.

    Keywords:
        rs (class) - reaction system object
        ID (str) - DB reaction ID

    """
    QUERY_URL = (
        'http://sabiork.h-its.org/testSabio/'
        'sabioRestWebServices/searchReactionParticipants'
    )
    # input: SabioReactionID
    # valid output fields: "fields[]":
    #    ["Name","Role","SabioCompoundID","ChebiID",
    #     "PubChemID","KeggCompoundID","InChI"]

    params = {
        "SabioReactionID": ID,
        "fields[]": [
            "Name", "Role", "SabioCompoundID", "ChebiID",
            "PubChemID", "KeggCompoundID", 'UniprotID'
        ]
    }
    request = requests.post(QUERY_URL, params=params)
    request.raise_for_status()
    if request.text == 'No results found for query':
        rs.skip_rxn = True
        rs.skip_reason = 'No results for SABIO R-ID'
        return rs
    # collate request output
    rs.components = []
    fail_list = fail_list_read(
        directory='/home/atarzia/psp/molecule_DBs/atarzia/',
        file_name='failures.txt'
    )
    for i in request.text.split('\n')[1:]:
        if len(i) > 1:
            mol, role, cID, chebiID, pubchemID, keggID, _ = i.split(
                '\t'
            )
            # due to a bug with SABIO - we need to skip all rxns with
            # 'DNA'
            if mol == 'DNA':
                rs.skip_rxn = True
                rs.skip_reason = 'DNA present - SABIO has a bug'
                break
            # if 'inhibitor' in role.lower() or 'activator' in
            # role.lower() or
            #       'unknown' in role.lower():
            if 'modifier' in role.lower():
                # as of 23/11/18 we are not including any modifiers
                # including co factors/activators/inhibitors and
                # unknown
                continue
            # check if component name should be changed to a common
            # name
            mol, role = check_arbitrary_names((mol, role))
            if mol in fail_list:
                rs.skip_rxn = True
                rs.skip_reason = 'one component failed resolution'
                print('one molecule in fail list - skipping...')
                print('>>> ', mol, 'failed')
                break
            new_mol = molecule(mol, role, 'SABIO', cID)
            new_mol.PubChemID = None
            new_mol.chebiID = chebiID.rstrip().lstrip()
            new_mol.KEGG_ID = keggID.rstrip().lstrip()
            print('name >', mol, '--- role >', new_mol.role)
            print('KEGG ID >', new_mol.KEGG_ID)
            print('CHEBI ID >', new_mol.chebiID)
            # attempt to find chebiID
            if new_mol.chebiID == '' or ' ' in new_mol.chebiID:
                # try with KEGG ID first
                test1 = new_mol.KEGG_ID != ''
                if test1 and ' ' not in new_mol.KEGG_ID:
                    # check for CHEBI ID in KEGG translations file
                    translated = check_translator(new_mol.KEGG_ID)
                    if translated is not None:
                        pkl = translated
                        print(
                            'collecting molecule using translator:',
                            new_mol.KEGG_ID
                        )
                        new_mol = load_molecule(pkl, verbose=True)
                        new_mol.KEGG_ID = keggID
                        # need to make sure tranlated molecule has the
                        # correct role
                        new_mol.role = role
                        new_mol.translated = True
                    else:
                        new_mol.chebiID = KEGGID_to_CHEBIID(
                            KEGG_ID=new_mol.KEGG_ID
                        )
                        new_mol.translated = False
                else:
                    # KEGG ID is not unambiguously defined.
                    # then try with libchebipy
                    print(
                        '>> re-searching for chebiID using libchebipy'
                    )
                    # search CHEBI using molecule name
                    chebiID = CHEBI_IO.get_chebiID(mol)
                    new_mol.chebiID = chebiID
            # add new_mol to reaction system class
            rs.components.append(new_mol)
    # import sys
    # sys.exit()
    return rs
