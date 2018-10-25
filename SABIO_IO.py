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
import rxn_syst
from molecule import molecule, iterate_rs_components, check_arbitrary_names
import CHEBI_IO
import PUBCHEM_IO
from molvs import standardize_smiles


def get_cmpd_information(molec):
    """Get information from SABIO Database of a compound with ID cID.

    """
    QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchCompoundDetails'

    # input: SabioCompoundID
    # valid output fields: "fields[]":["Name","ChebiID",
    #                                  "PubChemID","InChI",
    #                                  "SabioCompoundID","KeggCompoundID"]
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
        smiles = Chem.MolToSmiles(Chem.RemoveHs(molec.mol))
        molec.SMILES = smiles
        try:
            molec.SMILES = standardize_smiles(molec.SMILES)
        except ValueError:
            print('standardization failed - therefore assume')
            print('SMILES were invalid - skip')
            molec.SMILES = None
            molec.mol = None
            import sys
            sys.exit()
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

    mol = Chem.MolFromInchi(InChi)
    if AddHs is True:
        mol = Chem.AddHs(mol)

    return mol


def get_entries_per_EC(EC):
    """Collect SABIO entry IDs associated with an EC no.

    Arguments:
        EC (str) - format: X.X.X.X

    Returns:
        entries (list) - list of entries

    """
    ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'

    # ask SABIO-RK for all EntryIDs matching a query
    query_dict = {"ECNumber": EC}
    query_string = ' AND '.join(['%s:%s' % (k, v)
                                 for k, v in query_dict.items()])
    query = {'format': 'txt', 'q': query_string}
    # make GET request
    request = requests.get(ENTRYID_QUERY_URL, params=query)
    request.raise_for_status()  # raise if 404 error
    if request.text == 'No results found for query':
        print(request.text)
        return []
    # each entry is reported on a new line
    entries = [int(x) for x in request.text.strip().split('\n')]
    print('%d matching entries found.' % len(entries))

    return entries


def get_rxnID_from_eID(eID):
    """Collect SABIO reaction ID and properties from entry ID.

    """
    PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'
    # encode next request, for parameter data given entry IDs
    data_field = {'entryIDs[]': eID}
    query = {'format': 'tsv',
             'fields[]': ['EntryID',
                          'Organism',
                          'ECNumber',
                          'SabioReactionID',
                          'UniprotID']}
    # make POST request
    request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
    request.raise_for_status()
    if request.text == 'No results found for query':
        print(request.text)
        return None, None, None
    # results
    _, organism, _, rxn_id, UniprotID = request.text.split('\n')[1].split('\t')
    return organism, rxn_id, UniprotID


def get_rxn_systems(EC, output_dir, molecule_dataset,
                    clean_system=False, verbose=False):
    """Get reaction systems from SABIO entries in one EC and output to Pickle.

    """
    # get all SABIO entries
    entries = get_entries_per_EC(EC)
    # iterate over entries
    count = 0
    for eID in entries:
        if verbose:
            print('DB: SABIO - EC:', EC, '-',
                  'DB ID:', eID, '-', count, 'of', len(entries))
        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'SABIO', eID)
        if os.path.isfile(output_dir+rs.pkl) is True and clean_system is False:
            count += 1
            continue
        # get reaction ID
        # SABIO specific properties
        rs.organism, rs.rID, rs.UniprotID = get_rxnID_from_eID(eID)
        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.rID)
        if rs.skip_rxn is False:
            # append compound information
            iterate_rs_components(rs, molecule_dataset=molecule_dataset)
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
    QUERY_URL = 'http://sabiork.h-its.org/testSabio/sabioRestWebServices/searchReactionParticipants'
    # input: SabioReactionID
    # valid output fields: "fields[]":
    #    ["Name","Role","SabioCompoundID","ChebiID",
    #     "PubChemID","KeggCompoundID","InChI"]

    params = {"SabioReactionID": ID,
              "fields[]": ["Name", "Role", "SabioCompoundID", "ChebiID",
                           "PubChemID", "KeggCompoundID", 'UniprotID']}
    request = requests.post(QUERY_URL, params=params)
    request.raise_for_status()
    if request.text == 'No results found for query':
        rs.skip_rxn = True
        return rs
    # collate request output
    rs.components = []
    for i in request.text.split('\n')[1:]:
        if len(i) > 1:
            mol, role, cID, chebiID, pubchemID, keggID, _ = i.split('\t')
            print(mol, role)
            # check if component name should be changed to a common name
            print('original name', mol)
            mol, role = check_arbitrary_names((mol, role))
            print('new name', mol)
            new_mol = molecule(mol, role, 'SABIO', cID)
            new_mol.PubChemID = None
            new_mol.chebiID = chebiID
            new_mol.KEGG_ID = keggID
            print(mol, new_mol.role)
            # check for multiple entries from SABIO
            if new_mol.chebiID == '' or ' ' in new_mol.chebiID:
                print('researching for chebiID using libchebipy...')
                # search CHEBI using molecule name
                chebiID = CHEBI_IO.get_chebiID(mol)
                new_mol.chebiID = chebiID
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
                        print('all failed - skipping...')
                        continue
            # add new_mol to reaction system class
            rs.components.append(new_mol)
    return rs
