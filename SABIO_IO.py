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
# ensure cpickle usage
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


class SABIO_reaction:
    """Class that defines a reaction system within SABIO Database.

    """

    def __init__(self, EC, eID):
        self.EC = EC
        self.eID = eID  # SABIO entry ID
        self.pkl_name = 'sRS-'+EC.replace('.', '_')+'-'+str(eID)+'.pkl'
        self.rID = None  # SABIO reaction ID
        self.organism = None
        self.components = None  # molecular components
        self.UniprotID = None  # Uniprot ID for protein sequence of reaction
        self.all_fit = None  # do all the components fit?
        # None implies that it is ambiguous.
        self.seed_MOF = None  # will the protein sequence seed MOF growth
        self.req_mod = None  # does it require modification to seed MOF growth

    def check_all_fit(self, threshold, molecule_output):
        """Check if all components of reaction system fit.

        """
        all_fit = True
        max_comp_size = 0
        for r in self.components:
            # is molecule in molecule_output?
            if r.name in list(molecule_output['name']):
                r_diam = float(molecule_output[molecule_output['name'] == r.name]['mid_diam'])
                max_comp_size = max([r_diam, max_comp_size])
                if r_diam > threshold or r_diam == 0:
                    all_fit = False
            else:
                # molecule not in output - assume it wasn't in database
                # dont report!
                all_fit = False
                break

        if all_fit is True:
            self.all_fit = True
            self.max_comp_size = max_comp_size
            print('------------------------------------')
            self.print_rxn_system
            print("This reaction should diffuse!")
            print('------------------------------------')
        else:
            self.all_fit = False
            self.max_comp_size = max_comp_size

    def save_object(self, filename):
        """Pickle reaction system object to file.

        """
        # Overwrites any existing file.
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(self, filename):
        """unPickle reaction system object from file.

        """
        print('loading:', filename)
        with open(filename, 'rb') as input:
            self = pickle.load(input)
            return self

    def print_rxn_system(self):
        """Fancy print of reaction system.

        """
        print('--------------------------')
        print('EC:', self.EC)
        print('Organism:', self.organism)
        print('SABIO entry ID:', self.eID)
        print('SABIO reaction ID:', self.rID)
        print('Uniprot ID:', self.UniprotID)
        print('--------------------------')
        if self.components is not None:
            for i in self.components:
                print(i.name, ' (ID:', i.cID+') as', i.role)
            print('--------------------------')
        if self.all_fit is True:
            print('all components will diffuse through')
            print('--------------------------')

    def get_reaction_system(self):
        """Get reaction system from SABIO reaction ID (rID).

        """
        QUERY_URL = 'http://sabiork.h-its.org/testSabio/sabioRestWebServices/searchReactionParticipants'
        # input: SabioReactionID
        # valid output fields: "fields[]":
        #    ["Name","Role","SabioCompoundID","ChebiID",
        #     "PubChemID","KeggCompoundID","InChI"]

        params = {"SabioReactionID": self.rID,
                  "fields[]": ["Name", "Role", "SabioCompoundID"]}
        # , "ChebiID",
        # "PubChemID", "KeggCompoundID", 'UniprotID']}
        request = requests.post(QUERY_URL, params=params)
        request.raise_for_status()
        # collate request output
        self.components = []
        for i in request.text.split('\n')[1:]:
            if len(i) > 1:
                mol, role, cID = i.split('\t')
                new_mol = SABIO_molecule(mol, role, cID)
                # add new_mol to reaction system class
                self.components.append(new_mol)


class SABIO_molecule:
    """Class that defines the molecules extracted from SABIO database.

    """

    def __init__(self, name, role, cID):
        self.name = name
        self.role = role
        self.cID = cID
        self.InChi = None

    def get_cmpd_information(self):
        """Get information from SABIO Database of a compound with ID cID.

        """
        QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchCompoundDetails'

        # input: SabioCompoundID
        # valid output fields: "fields[]":["Name","ChebiID",
        #                                  "PubChemID","InChI",
        #                                  "SabioCompoundID","KeggCompoundID"]

        params = {"SabioCompoundID": self.cID,
                  "fields[]": ["Name", "ChebiID", "PubChemID", "InChI"]}
        if self.InChi is None:
            request = requests.post(QUERY_URL, params=params)
            request.raise_for_status()
            # results
            txt = request.text.split('\n')[1].split('\t')
            _, self.ChebiID, self.PubChemId, self.InChi = txt

        if self.InChi != 'null':
            self.mol = get_rdkit_mol_from_InChi(self.InChi)
        else:
            self.mol = None


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
    # results
    _, organism, _, rxn_id, UniprotID = request.text.split('\n')[1].split('\t')
    return organism, rxn_id, UniprotID
