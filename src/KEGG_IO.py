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
from os.path import exists, join

from rdkit_functions import MolBlock_to_SMILES
from reaction import Reaction
import molecule
import IO


class KEGG_Reaction(Reaction):
    """
    Class that defines a reaction system for KEGG database.

    """
    def __init__(self, EC, DB, DB_ID, params):

        super().__init__(EC, DB, DB_ID, params)

    def get_equation(self):
        # get Reaction information from KEGG API
        URL = 'http://rest.kegg.jp/get/reaction:'+self.DB_ID
        request = requests.post(URL)
        if request.text != '':
            request.raise_for_status()
        else:
            self.skip_rxn = True
            self.skip_reason = (
                'No result for KEGG URL search - likely outdated'
            )
        # because of the formatting of KEGG text - this is trivial
        self.equation = request.text.split('EQUATION    ')[1]
        self.equation = self.equation.split('\n')[0].rstrip()

    def get_components(self):
        """
        Get reactants and products from a KEGG equation string.

        Returns bool (reversible), list of components.

        """

        if '<=>' in self.equation:
            # implies it is reversible
            reactants, products = self.equation.split("<=>")
            self.reversible = True
        else:
            reactants, products = self.equation.split("=")
            self.reversible = False
        reactants = [i.lstrip().rstrip() for i in reactants.split("+")]
        products = [i.lstrip().rstrip() for i in products.split("+")]

        # collect KEGG Compound/Glycan ID
        # check if reactant or product are compound or glycan
        # remove stoichiometry for now
        self.comp_list = []
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
            self.comp_list.append((KID, 'reactant'))
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
            self.comp_list.append((KID, 'product'))

    def get_rxn_system(self):
        """
        Get reaction system from KEGG reaction ID (rID).

        Use KEGG API - online.

        Keywords:
            rs (class) - reaction system object
            ID (str) - DB reaction ID

        """

        # collate request output
        fail_list = IO.fail_list_read(
            directory=self.params['molec_dir'],
            file_name='failures.txt'
        )
        print(fail_list)

        self.get_equation()
        print(self.equation)
        self.components = []
        self.get_components()
        print(self.reversible, self.comp_list)
        for comp in self.comp_list:
            print(comp)
            molec_struct_name = (
                f"{self.params['molec_dir']}"
                f'/{comp[0]}_unopt.mol'
            )
            print(molec_struct_name)

            # handle polymeric species:
            if '(n)' in comp[0]:
                # implies polymer
                self.fail_polymer()
                # write KEGG ID to fail list
                IO.fail_list_write(
                    new_name=comp[0],
                    directory=self.params['molec_dir'],
                    file_name='failures.txt'
                )
                print('failed')
                import sys
                sys.exit()
            elif comp[0] in fail_list:
                self.fail_fail_list(comp[0])
                break

            # No failures yet. Collect structure.
            new_mol = molecule.Molecule(
                name=comp[0],
                role=comp[1],
                DB='KEGG',
                DB_ID=comp[0],
                params=self.params,
                structure_file=molec_struct_name
            )
            new_mol.KEGG_ID = comp[0]
            print(new_mol)
            print(new_mol.__dict__)
            if exists(new_mol.structure_file):
                new_mol.SMILES = new_mol.read_structure_to_smiles()
                self.components.append(new_mol)
            else:
                result = self.component_KEGGID_to_MOL(KEGG_ID=comp[0])
                if result is None:
                    self.fail_conversion(comp[0])
                    # write KEGG ID to fail list
                    IO.fail_list_write(
                        new_name=comp[0],
                        directory=self.params['molec_dir'],
                        file_name='failures.txt'
                    )
                elif result == 'generic':
                    self.fail_generic(comp[0])
                    # write KEGG ID to fail list
                    IO.fail_list_write(
                        new_name=comp[0],
                        directory=self.params['molec_dir'],
                        file_name='failures.txt'
                    )
                else:
                    # add new_mol to reaction system class
                    new_mol.SMILES = result[1]
                    new_mol.write_structure(result[0])
                    self.components.append(new_mol)
            print(new_mol)
            input()

    def component_KEGGID_to_MOL(self, KEGG_ID):
        """
        Use KEGG API to collect MOL file using KEGG API.

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
            result = MolBlock_to_SMILES(output)
            return result
        elif request.status_code == 400:
            # implies a bad request/syntax
            print(KEGG_ID, 'has bad syntax')
            return None
        else:
            print("haven't come across this yet.")
            raise ValueError(
                f'KEGG ID: {KEGG_ID} gives this error: '
                f'{request.status_code}'
            )

    def iterate_rs_components(self):
        """
        Iterate over all components in a reaction system and collect
        the component structures/properties.

        > This is a version of the same function in molecule.py
        (iterate_rs_components) that was designed to handle
        KEGG collected structures implemented 12/12/18.

        All changes to the reaction system should be done in-place.

        Arguments:
            rs (rxn_syst.reaction) - reaction system being tested
            molecule_dataset (Pandas DataFrame) -
                look up for known molecules

        """
        for m in self.components:
            print('>', m)
            # check for wildcard in SMILES
            if '*' in m.SMILES:
                self.fail_generic_smiles(m.KEGG_ID)
                # write KEGG ID to fail list
                IO.fail_list_write(
                    new_name=m.KEGG_ID,
                    directory=self.params['molec_dir'],
                    file_name='failures.txt'
                )
                break

            m.get_properties()
            print('lp', m.logP, 'ls', m.logS, 'SA', m.Synth_score)

        print([j.logP for j in self.components])
        input()


def check_translator(ID, params):
    """
    Check for KEGG ID in KEGG translation file.

    Links to molecules already collected from online.

    """

    # read translator
    with gzip.GzipFile(params['translator_kegg'], 'rb') as output:
        translation = pickle.load(output)
    try:
        return translation[ID]
    except KeyError:
        return None


def get_EC_rxns_from_JSON(JSON_DB, EC):
    """
    Get reactions associated with EC number in KEGG.

    """
    try:
        rxn_list = JSON_DB[EC]
        return rxn_list
    except KeyError:
        return None


def get_rxn_systems(
    EC,
    output_dir,
    params,
    clean_system=False,
    verbose=False
):
    """
    Get reaction systems from KEGG entries in one EC and output to
    Pickle.

    """

    # read in JSON file of whole DB
    rxn_DB_file = params['KEGG_json_ec_file']
    with open(rxn_DB_file, 'r') as data_file:
        rxn_DB = json.load(data_file)

    # get EC specific entries
    EC_rxns = get_EC_rxns_from_JSON(rxn_DB, EC)
    if EC_rxns is None:
        return None

    # iterate over reactions
    count = 0
    for rxn in EC_rxns:
        print(rxn)
        # get KEGG rxn id
        string = rxn['name']
        K_Rid = string.split(' ')[0].rstrip()
        # initialise reaction system object
        rs = KEGG_Reaction(EC, 'KEGG', K_Rid, params)
        print(rs)

        if exists(join(output_dir, rs.pkl)) and clean_system is False:
            count += 1
            continue
        if verbose:
            print('=================================================')
            print(
                'DB: KEGG - EC:', EC, '-',
                'DB ID:', K_Rid, '-', count, 'of', len(EC_rxns)
            )

        # get reaction system using DB specific function
        rs.get_rxn_system()
        print('-----------------')
        print(rs)
        print(rs.components)
        print(rs.skip_rxn)
        if not rs.skip_rxn:
            # append compound information
            rs.iterate_rs_components()
        # pickle reaction system object to file
        rs.save_reaction(join(output_dir, rs.pkl))
        count += 1

        print('--- updating done file ---')
        done_file = 'done_RS.txt'
        if exists(done_file):
            with open(done_file, 'a') as f:
                f.write(self.pkl+'\n')
        else:
            with open(done_file, 'w') as f:
                f.write('pkl\n'+self.pkl+'\n')
