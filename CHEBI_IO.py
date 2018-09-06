#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of CHEBI DB.

Using downloaded flat files. Offline.

Author: Andrew Tarzia

Date Created: 30 Aug 2018


"""

import pandas as pd
import DB_functions
from rdkit.Chem import AllChem as Chem


def search_for_compound_by_name(file, cmpd):
    """Search compounds.tsv file for matching name.

    """
    # read file line by line
    # header = 'ID', 'STATUS', 'CHEBI_ACCESSION',
    #          'SOURCE', 'PARENT_ID', 'NAME',
    #          'DEFINITION', 'MODIFIED_ON', 'CREATED_BY', 'STAR'
    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if name == cmpd:
                return ID, parent_id, name, star
            elif name.lower() == cmpd.lower():
                return ID, parent_id, name, star

    return None


def search_for_compound_by_id(file, chebiID):
    """Search compounds.tsv file for matching name.

    """
    # read file line by line
    # header = 'ID', 'STATUS', 'CHEBI_ACCESSION',
    #          'SOURCE', 'PARENT_ID', 'NAME',
    #          'DEFINITION', 'MODIFIED_ON', 'CREATED_BY', 'STAR'
    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if ID == str(chebiID):
                return ID, parent_id, name, star
    return None


def search_for_name_by_id(file, chebiID):
    """Search names.tsv file for matching CHEBI ID to get name.

    """
    # read file line by line
    # header = 'ID', 'COMPOUND_ID', 'TYPE',
    #          'SOURCE', 'NAME', 'ADAPTED', 'LANGUAGE'

    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, C_ID, _, _, name, _, _ = line_split
            if ID == 'ID':
                continue
            if C_ID == str(chebiID):
                return C_ID, name

    return None


def search_for_name_by_name(file, cmpd):
    """Search names.tsv file for matching name to get CHEBI ID.

    """
    # read file line by line
    # header = 'ID', 'COMPOUND_ID', 'TYPE',
    #          'SOURCE', 'NAME', 'ADAPTED', 'LANGUAGE'

    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, C_ID, _, _, name, _, _ = line_split
            if ID == 'ID':
                continue
            if name == cmpd:
                return C_ID, name
            elif name.lower() == cmpd.lower():
                return C_ID, name

    return None


def convert_nameID_to_parent(file, nameID):
    """Make sure ID extracted for the name is parent - get CHEBI ID of parent
    if not.

    """
    # read file line by line
    # header = 'ID', 'STATUS', 'CHEBI_ACCESSION',
    #          'SOURCE', 'PARENT_ID', 'NAME',
    #          'DEFINITION', 'MODIFIED_ON', 'CREATED_BY', 'STAR'
    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if nameID == ID:
                if parent_id != 'null':
                    desired_parent_id = parent_id
                    break
                else:
                    return ID, parent_id, name, star

    # now search for parent if necessary
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if ID == desired_parent_id:
                if parent_id == 'null':
                    return ID, parent_id, name, star

    return None


def get_structure(file, C_ID):
    """Get structure from CHEBI structures.csv file.

    """
    # use pandas (in chunks to avoid memory issues)
    # because it is a nice CSV file.
    structs = pd.read_csv(file, chunksize=1000)

    for chunk in structs:
        CID_frame = chunk[chunk['COMPOUND_ID'] == int(C_ID)]
        if len(CID_frame) > 0:
            # may be a 2D or 3D structure
            # output in hierarchy
            types = list(CID_frame['TYPE'])
            if 'SMILES' in types:
                row = CID_frame[CID_frame['TYPE'] == 'SMILES']
                structure = row['STRUCTURE'].iloc[0]
                s_type = row['TYPE'].iloc[0]
            elif 'InChI' in types:
                row = CID_frame[CID_frame['TYPE'] == 'InChI']
                structure = row['STRUCTURE'].iloc[0]
                s_type = row['TYPE'].iloc[0]
            elif 'mol' in types:
                row = CID_frame[CID_frame['TYPE'] == 'mol']
                structure = row['STRUCTURE'].iloc[0]
                s_type = row['TYPE'].iloc[0]
            else:
                structure = None
                s_type = None
            return structure, s_type

    return None, None


def get_cmpd_information(molec):
    """Get information from CHEBI Database of a compound from CHEBI ID.

    Done Offline. molec must have attribute 'chebiID' as integer.

    """

    DB_prop = DB_functions.get_DB_prop('CHEBI')
    compounds_file = DB_prop[0]+DB_prop[1]['cmpds_file']
    names_file = DB_prop[0]+DB_prop[1]['names_file']
    structures_file = DB_prop[0]+DB_prop[1]['strct_file']

    # set name by searching compound file
    res = search_for_compound_by_id(compounds_file, molec.chebiID)
    if res is None:
        print('no match in DB - this should not happen for CHEBI ID search')
        print('check this!')
        print('Exitting....')
        import sys
        sys.exit()
    else:
        ID, parent_id, name, star = res
        molec.name = name
        molec.change_name = False

    # make sure is parent compound
    if parent_id != 'null':
        res = convert_nameID_to_parent(compounds_file, nameID=ID)
        if res is None:
            print("this should not happen - error with cross reference")
            print('check this!')
            print('Exitting....')
            import sys
            sys.exit()
        ID, parent_id, name, star = res
        molec.name = name
        molec.change_name = False
        molec.chebiID = int(ID)

    # get structure using CHEBI ID
    # structures.csv - read in, get COMPOUND ID match then extract the
    # get SMILES
    structure, s_type = get_structure(structures_file, molec.chebiID)
    if structure is not None:
        # is structure a MolBlock or Smiles
        if s_type == 'mol':
            # convert structure to SMILEs
            rdkitmol = Chem.MolFromMolBlock(structure)
            if rdkitmol is None:
                print('structure could not be deciphered')
                smile = '-'
                molec.SMILES = smile
                molec.mol = None
                print('probably a polymeric structure - skipping.')
            else:
                rdkitmol.Compute2DCoords()
                smile = Chem.MolToSmiles(rdkitmol)
                molec.mol = rdkitmol
                molec.SMILES = smile
        elif s_type == 'SMILES':
            smile = structure
            rdkitmol = Chem.MolFromSmiles(smile)
            rdkitmol.Compute2DCoords()
            molec.SMILES = smile
            molec.mol = rdkitmol
        elif s_type == 'InChIKey' or s_type == 'InChI':
            rdkitmol = Chem.MolFromInchi(structure)
            rdkitmol.Compute2DCoords()
            smile = Chem.MolToSmiles(rdkitmol)
            molec.mol = rdkitmol
            molec.SMILES = smile
    else:
        smile = '-'
        molec.SMILES = smile
        molec.mol = None
        print('molecule does not have recorded structure in CHEBI DB')
        print('probably a generic structure - skipping.')
