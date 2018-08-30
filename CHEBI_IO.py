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


def search_for_compound(file, cmpd):
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


def search_for_name(file, cmpd):
    """Search names.tsv file for matching name.

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
            structure = CID_frame['STRUCTURE'].iloc[0]
            s_type = CID_frame['TYPE'].iloc[0]
            return structure, s_type

    return None, None
