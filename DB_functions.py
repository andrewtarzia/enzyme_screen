#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O with molecule output files and DB handling.

Author: Andrew Tarzia

Date Created: 30 Aug 2018
"""

import pandas as pd


def initialize_mol_output_DF(filename, overwrite=False):
    """Read in or overwrite molecule output file to Pandas DF.

    """
    if overwrite is True:
        # pandas data frame for output
        molecule_output = pd.DataFrame(columns=['EC', 'name', 'iupac_name',
                                                'DB', 'DB_ID', 'SMILE', 'role',
                                                'min_diam', 'mid_diam',
                                                'max_diam', 'ratio_1',
                                                'ratio_2'])
        molecule_output.to_csv(filename, index=False)
    else:
        molecule_output = pd.read_csv(filename)

    return molecule_output


def save_mol_output_DF(filename, molecule_output):
    """Save updated version of molecule output file from PANDAS DF.

    """
    molecule_output.to_csv(filename, index=False)


def get_DB_prop(DB=None):
    """Returns properties of available DBs.

    {DB: (directory, class initalized object specific to DB)}

    """

    DBs = {
         'CHEBI': ('/home/atarzia/psp/molecule_DBs/chebi/', '_'),

          }
    if DB is None:
        return DBs
    else:
        try:
            return {DB: DBs[DB]}
        except KeyError:
            print("Error: This DB does not exist.")
            print('available DBs:')
            print(list(DBs.keys()))
