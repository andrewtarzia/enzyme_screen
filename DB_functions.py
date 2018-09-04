#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O with molecule output files and DB handling.

Author: Andrew Tarzia

Date Created: 30 Aug 2018
"""

import pandas as pd
from numpy import average
import rdkit_functions


def initialize_mol_output_DF(filename, overwrite=False):
    """Read in or overwrite molecule output file to Pandas DF.

    """
    if overwrite is True:
        # pandas data frame for output
        molecule_output = pd.DataFrame(columns=['name', 'iupac_name',
                                                'DB', 'DB_ID', 'SMILE', 'role',
                                                'min_diam', 'mid_diam',
                                                'max_diam', 'ratio_1',
                                                'ratio_2'])
        molecule_output.to_csv(filename, index=False)
    else:
        try:
            molecule_output = pd.read_csv(filename)
        except FileNotFoundError:
            # pandas data frame for output
            molecule_output = pd.DataFrame(columns=['name', 'iupac_name',
                                                    'DB', 'DB_ID', 'SMILE', 'role',
                                                    'min_diam', 'mid_diam',
                                                    'max_diam', 'ratio_1',
                                                    'ratio_2'])
            molecule_output.to_csv(filename, index=False)

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


def get_molecule_diameters(mol_dict, molecule_output, mol_output_file, db_dir,
                           vdwScale=1.0,
                           boxMargin=4.0,
                           spacing=1.0,
                           N_conformers=10,
                           MW_tresh=130):
    """Get the molecule diameters of molecules in a dictionary.

    Role is reactant or product.

    """
    for key, val in mol_dict.items():
        if val[0] == '-':
            # check if key already output
            if key in list(molecule_output['name']):
                continue
            out_row = pd.DataFrame([
                [key, val[3], val[1], val[2], val[0],  val[4],
                 0, 0, 0, 0, 0]],
                columns=molecule_output.columns)
            # append row to molecule_output
            molecule_output = molecule_output.append(out_row,
                                                     ignore_index=True)
        # check if calculation already done
        # collect results if so
        if key in list(molecule_output['name']):
            res_line = molecule_output[molecule_output['name'] == key]
            old_role = res_line['role'].iloc[0]
            # if previous calculation was for a different role
            # then modify the existing role to be 'both'
            if val[4] != old_role and old_role != 'both':
                res_line['role'] = 'both'
            # update line
            molecule_output[molecule_output['name'] == key] = res_line

        # check IUPAC name column also
        elif key in list(molecule_output['iupac_name']):
            res_line = molecule_output[molecule_output['iupac_name'] == key]
            old_role = res_line['role'].iloc[0]
            # if previous calculation was for a different role
            # then modify the existing role to be 'both'
            if val[4] != old_role and old_role != 'both':
                res_line['role'] = 'both'
            # update line
            molecule_output[molecule_output['iupac_name'] == key] = res_line

        else:
            print('doing calculation...')
            # name: smiles
            molecule = {key: val[0]}
            res = rdkit_functions.calc_molecule_diameters(molecule,
                                                          out_dir=db_dir,
                                                          vdwScale=vdwScale,
                                                          boxMargin=boxMargin,
                                                          spacing=spacing,
                                                          N_conformers=N_conformers,
                                                          MW_tresh=MW_tresh)
            if res is None:
                out_row = pd.DataFrame([
                    [key, val[3], val[1], val[2], val[0], val[4],
                     0, 0, 0, 0, 0]],
                    columns=molecule_output.columns)
            else:
                # get the min values of all diameters of all conformers
                min_diam = min(res['diam1'])
                mid_diam = min(res['diam2'])
                max_diam = min(res['diam3'])
                # get avg values of all ratios of all conformers
                ratio_1 = average(res['ratio_1'])
                ratio_2 = average(res['ratio_2'])

                out_row = pd.DataFrame([
                    [key, val[3], val[1], val[2], val[0], val[4],
                     min_diam, mid_diam, max_diam, ratio_1, ratio_2]],
                    columns=molecule_output.columns)

            # append row to molecule_output
            molecule_output = molecule_output.append(out_row,
                                                     ignore_index=True)

        # update molecule output file
        save_mol_output_DF(mol_output_file, molecule_output)
