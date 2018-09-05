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


def get_DB_prop(DB):
    """Returns properties of specific database.

    {DB: (directory, DB specific dictionary)}

    """

    DBs = {
         'CHEBI': (
            '/home/atarzia/psp/molecule_DBs/chebi/',
            {'cmpds_file': 'compounds.tsv',
             'names_file': 'names.tsv',
             'strct_file': 'structures.csv'}
                  ),
         'BKMS': (
            '/home/atarzia/psp/molecule_DBs/BKMS_react/', {}
                  ),
         'SABIO': (
            '/home/atarzia/psp/molecule_DBs/SABIO/', {}
                   ),
         'KEGG': (
            '/home/atarzia/psp/molecule_DBs/KEGG/',
            {'JSON_file': 'br08201.json'}
                  ),

          }
    try:
        return DBs[DB]
    except KeyError:
        print("Error: This DB does not exist.")
        print('available DBs:')
        print(list(DBs.keys()))


def get_molecule_diameters(mol_dict, molecule_output, mol_output_file, db_dir,
                           vdwScale=1.0,
                           boxMargin=4.0,
                           spacing=1.0,
                           N_conformers=10,
                           MW_thresh=130):
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
        # check by SMILES (not name) because they should be specific.
        # collect results if so
        if val[0] in list(molecule_output['SMILE']):
            res_line = molecule_output[molecule_output['SMILE'] == val[0]]
            old_role = res_line['role'].iloc[0]
            # if previous calculation was for a different role
            # then modify the existing role to be 'both'
            if val[4] != old_role and old_role != 'both':
                res_line['role'] = 'both'
            # update line
            molecule_output[molecule_output['SMILE'] == val[0]] = res_line

        # # check IUPAC name column also
        # elif key in list(molecule_output['iupac_name']):
        #     res_line = molecule_output[molecule_output['iupac_name'] == key]
        #     old_role = res_line['role'].iloc[0]
        #     # if previous calculation was for a different role
        #     # then modify the existing role to be 'both'
        #     if val[4] != old_role and old_role != 'both':
        #         res_line['role'] = 'both'
        #     # update line
        #     molecule_output[molecule_output['iupac_name'] == key] = res_line

        else:
            print('doing calculation...')
            # name: smiles
            res = rdkit_functions.calc_molecule_diameter(key, val[0],
                                                         out_dir=db_dir,
                                                         vdwScale=vdwScale,
                                                         boxMargin=boxMargin,
                                                         spacing=spacing,
                                                         N_conformers=N_conformers,
                                                         MW_thresh=MW_thresh)
            if res is None:
                out_row = pd.DataFrame([
                    [key, val[3], val[1], val[2], val[0], val[4],
                     0, 0, 0, 0, 0]],
                    columns=molecule_output.columns)
            else:
                # get the min values of all diameters of all conformers
                min_diam = round(min(res['diam1']), 3)
                mid_diam = round(min(res['diam2']), 3)
                max_diam = round(min(res['diam3']), 3)
                # get avg values of all ratios of all conformers
                ratio_1 = round(average(res['ratio_1']), 3)
                ratio_2 = round(average(res['ratio_2']), 3)

                out_row = pd.DataFrame([
                    [key, val[3], val[1], val[2], val[0], val[4],
                     min_diam, mid_diam, max_diam, ratio_1, ratio_2]],
                    columns=molecule_output.columns)

            # append row to molecule_output
            molecule_output = molecule_output.append(out_row,
                                                     ignore_index=True)

        # update molecule output file
        save_mol_output_DF(mol_output_file, molecule_output)
