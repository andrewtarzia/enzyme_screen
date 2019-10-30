#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Define globally used parameters.
"""


def get_parameters():
    """
    Get parameter dictionary
    """
    param_dict = {
        # where all the brenda data files are
        'BRENDA_DB_loc': '/home/atarzia/psp/molecule_DBs/brenda_details/',
        'out_CSV_pi': "output_data_pi.csv",
        'out_columns_pi': ['fasta_file', 'acc.code',
                           'organism', 'EC.code', 'species',
                           'note', 'pi', 'modification', 'category'],
        'out_CSV_br': "output_data_br.csv",
        'out_columns_br': ['fasta_file'],
        # cut off for ZIF growth from pI
        'cutoff_pi': 6,
        # modification types + colours
        # just implement succinylation for now
        # why would you want to do acetylation if you can succinylate??
        # currently succinylation == LYS swapped with GLU
        'modifications': {
            '0': {
                'colour': 'k',
                'name': 'unmodified',
            },
            '1': {
                'colour': 'firebrick',
                'name': 'succinylated',
                'target_res': 'LYS',
                'replace_res': 'GLU',
            }
        },
        'diffuse_threshold': 4.2,

    }

    return param_dict
