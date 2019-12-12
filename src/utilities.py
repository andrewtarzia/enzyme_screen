#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utility functions.

Author: Andrew Tarzia

Date Created: 19 Nov 2019

"""

import pandas as pd


def get_ECs_from_file(EC_file):
    """
    Read in ECs to search from a file.

    """
    # get search EC numbers from file:
    # set EC numbers of interest
    # get from a data file - manually made from
    # https://enzyme.expasy.org/enzyme-byclass.html
    EC_DF = pd.read_table(
        EC_file,
        delimiter='__',
        names=['EC_no'],
        engine='python'
    )
    search_ECs = list(EC_DF['EC_no'])

    # remove all spaces within EC numbers
    search_ECs = [i.replace(' ', '') for i in search_ECs]

    # add check for '1' from '1.-.-.-'
    new_search_ECs = []
    for EC in search_ECs:
        if '-' in EC:
            new_search_ECs.append(EC.replace('.-', ''))
            new_search_ECs.append(EC)
        else:
            new_search_ECs.append(EC)

    print('>>', len(search_ECs), 'EC numbers to test')
    print('>> first EC:', search_ECs[0], '-> last EC:', search_ECs[-1])
    return new_search_ECs


def read_params(file):
    """
    Read parametrs for screening.

    Returns
    -------
    params : :class:`dict`
        Dictionary of parameters.

    """
    params = {}
    with open(file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        key, val = line.rstrip().split(':')
        if val == 'True' or val == 'False':
            val = True if val == 'True' else False
        else:
            try:
                val = float(val)
            except ValueError:
                pass
        params[key] = val

    print_params(params)

    return params


def print_params(pars):
    print('---------------------------------------------------------')
    print('run parameters:')
    print('VDW scale:', pars['vdwScale'])
    print('Box Margin:', pars['boxMargin'], 'Angstrom')
    print('Grid spacing:', pars['spacing'], 'Angstrom')
    print('show VDW?:', pars['show_vdw'])
    print('Plot Ellipsoid?:', pars['plot_ellip'])
    print('No Conformers:', pars['N_conformers'])
    print('MW threshold:', pars['MW_thresh'], 'g/mol')
    print('pI threshold:', pars['pI_thresh'])
    print('Diffusion threshold:', pars['size_thresh'], 'Angstrom')
    print('DBs:', pars['DBs'])
    print('---------------------------------------------------------')
