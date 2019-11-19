#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Utility functions.

Author: Andrew Tarzia

Date Created: 19 Nov 2019

"""


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
            val = True if val =='True' else False
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
    print('---------------------------------------------------------')
