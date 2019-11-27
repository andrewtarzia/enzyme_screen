#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to populate the properties of all molecules in database.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os.path import exists
import sys
import glob
import json
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

import IO
import rdkit_functions as rdkf
import plotting_fn
import utilities


def populate_all_molecules(params, redo, mol_file=None):
    """
    Populate all molecules in pickle files in directory.

    """

    vdwScale = params['vdwScale']
    boxMargin = params['boxMargin']
    spacing = params['spacing']
    show_vdw = params['show_vdw']
    plot_ellip = params['plot_ellip']
    N_conformers = int(params['N_conformers'])
    MW_thresh = params['MW_thresh']
    seed = int(params['seed'])

    fail_list = IO.fail_list_read(
        directory=params['molec_dir'],
        file_name='failures.txt'
    )

    if mol_file is None:
        molecule_list = glob.glob('*_unopt.mol')
    else:
        molecule_list = IO.read_molecule_list(mol_file)

    print(f'{len(molecule_list)} molecules in DB.')

def main():
    if (not len(sys.argv) == 5):
        print("""
Usage: molecule_population.py param_file redo mol_file
    param_file:
    redo:
        t to overwrite all rxn systems.
    plot:
        t to plot distributions of molecule properties.
    mol_file :
        file name of list of molecules to allow for trivial
        parallelisation, `f` if not specified, where all `mol`
        files are populated.
""")
        sys.exit()
    else:
        params = utilities.read_params(sys.argv[1])
        redo = True if sys.argv[2] == 't' else False
        plot = True if sys.argv[3] == 't' else False
        mol_file = None if sys.argv[4] == 'f' else sys.argv[4]

    if input('do plots? (t/f)') == 't':
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        #######
        # molecule distributions
        #######
        # plot synthetic accessibility VS no heavy atoms
        mol_logP_vs_logS(output_dir=directory, plot_suffix='mol_cf')
        mol_logP_vs_XlogP(output_dir=directory, plot_suffix='mol_cf')
        # plot synthetic accessibility VS complexity
        mol_SA_vs_compl(output_dir=directory, plot_suffix='mol_cf')
        # plot synthetic accessibility VS no rotatable bonds
        mol_SA_vs_NRB(output_dir=directory, plot_suffix='mol_cf')
        # # plot synthetic accessibility VS no heavy atoms
        mol_SA_vs_NHA(output_dir=directory, plot_suffix='mol_cf')
        # # plot distributions of attributes
        mol_all_dist(output_dir=directory, plot_suffix='mol_cf')
    sys.exit()
    print('settings:')
    print('    Molecule file:', mol_file)
    print(
        'populate the properties attributes for all '
        'molecules in DB...'
    )

    populate_all_molecules(params=params, mol_file=mol_file, redo=redo)



if __name__ == "__main__":
    main()
