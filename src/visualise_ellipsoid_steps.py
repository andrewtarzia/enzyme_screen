#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to populate the properties of all molecules in database.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import sys

import rdkit_functions as rdkf
import utilities


def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: molecule_population.py param_file redo mol_file
    param_file:
    molecule :
        molecule file (_unopt.mol) to show the ellipsoid of.
""")
        sys.exit()
    else:
        params = utilities.read_params(sys.argv[1])
        molecule = sys.argv[2]

    vdwScale = params['vdwScale']
    boxMargin = params['boxMargin']
    spacing = params['spacing']
    show_vdw = params['show_vdw']
    plot_ellip = params['plot_ellip']
    # Set conformers to 1 to speed up.
    N_conformers = 1
    MW_thresh = params['MW_thresh']
    seed = int(params['seed'])

    name = molecule.replace('_unopt.mol', '')
    diam_file = f'ellips_plots/{name}_diam.csv'
    smiles = rdkf.read_structure_to_smiles(molecule)

    print('>> getting molecular size')
    _ = rdkf.calc_molecule_diameter(
        name,
        smiles,
        out_file=diam_file,
        vdwScale=vdwScale,
        boxMargin=boxMargin,
        spacing=spacing,
        MW_thresh=MW_thresh,
        show_vdw=show_vdw,
        plot_ellip=plot_ellip,
        N_conformers=N_conformers,
        rSeed=seed,
        do_step_plot=True
    )
    del _


if __name__ == "__main__":
    main()
