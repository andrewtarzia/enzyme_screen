#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for linB reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import time
import sys

import utilities
import rdkit_functions as rdkf
import plotting_fn as pfn


def main():
    if (not len(sys.argv) == 4):
        print("""
    Usage: linB_screening.py

        molecule_file

        rerun_diameter_calc

        param_file

        """)
        sys.exit()
    else:
        molecule_file = sys.argv[1]
        rerun_diameter_calc = True if sys.argv[2] == 't' else False
        pars = utilities.read_params(sys.argv[3])

    start = time.time()

    df, molecules, diameters = rdkf.read_mol_txt_file(
        molecule_file
    )

    # draw 2D structures
    print('--- draw 2D structures...')
    rdkf.draw_svg_for_all_molecules(molecules)

    # calculate the size of the ellipsoid surroudning all molecules
    # using input pars
    if rerun_diameter_calc:
        print('--- calculating molecular diameters...')
        rdkf.calc_molecule_diameters(
            molecules,
            pars=pars,
            out_dir='linB_pars',
        )

    # print results for each molecule
    print('--- print results and plot...')
    pfn.print_results(
        molecules,
        threshold=pars['size_thresh'],
        output_dir='linB_pars'
    )

    pfn.categorical(
        molecules,
        threshold=pars['size_thresh'],
        output_dir='linB_pars',
        plot_suffix='linB'
    )
    pfn.shapes(
        molecules,
        threshold=pars['size_thresh'],
        output_dir='linB_pars',
        plot_suffix='linB'
    )

    end = time.time()
    print(f'---- total time taken = {round(end-start, 2)} s')


if __name__ == "__main__":
    main()
