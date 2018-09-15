#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for linB reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import rdkit_functions
import plotting

# script/data set specific functions


if __name__ == "__main__":
    # set parameters
    # molecule file dir
    molecule_file = '/home/atarzia/psp/linBmolecules/linbmolecules.txt'
    # output dir
    output_dir = '/home/atarzia/psp/linBmolecules/'
    vdwScale = 0.8
    boxMargin = 4.0
    spacing = 0.6
    show_vdw = False
    plot_ellip = False
    N_conformers = 50
    MW_thresh = 2000
    pI_thresh = 6
    size_thresh = 4.2
    rerun_diameter_calc = True
    print('------------------------------------------------------------------')
    print('run parameters:')
    print('molecule database file:', molecule_file)
    print('output dir:', output_dir)
    print('VDW scale:', vdwScale)
    print('Box Margin:', boxMargin, 'Angstrom')
    print('Grid spacing:', spacing, 'Angstrom')
    print('show VDW?:', show_vdw)
    print('Plot Ellipsoid?:', plot_ellip)
    print('No Conformers:', N_conformers)
    print('MW threshold:', MW_thresh, 'g/mol')
    print('pI threshold:', pI_thresh)
    print('Diffusion threshold:', size_thresh, 'Angstrom')
    print('Rerun diameter calculation?:', rerun_diameter_calc)
    print('------------------------------------------------------------------')

    df, molecules, diameters = rdkit_functions.read_mol_txt_file(molecule_file)
    # draw 2D structures
    print('--- draw 2D structures...')
    rdkit_functions.draw_svg_for_all_molecules(molecules,
                                               output_dir=output_dir)

    # calculate all Molecule Weights
    print('--- calculate MWs...')
    rdkit_functions.calculate_all_MW(molecules)

    # calculate the size of the ellipsoid surroudning all molecules
    print('--- calculate molecular diameters...')
    rdkit_functions.calc_molecule_diameters(
                    molecules, out_dir=output_dir, vdwScale=vdwScale,
                    boxMargin=boxMargin, spacing=spacing,
                    show_vdw=show_vdw, plot_ellip=plot_ellip,
                    N_conformers=N_conformers, MW_thresh=MW_thresh,
                    rerun=rerun_diameter_calc)

    # print results for each molecule
    print('--- print results and plot...')
    plotting.print_results(molecules,
                           threshold=size_thresh,
                           output_dir=output_dir)

    # plotting
    plotting.categorical(molecules,
                         threshold=size_thresh,
                         output_dir=output_dir)
    plotting.shapes(molecules,
                    threshold=size_thresh,
                    output_dir=output_dir)
