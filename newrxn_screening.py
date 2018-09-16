#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for new reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import pandas as pd
import glob
import time
import pi_fn
import rdkit_functions
import plotting
import DB_functions
import rxn_syst

# script/data set specific functions


def get_molecule_DB(EC_mol_set, output_dir):
    """Get molecule dictionary + output 2D structures.

    """
    molecules = {}
    diameters = {}
    for i in EC_mol_set.keys():
        for j in EC_mol_set[i].keys():
            for mol in EC_mol_set[i][j]:
                molecules[mol[0]] = mol[1]
                diameters[mol[0]] = 0
                rdkit_functions.draw_smiles_to_svg(
                    mol[1], output_dir+mol[0].replace(' ', '_')+'_2d.svg')
    return molecules, diameters


if __name__ == "__main__":
    start = time.time()
    # set parameters
    EC_set, EC_mol_set, EC_descriptors = EC_sets()
    # pI
    pI_DB_dir = '/home/atarzia/psp/sequence_db/bio_min_dataset/'
    pI_output_dir = pI_DB_dir
    pI_csv = "output_data_pi.csv"
    redo_pI = False
    redo_pI_plots = False
    pI_thresh = 6
    # known molecule screening
    mol_DB_dir = '/home/atarzia/psp/screening_results/biomin_known/'
    mol_output_dir = mol_DB_dir
    vdwScale = 0.8
    boxMargin = 4.0
    spacing = 0.6
    show_vdw = False
    plot_ellip = False
    N_conformers = 50
    MW_thresh = 2000
    size_thresh = 4.2
    rerun_diameter_calc = False
    # reaction search
    search_DBs = ['BRENDA', 'SABIO', 'KEGG', 'BKMS', ]
    search_output_dir = '/home/atarzia/psp/screening_results/biomin_search/'
    search_ECs = ['1.11.1.5', '1.11.1.6', '1.11.1.7', '1.9.3.1',
                  '1.1.5.2', '3.5.1.5', '1.1.3.4', '1.13.12.4',
                  '3.2.1.26', '3.1.1.3', '3.1.1.6', '3.5.1.11']
    search_mol_output_file = search_output_dir+'screening_output.csv'
    search_MW_thresh = 250
    print('------------------------------------------------------------------')
    print('run parameters:')
    print('pI database dir:', pI_DB_dir)
    print('pI output dir:', pI_output_dir)
    print('Redo pI screening?:', redo_pI)
    print('Redo pI plotting?:', redo_pI_plots)
    print('molecule database dir:', mol_DB_dir)
    print('molecule output dir:', mol_output_dir)
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
    print('Search output dir:', search_output_dir)
    print('Search molecule output file:', search_mol_output_file)
    print('Search MW threshold:', search_MW_thresh, 'g/mol')
    print('------------------------------------------------------------------')

    print('------------------------------------------------------------------')
    print('Screen pIs')
    print('------------------------------------------------------------------')
    temp_time = time.time()
    # prepare pI calculations
    database_names = pi_fn.prepare_pI_calc(database_directory=pI_DB_dir,
                                           redo_pi=redo_pI,
                                           output_dir=pI_output_dir,
                                           csv=pI_csv)
    # screen protein sequence from EC numbers
    print('--- calculate all pIs for target EC sequences...')
    pi_fn.screen_pIs(database_names, redo_pI_plots=redo_pI_plots,
                     redo_pI=redo_pI, pI_csv=pI_csv,
                     pI_output_dir=pI_output_dir, cutoff_pi=pI_thresh,
                     descriptors=EC_descriptors)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    print('------------------------------------------------------------------')
    print('Screen new reactions')
    print('------------------------------------------------------------------')
    temp_time = time.time()
    print('collect all reaction systems (ONLINE)...')
    for DB in search_DBs:
        # get database specific information
        DB_prop = DB_functions.get_DB_prop(DB)
        db_dir = DB_prop[0]
        # iterate over EC numbers of interest
        for EC in search_ECs:
            rxn_syst.get_reaction_systems(EC, DB,
                                          search_output_dir,
                                          clean_system=False)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    rxn_syst.percent_skipped(search_output_dir)
    print('check all reaction systems for diffusion of components (ONLINE)...')
    rxn_syst.check_all_RS_diffusion(output_dir=search_output_dir,
                                    mol_output_file=search_mol_output_file,
                                    threshold=size_thresh,
                                    vdwScale=vdwScale,
                                    boxMargin=boxMargin, spacing=spacing,
                                    N_conformers=N_conformers,
                                    MW_thresh=search_MW_thresh)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    rxn_syst.percent_skipped(search_output_dir)
    print('get subset of reactions with known protein sequences...')
    rxn_syst.check_all_seedMOF(search_output_dir, pI_thresh)
    rxn_syst.percent_w_sequence(search_output_dir)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    print('--- print results and plot...')
    # plot distribution of pI of all known sequences
    plotting.rs_pI_distribution(output_dir=search_output_dir,
                                cutoff_pI=pI_thresh)
    # plot max component size vs pI
    plotting.rs_size_vs_pI(output_dir=search_output_dir,
                           cutoff_pI=pI_thresh,
                           size_thresh=size_thresh)
    # plot number of new reactions as a function of size threshold
    plotting.number_rxns_vs_size(output_dir=search_output_dir,
                                 size_thresh=size_thresh)
    # categorize all molecules in mol output file
    plotting.categorical_moloutput(mol_output_file=search_mol_output_file,
                                   threshold=size_thresh,
                                   output_dir=search_output_dir)
    # print new reactions
    plotting.print_new_rxns(output_dir=search_output_dir)

    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    end = time.time()
    print('---- total time taken =', '{0:.2f}'.format(end-start), 's')
