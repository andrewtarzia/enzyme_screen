#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for new reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import pandas as pd
import time
import pi_fn
import plotting
import DB_functions
import rxn_syst

# script/data set specific functions


def get_ECs_from_file(EC_file):
    """Read in ECs to search from a file.

    """
    # get search EC numbers from file:
    # set EC numbers of interest
    # get from a data file - manually made from
    # https://enzyme.expasy.org/enzyme-byclass.html
    EC_DF = pd.read_table(EC_file, delimiter='__',
                          names=['EC_no', 'description'], engine='python')
    search_ECs = list(EC_DF['EC_no'])

    print(len(search_ECs), 'EC numbers to test')
    print('first EC:', search_ECs[0], '---- last EC:', search_ECs[-1])
    print('collect all reaction systems (ONLINE)...')
    return search_ECs


if __name__ == "__main__":
    start = time.time()
    # set parameters
    # pI
    pI_DB_dir = '/home/atarzia/psp/screening_results/new_reactions/sequences/'
    pI_output_dir = pI_DB_dir
    pI_csv = "output_data_pi.csv"
    redo_pI = False
    redo_pI_plots = False
    pI_thresh = 6
    vdwScale = 0.8
    boxMargin = 4.0
    spacing = 0.6
    show_vdw = False
    plot_ellip = False
    N_conformers = 50
    size_thresh = 4.2
    # reaction search
    search_DBs = ['BRENDA', 'SABIO', 'KEGG', 'BKMS', ]
    search_output_dir = '/home/atarzia/psp/screening_results/new_reactions/'
    search_EC_file = search_output_dir+'desired_EC.txt'
    search_mol_output_file = search_output_dir+'screening_output.csv'
    search_MW_thresh = 250
    search_redo = False
    print('------------------------------------------------------------------')
    print('run parameters:')
    print('pI database dir:', pI_DB_dir)
    print('pI output dir:', pI_output_dir)
    print('Redo pI screening?:', redo_pI)
    print('Redo pI plotting?:', redo_pI_plots)
    print('VDW scale:', vdwScale)
    print('Box Margin:', boxMargin, 'Angstrom')
    print('Grid spacing:', spacing, 'Angstrom')
    print('show VDW?:', show_vdw)
    print('Plot Ellipsoid?:', plot_ellip)
    print('No Conformers:', N_conformers)
    print('pI threshold:', pI_thresh)
    print('Diffusion threshold:', size_thresh, 'Angstrom')
    print('Rerun diameter calculation?:', rerun_diameter_calc)
    print('Search output dir:', search_output_dir)
    print('EC file:', search_EC_file)
    print('Search molecule output file:', search_mol_output_file)
    print('Search MW threshold:', search_MW_thresh, 'g/mol')
    print('Redo Search?:', search_redo)
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
                     descriptors=None)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    print('------------------------------------------------------------------')
    print('Screen new reactions')
    print('------------------------------------------------------------------')
    temp_time = time.time()
    search_ECs = get_ECs_from_file(EC_file=search_EC_file)
    for DB in search_DBs:
        # get database specific information
        DB_prop = DB_functions.get_DB_prop(DB)
        db_dir = DB_prop[0]
        # iterate over EC numbers of interest
        for EC in search_ECs:
            print('doing:', DB, 'EC:', EC)
            rxn_syst.get_reaction_systems(EC, DB,
                                          search_output_dir,
                                          clean_system=search_redo,
                                          verbose=True)
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
