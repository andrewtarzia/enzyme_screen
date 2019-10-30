#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for new reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import glob
import pandas as pd
import time
from ercollect import pi_fn
from ercollect import plotting
from ercollect import rxn_syst
from multiprocessing import Pool


# script/data set specific functions
def get_ECs_from_file(EC_file):
    """Read in ECs to search from a file.

    """
    # get search EC numbers from file:
    # set EC numbers of interest
    # get from a data file - manually made from
    # https://enzyme.expasy.org/enzyme-byclass.html
    EC_DF = pd.read_table(
        EC_file,
        delimiter='__',
        names=['EC_no', 'description'],
        engine='python'
    )
    search_ECs = list(EC_DF['EC_no'])

    print(len(search_ECs), 'EC numbers to test')
    print('first EC:', search_ECs[0], '---- last EC:', search_ECs[-1])
    print('collect all reaction systems (ONLINE)...')
    return search_ECs


if __name__ == "__main__":
    start = time.time()
    # set parameters
    # pI
    pI_DB_dir = (
        '/home/atarzia/psp/screening_results/new_reactions/sequences/'
    )
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
    search_DBs = ['KEGG', 'BKMS', 'BRENDA', 'SABIO', ]
    search_output_dir = (
        '/home/atarzia/psp/screening_results/new_reactions/'
    )
    search_EC_file = search_output_dir+'desired_EC.txt'
    search_mol_output_file = search_output_dir+'screening_output.csv'
    search_MW_thresh = 250
    search_run = True
    search_redo = False
    collect_mol_prop = True
    NP = 2  # number of processes
    print('----------------------------------------------------------')
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
    print('Search output dir:', search_output_dir)
    print('EC file:', search_EC_file)
    print('Search molecule output file:', search_mol_output_file)
    print('Search MW threshold:', search_MW_thresh, 'g/mol')
    print('Run Search?:', search_run)
    print('Redo Search?:', search_redo)
    print('Collect molecule properties?:', collect_mol_prop)
    print('----------------------------------------------------------')

    print('----------------------------------------------------------')
    print('Screen pIs')
    print('----------------------------------------------------------')
    temp_time = time.time()
    # prepare pI calculations
    database_names = pi_fn.prepare_pI_calc(
        database_directory=pI_DB_dir,
        redo_pi=redo_pI,
        output_dir=pI_output_dir,
        csv=pI_csv
    )
    # screen protein sequence from EC numbers
    print('--- calculate all pIs for target EC sequences...')
    pi_fn.screen_pIs(database_names, redo_pI_plots=redo_pI_plots,
                     redo_pI=redo_pI, pI_csv=pI_csv,
                     pI_output_dir=pI_output_dir, cutoff_pi=pI_thresh,
                     descriptors=None)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    print('----------------------------------------------------------')
    print('Screen new reactions')
    print('----------------------------------------------------------')
    temp_time = time.time()
    search_ECs = get_ECs_from_file(EC_file=search_EC_file)
    if search_run is True:
        for DB in search_DBs:
            # iterate over EC numbers of interest
            # Create a multiprocessing Pool
            with Pool(NP) as pool:
                # process data_inputs iterable with pool
                # func(EC, DB, search_output_dir, search_redo, verbose)
                args = [(EC, DB, search_output_dir, search_redo, True)
                        for EC in search_ECs]
                pool.starmap(rxn_syst.process_collection, args)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    rxn_syst.percent_skipped(output_dir=search_output_dir)
    react_syst_files = glob.glob(search_output_dir+'sRS-*.pkl')
    if collect_mol_prop is True:
        print('collect all molecule properties (ONLINE)...')
        print('if check is False - overwrite previous settings.')
        # iterate over reaction systems
        # Create a multiprocessing Pool
        with Pool(NP) as pool:
            # process data_inputs iterable with pool
            # func(EC, DB, search_output_dir, search_redo, verbose)
            args = [
                (rs, search_output_dir, False, i, react_syst_files)
                for i, rs in enumerate(
                    rxn_syst.yield_rxn_syst(search_output_dir)
                )
            ]
            pool.starmap(rxn_syst.process_molecule_collection, args)
        # rxn_syst.collect_all_molecule_properties(output_dir=search_output_dir,
        #                                          check=False)
    print(
        'check all reaction systems for diffusion of components '
        '(ONLINE)...'
    )
    # iterate over reaction systems
    # Create a multiprocessing Pool
    with Pool(NP) as pool:
        # process data_inputs iterable with pool
        # func(EC, DB, search_output_dir, search_redo, verbose)
        args = [
            (
                rs,
                i,
                react_syst_files,
                search_output_dir,
                search_mol_output_file,
                size_thresh,
                vdwScale,
                boxMargin,
                spacing,
                N_conformers,
                search_MW_thresh
            )
            for i, rs in enumerate(
                rxn_syst.yield_rxn_syst(search_output_dir)
            )
        ]
        pool.starmap(rxn_syst.process_RS_diffusion, args)
    # rxn_syst.check_all_RS_diffusion(
    #     output_dir=search_output_dir,
    #     mol_output_file=search_mol_output_file,
    #     threshold=size_thresh,
    #     vdwScale=vdwScale,
    #     boxMargin=boxMargin, spacing=spacing,
    #     N_conformers=N_conformers,
    #     MW_thresh=search_MW_thresh
    # )
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    rxn_syst.percent_skipped(search_output_dir)
    print('get subset of reactions with known protein sequences...')
    rxn_syst.check_all_seedMOF(search_output_dir, pI_thresh)
    rxn_syst.percent_w_sequence(search_output_dir)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    print('determine solubility range of all reactions using logP...')
    rxn_syst.check_all_solubility(output_dir=search_output_dir)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    print('determine change in synthetic accessibility all reactions.')
    rxn_syst.delta_sa_score(output_dir=search_output_dir)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    print('determine solubility range of all reactions using XlogP...')
    rxn_syst.check_all_solubility_X(output_dir=search_output_dir)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    print('determine change in molecular complexity all reactions...')
    rxn_syst.delta_complexity_score(output_dir=search_output_dir)
    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    print('--- print results and plot...')
    # plot a distribution of the number of reactnts in each RS
    plotting.rs_no_reactants(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot a distribution of the number of products in each RS
    plotting.rs_no_products(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot a distribution of the change in synthetic accesibility
    plotting.rs_dist_deltaSA(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot a distribution of all molecule complexity
    plotting.rs_dist_complexity(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot a distribution of the change in complexity
    plotting.rs_dist_deltacomplexity(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot max component size vs synthetic accessibility vs logP
    plotting.rs_size_vs_SA_vs_logP(
        output_dir=search_output_dir,
        size_thresh=size_thresh,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot max component size vs complexity vs XlogP
    plotting.rs_size_vs_complexity_vs_XlogP(
        output_dir=search_output_dir,
        size_thresh=size_thresh,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot number of new reactions as a function of size threshold
    plotting.rs_number_rxns_vs_size(
        output_dir=search_output_dir,
        size_thresh=size_thresh,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot distribution of pI of all known sequences
    plotting.rs_pI_distribution(
        output_dir=search_output_dir,
        cutoff_pI=pI_thresh,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot max component size vs pI
    plotting.rs_size_vs_pI(
        output_dir=search_output_dir,
        cutoff_pI=pI_thresh,
        size_thresh=size_thresh,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # categorize all molecules in mol output file
    plotting.categorical_moloutput(
        mol_output_file=search_mol_output_file,
        threshold=size_thresh,
        output_dir=search_output_dir
    )
    # print new reactions
    plotting.print_new_rxns(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )
    # plot a distribution of the change in molecule size due to
    # reaction
    plotting.rs_delta_size(
        output_dir=search_output_dir,
        generator=rxn_syst.yield_rxn_syst(search_output_dir)
    )

    print(
        '---- step time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )
    temp_time = time.time()
    end = time.time()
    print('---- total time taken =', '{0:.2f}'.format(end-start), 's')
