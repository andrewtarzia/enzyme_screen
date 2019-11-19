#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting functions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import os
import sys

from reaction import yield_rxn_syst
import plotting_fn as pfn


if __name__ == "__main__":
    if (not len(sys.argv) == 2):
        print('Usage: plotting.py plot_suffix\n')
        print(
            'plot_suffix: string to put at the end of plot file names.'
        )
        sys.exit()
    else:
        plot_suffix = sys.argv[1]

    search_output_dir = os.getcwd()+'/'
    pI_thresh = 6
    size_thresh = 4.2  # angstroms
    print('settings:')
    print('    pI threshold:', pI_thresh)
    print('    Diffusion threshold:', size_thresh, 'Angstrom')
    DB_switch = input('biomin (1) or new (2) or KEGG/ATLAS (3)?')
    if DB_switch == '1':
        DB_switch = 1
    elif DB_switch == '2':
        DB_switch = 2
    elif DB_switch == '3':
        DB_switch = 3
    else:
        print('answer correctly...')
        sys.exit()

    print('------------------------------------------------------')
    print('Plotting all the plots')
    print('------------------------------------------------------')
    #######
    # RS property plots
    #######
    # plot number of new reactions as a function of size threshold
    if input("do % skipped? (t/f)") == 't':
        print('doing...')
        from ercollect.rxn_syst import percent_skipped
        percent_skipped(output_dir=search_output_dir)
    if input('do no. rxns w size? (t/f)') == 't':
        print('doing....')
        pfn.rs_number_rxns_vs_size(output_dir=search_output_dir,
                               size_thresh=size_thresh,
                               generator=yield_rxn_syst(search_output_dir),
                               plot_suffix=plot_suffix)
    if input('do dist_logPs? (t/f)') == 't':
        print('doing....')
        pfn.rs_dist_logP(output_dir=search_output_dir,
                     generator=yield_rxn_syst(search_output_dir),
                     plot_suffix=plot_suffix,
                     extreme='max')
        # rs_dist_logP(output_dir=search_output_dir,
        #              generator=yield_rxn_syst(search_output_dir),
        #              plot_suffix=plot_suffix,
        #              extreme='min')
    if input('do dist_logS? (t/f)') == 't':
        print('doing....')
        # rs_dist_logS(output_dir=search_output_dir,
        #              generator=yield_rxn_syst(search_output_dir),
        #              plot_suffix=plot_suffix,
        #              extreme='max')
        pfn.rs_dist_logS(output_dir=search_output_dir,
                     generator=yield_rxn_syst(search_output_dir),
                     plot_suffix=plot_suffix,
                     extreme='min')
    if input('do dist_logPs -- no EC? (t/f)') == 't':
        print('doing....')
        pfn.rs_dist_logP_noEC(output_dir=search_output_dir,
                          generator=yield_rxn_syst(search_output_dir),
                          plot_suffix=plot_suffix,
                          extreme='max')
        # rs_dist_logP(output_dir=search_output_dir,
        #              generator=yield_rxn_syst(search_output_dir),
        #              plot_suffix=plot_suffix,
        #              extreme='min')
    if input('do dist_logS -- no EC? (t/f)') == 't':
        print('doing....')
        # rs_dist_logS(output_dir=search_output_dir,
        #              generator=yield_rxn_syst(search_output_dir),
        #              plot_suffix=plot_suffix,
        #              extreme='max')
        pfn.rs_dist_logS_noEC(output_dir=search_output_dir,
                          generator=yield_rxn_syst(search_output_dir),
                          plot_suffix=plot_suffix,
                          extreme='min')
    # rs_dist_delta_complexity_vs_size(
    #                 output_dir=search_output_dir,
    #                 generator=yield_rxn_syst(search_output_dir),
    #                 plot_suffix=plot_suffix)
    if DB_switch == 1:
        # print new reactions
        pfn.print_new_rxns(output_dir=search_output_dir,
                       generator=yield_rxn_syst(search_output_dir))
    if input('do dist_no prod and reacts? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the number of reactnts in each reaction system
        pfn.rs_dist_no_reactants(output_dir=search_output_dir,
                             generator=yield_rxn_syst(search_output_dir),
                             plot_suffix=plot_suffix)
        # plot a distribution of the number of products in each reaction system
        pfn.rs_dist_no_products(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)
    if input('do dist max SIZE? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        pfn.rs_dist_max_size(output_dir=search_output_dir,
                         generator=yield_rxn_syst(search_output_dir),
                         plot_suffix=plot_suffix)
    if input('do dist max SIZE -- no EC? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        pfn.rs_dist_max_size_noEC(output_dir=search_output_dir,
                              generator=yield_rxn_syst(search_output_dir),
                              plot_suffix=plot_suffix)
    if input('do violin max SIZE? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        pfn.violin_max_size(output_dir=search_output_dir,
                        generator=yield_rxn_syst(search_output_dir),
                        plot_suffix=plot_suffix)
    if input('do dist_delta SIZE? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        pfn.rs_dist_delta_size(output_dir=search_output_dir,
                           generator=yield_rxn_syst(search_output_dir),
                           plot_suffix=plot_suffix)
    if input('do dist_deltaSA? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in synthetic accesibility
        pfn.rs_dist_delta_SA(output_dir=search_output_dir,
                         generator=yield_rxn_syst(search_output_dir),
                         plot_suffix=plot_suffix)
    if input('do dist_deltaSA? -- no EC (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in synthetic accesibility
        pfn.rs_dist_delta_SA_noEC(output_dir=search_output_dir,
                              generator=yield_rxn_syst(search_output_dir),
                              plot_suffix=plot_suffix)

    if input('do dist_deltacomplexity? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in synthetic accesibility
        pfn.rs_dist_delta_complexity(output_dir=search_output_dir,
                                 generator=yield_rxn_syst(search_output_dir),
                                 plot_suffix=plot_suffix)

    # plot distributions of protein sequence properties
    if DB_switch != 3:
        if input('do dist_GRAVY? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_GRAVY(output_dir=search_output_dir,
                          generator=yield_rxn_syst(search_output_dir),
                          plot_suffix=plot_suffix)
        if input('do dist_I index? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_I_index(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)
        if input('do dist_A index? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_A_index(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)

        if input('do dist_TM index? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_TM_index(output_dir=search_output_dir,
                             generator=yield_rxn_syst(search_output_dir),
                             plot_suffix=plot_suffix)
        if input('do dist_pI? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_pI(output_dir=search_output_dir,
                       generator=yield_rxn_syst(search_output_dir),
                       plot_suffix=plot_suffix)
    sys.exit()
    # plot max component size vs synthetic accessibility vs logP
    # rs_size_vs_SA_vs_logP(output_dir=search_output_dir,
    #                       size_thresh=size_thresh,
    #                       generator=yield_rxn_syst(search_output_dir),
    #                       plot_suffix=plot_suffix)
    # # plot max component size vs complexity vs XlogP
    # rs_size_vs_complexity_vs_XlogP(output_dir=search_output_dir,
    #                                size_thresh=size_thresh,
    #                                generator=yield_rxn_syst(search_output_dir),
    #                                plot_suffix=plot_suffix)
    # plot max component size vs SA score vs XlogP vs aliphatic index
    # rs_size_vs_SA_vs_XlogP_vs_aindex(output_dir=search_output_dir,
    #                                  size_thresh=size_thresh,
    #                                  generator=yield_rxn_syst(search_output_dir),
    #                                  plot_suffix=plot_suffix)
