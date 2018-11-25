#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening of parameterisation molecule data set.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
from ercollect import rdkit_functions
import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ercollect import plotting
import pickle


def print_results_cf_known(molecules, known_df, threshold, output_dir):
    """Print the comparison of calculated and literature diffusion.

    """
    diffuse = {}
    no_diffuse = {}
    for name, smile in molecules.items():
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        print(out_file)
        if os.path.isfile(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        min_diam = min(results['diam1'])
        mid_diam = min(results['diam2'])
        lit_d = known_df[known_df['molecule'] == name]['diffuse'].iloc[0]
        if lit_d == 't':
            if mid_diam <= threshold:
                print(name+':')
                print('can diffuse')
                print('min diameter =', round(min_diam, 3), 'angstrom')
                print('mid diameter =', round(mid_diam, 3), 'angstrom')
                diffuse[name] = smile
            else:
                print(name+':')
                print('cannot diffuse - lit says it can!')
                print('min diameter =', round(min_diam, 3), 'angstrom')
                print('mid diameter =', round(mid_diam, 3), 'angstrom')
                no_diffuse[name] = smile
        else:
            if mid_diam <= threshold:
                print(name+':')
                print('can diffuse - lit says it cannot!')
                print('min diameter =', round(min_diam, 3), 'angstrom')
                print('mid diameter =', round(mid_diam, 3), 'angstrom')
                diffuse[name] = smile
            else:
                print(name+':')
                print('cannot diffuse')
                print('min diameter =', round(min_diam, 3), 'angstrom')
                print('mid diameter =', round(mid_diam, 3), 'angstrom')
                no_diffuse[name] = smile
        print('-')


def parity_with_known(molecules, diameters, known_df, threshold, output_dir):
    """Parity plot of calculated diameters and known kinetic diameters.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name, smile in molecules.items():
        try:
            kin_diam = float(diameters[name])
        except ValueError:
            print('no radius given for this molecule - skipped')
            continue
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        if os.path.isfile(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        C = 'none'
        M = 'o'
        E = 'k'
        print(name, kin_diam, mid_diam)
        ax.scatter(kin_diam, mid_diam, c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=80)

    # ax.axhspan(ymin=3.2, ymax=threshold, facecolor='k', alpha=0.2)
    # ax.axvspan(xmin=3.2, xmax=threshold, facecolor='k', alpha=0.2)
    ax.plot(np.linspace(-1, 12, 2), np.linspace(-1, 12, 2), c='k', alpha=0.4)
    # plot the limit from the two Sholl papers on diffusion
    # ax.axvspan(4.0, 4.2, facecolor='r', alpha=0.5)

    plotting.define_standard_plot(
                        ax,
                        title='',
                        xtitle='kinetic diameter [$\mathrm{\AA}$]',
                        ytitle='intermediate diameter [$\mathrm{\AA}$]',
                        xlim=(1, 10),
                        ylim=(1, 10))
    fig.tight_layout()
    fig.savefig(output_dir+"parity.pdf", dpi=720,
                bbox_inches='tight')


def categorical_with_known(molecules, known_df, threshold, output_dir):
    """Categorical scatter plot considering experimental results.

    """
    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))
    for name, smile in molecules.items():
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        if os.path.isfile(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        lit_d = known_df[known_df['molecule'] == name]['diffuse'].iloc[0]
        if lit_d == 't':
            if mid_diam <= threshold:
                C = 'b'
                M = 'o'
                E = 'k'
                D = 0.25
            else:
                C = 'b'
                M = 'o'
                E = 'k'
                D = 0.25
        elif lit_d == 'f':
            if mid_diam <= threshold:
                C = 'r'
                M = 'o'
                E = 'k'
                D = 0.75
            else:
                C = 'r'
                M = 'o'
                E = 'k'
                D = 0.75
        else:
            continue
        ax.scatter(D+(dx*(np.random.random() - 0.5) * 2),
                   mid_diam, c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=80)

    ax.axhspan(ymin=3.2, ymax=threshold, facecolor='k', alpha=0.2)

    plotting.define_diff_categ_plot(
                        ax,
                        title='',
                        xtitle='',
                        ytitle='intermediate diameter [$\mathrm{\AA}$]',
                        xlim=(0, 1),
                        ylim=(0, 10))
    fig.tight_layout()
    fig.savefig(output_dir+"categorical.pdf", dpi=720,
                bbox_inches='tight')


def parameter_tests(molecules, output_dir):
    """Function that scans over multiple parameters for a set of molecules
    and outputs figures.

    """
    inp = input('collect results already done? (T/F)')
    if inp == 'T':
        # read pickle
        rerun = False
        pass
    elif inp == 'F':
        rerun = True

    test_mol = [
                # 'carbon dioxide',
                'n-butane',
                # 'para-xylene',
                'meta-xylene',
                # 'n-hexane',
                'n-octane',
                # 'ethanol',
                'toluene'
                ]
    markers = {'carbon dioxide': 'o',
               'n-butane': 'X',
               'n-hexane': 'D',
               'n-octane': 'P',
               'ethanol': 'o',
               'para-xylene': '^',
               'meta-xylene': 'o',
               'toluene': '<'}
    colours = {'carbon dioxide': 'k',
               'n-butane': 'r',
               'n-hexane': 'purple',
               'n-octane': 'k',
               'ethanol': 'k',
               'para-xylene': 'b',
               'meta-xylene': 'b',
               'toluene': 'green'}

    values = {'space': [0.3, 0.4, 0.6, 0.8, 1.0],
              'conf': [10, 50, 100, 200],
              'vdw': [0.5, 0.8, 1.0],
              'box': [4, 6, 8]
              }
    test = ['space', 'conf', 'vdw', 'box']

    if rerun is True:
        full_results = {}
        for t in test:
            full_results[t] = {}
            for name, smile in molecules.items():
                if name in test_mol:
                    full_results[t][name] = []

        for name, smile in molecules.items():
            if name not in test_mol:
                continue
            for t in test:
                for v in values[t]:
                    # remove CSV files
                    os.system('rm '+output_dir+'*diam*.csv')
                    N_conformers = 50
                    boxMargin = 4
                    spacing = 0.6
                    vdwScale = 0.8
                    plot_ellip = False
                    show_vdw = False
                    MW_thresh = 2000
                    if t == 'conf':
                        N_conformers = v
                    if t == 'space':
                        spacing = v
                    if t == 'vdw':
                        vdwScale = v
                    if t == 'box':
                        boxMargin = v
                    print('doing', name, '- spacing:', spacing,
                          'vdw:', vdwScale, 'conf:', N_conformers,
                          'box:', boxMargin)
                    res = rdkit_functions.calc_molecule_diameter(
                                        name, smile,
                                        out_dir=output_dir,
                                        vdwScale=vdwScale,
                                        boxMargin=boxMargin,
                                        spacing=spacing,
                                        MW_thresh=MW_thresh,
                                        show_vdw=show_vdw,
                                        plot_ellip=plot_ellip,
                                        N_conformers=N_conformers)
                    min_diam_avg = np.average(res['diam1'])
                    min_diam_std = np.std(res['diam1'])
                    mid_diam_avg = np.average(res['diam2'])
                    mid_diam_std = np.std(res['diam2'])
                    min_mid = min(res['diam2'])
                    result = (min_diam_avg, min_diam_std,
                              mid_diam_avg, mid_diam_std,
                              min_mid)
                    full_results[t][name].append(result)
        # save file
        pickle.dump(full_results, open("param_test.pkl", "wb"))
    # load results
    full_results = pickle.load(open("param_test.pkl", "rb"))

    for t in test:
        fig, ax = plt.subplots()
        for name, smile in molecules.items():
            if name not in test_mol:
                continue
            X = []
            Y = []
            Y_err = []
            for i, v in enumerate(values[t]):
                min_diam_avg, min_diam_std, mid_diam_avg, mid_diam_std, min_mid = full_results[t][name][i]
                avg = float(min_diam_avg)
                std = float(min_diam_std)
                # if i == 0:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #                 yerr=std, fmt=markers[name], label=name)
                # else:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #                 yerr=std, fmt=markers[name])
                X.append(float(v))
                Y.append(avg)
                Y_err.append(std)
            if i == 0:
                ax.plot(X, Y, c=colours[name], marker=markers[name], label=name)
            else:
                ax.plot(X, Y, c=colours[name], marker=markers[name])
        if t == 'conf':
            t_lim = (0, 220)
            t_name = 'no. conformers'
        if t == 'space':
            t_lim = (0, 1.2)
            t_name = 'grid spacing [$\mathrm{\AA}$]'
        if t == 'vdw':
            t_lim = (0.4, 1.2)
            t_name = 'vdW scale parameter'
        if t == 'box':
            t_lim = (2, 10)
            t_name = 'box margin [$\mathrm{\AA}$]'
        plotting.define_standard_plot(
                            ax,
                            title='',
                            xtitle=t_name,
                            ytitle='avg. minimum diameter [$\mathrm{\AA}$]',
                            xlim=t_lim,
                            ylim=(0, 10))
        ax.legend(loc=1, fontsize=16)
        fig.tight_layout()
        fig.savefig(output_dir+"min_"+t+".pdf", dpi=720,
                    bbox_inches='tight')
    for t in test:
        fig, ax = plt.subplots()
        for name, smile in molecules.items():
            if name not in test_mol:
                continue
            X = []
            Y = []
            Y_err = []
            for i, v in enumerate(values[t]):
                min_diam_avg, min_diam_std, mid_diam_avg, mid_diam_std, min_mid = full_results[t][name][i]
                avg = float(mid_diam_avg)
                std = float(mid_diam_std)
                # if i == 0:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #                 yerr=std, fmt=markers[name], label=name)
                # else:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #                 yerr=std, fmt=markers[name])
                X.append(float(v))
                Y.append(avg)
                Y_err.append(std)
            X = np.asarray(X)
            Y = np.asarray(Y)
            Y_err = np.asarray(Y_err)
            ax.plot(X, Y, c=colours[name], marker=markers[name],
                    label=name)
            ax.fill_between(X, Y-Y_err, Y+Y_err, alpha=0.2,
                            facecolor=colours[name])
        if t == 'conf':
            t_lim = (0, 220)
            t_name = 'no. conformers'
        if t == 'space':
            t_lim = (0.2, 1.1)
            t_name = 'grid spacing [$\mathrm{\AA}$]'
        if t == 'vdw':
            t_lim = (0.4, 1.1)
            t_name = 'vdW scale parameter'
        if t == 'box':
            t_lim = (3, 9)
            t_name = 'box margin [$\mathrm{\AA}$]'
        plotting.define_standard_plot(
                            ax,
                            title='',
                            xtitle=t_name,
                            ytitle='avg. intermediate diameter [$\mathrm{\AA}$]',
                            xlim=t_lim,
                            ylim=(3.5, 9))
        # ax.legend(fontsize=16, ncol=2)
        fig.tight_layout()
        fig.savefig(output_dir+"mid_"+t+".pdf", dpi=720,
                    bbox_inches='tight')
    for t in test:
        fig, ax = plt.subplots()
        for name, smile in molecules.items():
            if name not in test_mol:
                continue
            X = []
            Y = []
            Y_err = []
            for i, v in enumerate(values[t]):
                min_diam_avg, min_diam_std, mid_diam_avg, mid_diam_std, min_mid = full_results[t][name][i]
                # if i == 0:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #                 yerr=std, fmt=markers[name], label=name)
                # else:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #                 yerr=std, fmt=markers[name])
                X.append(float(v))
                Y.append(min_mid)
            X = np.asarray(X)
            Y = np.asarray(Y)
            Y_err = np.asarray(Y_err)
            ax.plot(X, Y, c=colours[name], marker=markers[name],
                    label=name)
        if t == 'conf':
            t_lim = (0, 220)
            t_name = 'no. conformers'
        if t == 'space':
            t_lim = (0.2, 1.1)
            t_name = 'grid spacing [$\mathrm{\AA}$]'
        if t == 'vdw':
            t_lim = (0.4, 1.1)
            t_name = 'vdW scale parameter'
        if t == 'box':
            t_lim = (3, 9)
            t_name = 'box margin [$\mathrm{\AA}$]'
        plotting.define_standard_plot(
                            ax,
                            title='',
                            xtitle=t_name,
                            ytitle='min. intermediate diameter [$\mathrm{\AA}$]',
                            xlim=t_lim,
                            ylim=(3.5, 9))
        # ax.legend(fontsize=16, ncol=2)
        fig.tight_layout()
        fig.savefig(output_dir+"min_of_mid_"+t+".pdf", dpi=720,
                    bbox_inches='tight')


def shapes_with_known(molecules, known_df, threshold, output_dir):
    """Plot molecule shapes considering experimental results.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name, smile in molecules.items():
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        if os.path.isfile(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        lit_d = known_df[known_df['molecule'] == name]['diffuse'].iloc[0]
        if lit_d == 't':
            if mid_diam <= threshold:
                C = 'b'
                M = 'o'
                E = 'k'
            else:
                C = 'b'
                M = 'X'
                E = 'k'
        elif lit_d == 'f':
            if mid_diam <= threshold:
                C = 'r'
                M = 'X'
                E = 'k'
                D = 0.75
            else:
                C = 'r'
                M = 'o'
                E = 'k'
                D = 0.75
        else:
            continue
        ax.scatter(np.average(results['ratio_1']),
                   np.average(results['ratio_2']),
                   c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=80)

    ax.plot([0, 0.5, 1, 0], [1, 0.5, 1, 1], c='k', lw=2)
    ax.text(0.75, 1.03, 'sphere', fontsize=20)
    ax.text(0.4, 0.45, 'oblate', fontsize=20)
    ax.text(-0.05, 1.03, 'prolate', fontsize=20)

    plotting.define_standard_plot(
                        ax,
                        title='',
                        xtitle='$I_1$ / $I_3$',
                        ytitle='$I_2$ / $I_3$',
                        xlim=(-0.1, 1.1),
                        ylim=(0.4, 1.1))
    fig.tight_layout()
    fig.savefig(output_dir+"shape.pdf", dpi=720,
                bbox_inches='tight')


if __name__ == "__main__":
    start = time.time()
    # set parameters
    # molecule file dir
    molecule_file = '/home/atarzia/psp/molecule_param/test_molecules.txt'
    # output dir
    output_dir = '/home/atarzia/psp/molecule_param/'
    vdwScale = 0.8
    boxMargin = 4.0
    spacing = 0.6
    show_vdw = False
    plot_ellip = False
    N_conformers = 50
    MW_thresh = 2000
    pI_thresh = 6
    size_thresh = 4.2
    rerun_diameter_calc = False
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
    print('Diffusion threshold:', size_thresh, 'Ansgtrom')
    print('Rerun diameter calculation?:', rerun_diameter_calc)
    print('------------------------------------------------------------------')

    print('------------------------------------------------------------------')
    print('Screen molecular size of compounds in known reactions')
    print('------------------------------------------------------------------')

    df, molecules, diameters = rdkit_functions.read_mol_txt_file(molecule_file)
    # draw 2D structures
    print('--- draw 2D structures...')
    rdkit_functions.draw_svg_for_all_molecules(molecules,
                                               output_dir=output_dir)

    # calculate all Molecule Weights
    print('--- calculate MWs...')
    rdkit_functions.calculate_all_MW(molecules)

    input('happy with MW cut-off?')
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
    print_results_cf_known(molecules,
                           known_df=df,
                           threshold=size_thresh,
                           output_dir=output_dir)

    # plotting
    parity_with_known(molecules, diameters,
                      known_df=df,
                      threshold=size_thresh,
                      output_dir=output_dir)
    categorical_with_known(molecules,
                           known_df=df,
                           threshold=size_thresh,
                           output_dir=output_dir)
    shapes_with_known(molecules,
                      known_df=df,
                      threshold=size_thresh,
                      output_dir=output_dir)
    parameter_tests(molecules,
                    output_dir=output_dir)
    end = time.time()
    print('---- total time taken =', '{0:.2f}'.format(end-start), 's')
