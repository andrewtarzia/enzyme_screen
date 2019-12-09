#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening of parameterisation molecule data set.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error
import pickle
import sys

import utilities
import rdkit_functions as rdkf
import plotting_fn as pfn


def print_results_cf_known(molecules, known_df, threshold, output_dir):
    """
    Print the comparison of calculated and literature diffusion.

    """
    diffuse = {}
    no_diffuse = {}
    for name in molecules:
        smile = molecules[name]
        out_file = (
            f"{output_dir}/{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        min_diam = min(results['diam1'])
        mid_diam = min(results['diam2'])
        lit_d = known_df[
            known_df['molecule'] == name
        ]['diffuse'].iloc[0]
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


def parity_with_known(
    molecules,
    diameters,
    output_dir
):
    """
    Parity plot of calculated diameters and known kinetic diameters.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name in molecules:
        try:
            kin_diam = float(diameters[name])
        except ValueError:
            print('no radius given for this molecule - skipped')
            continue
        out_file = (
            f"{output_dir}/{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        C = 'none'
        M = 'o'
        print(name, kin_diam, mid_diam)
        ax.scatter(kin_diam, mid_diam, c=C,
                   edgecolors='k', marker=M, alpha=1.0,
                   s=80)

    ax.plot(
        np.linspace(-1, 12, 2),
        np.linspace(-1, 12, 2),
        c='k',
        alpha=0.4
    )
    # plot the limit from the two Sholl papers on diffusion
    # ax.axvspan(4.0, 4.2, facecolor='r', alpha=0.5)

    pfn.define_standard_plot(
        ax,
        xtitle=r'kinetic diameter [$\mathrm{\AA}$]',
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xlim=(1, 10),
        ylim=(1, 10)
    )
    fig.tight_layout()
    fig.savefig("parity.pdf", dpi=720, bbox_inches='tight')


def parity_cf_scale_with_known(
    molecules,
    diameters,
    known_df,
    pars,
    scale_info
):
    """
    Produce a parity plot of calculated diameters and known kinetic
    diameters for multiple input parameters.

    """

    fig, ax = plt.subplots(figsize=(5, 5))
    for dir in scale_info:
        kin_diams = []
        mid_diams = []
        sc, C, M, A, E = scale_info[dir]
        scale_output = f'scale_sc_{dir}.txt'
        if os.path.exists(scale_output):
            with open(scale_output, 'r') as f:
                for line in f:
                    res = line.rstrip().split('__')
                    name, kin_diam, mid_diam = res
                    kin_diams.append(float(kin_diam))
                    mid_diams.append(float(mid_diam))
                    ax.scatter(
                        float(kin_diam),
                        float(mid_diam),
                        c=C,
                        edgecolors=E,
                        marker=M,
                        alpha=A,
                        s=40
                    )
        else:
            with open(scale_output, 'w') as f:
                for name in molecules:
                    try:
                        kin_diam = float(diameters[name])
                    except ValueError:
                        print(
                            'no radius given for this molecule '
                            '- skipped'
                        )
                        continue
                    out_file = (
                        f"{dir}/"
                        f"{name.replace(' ', '_').replace('/', '__')}"
                        '_diam_result.csv'
                    )
                    if os.path.exists(out_file) is False:
                        continue
                    results = pd.read_csv(out_file)
                    if len(results) == 0:
                        continue
                    mid_diam = min(results['diam2'])
                    kin_diams.append(float(kin_diam))
                    mid_diams.append(float(mid_diam))
                    ax.scatter(
                        float(kin_diam),
                        float(mid_diam),
                        c=C,
                        edgecolors=E,
                        marker=M,
                        alpha=A,
                        s=40
                    )
                    f.write(
                        name+'__'+str(kin_diam)+'__'+str(mid_diam)+'\n'
                    )
        corr = pearsonr(kin_diams, mid_diams)
        MAE = mean_absolute_error(kin_diams, mid_diams)
        print(f'{dir} R^2: {corr}, MAE: {MAE}')

    ax.plot(
        np.linspace(-1, 12, 2),
        np.linspace(-1, 12, 2),
        c='k',
        alpha=0.4
    )

    pfn.define_standard_plot(
        ax,
        xtitle=r'kinetic diameter [$\mathrm{\AA}$]',
        # ytitle='intermediate diameter [$\mathrm{\AA}$]',
        ytitle=r'$d$ [$\mathrm{\AA}$]',
        xlim=(1, 10),
        ylim=(1, 10)
    )

    # legend
    for dir in scale_info:
        sc, C, M, A, E = scale_info[dir]
        ax.scatter(
            -100, -100,
            c=C,
            edgecolors=E,
            marker=M,
            alpha=A,
            s=40,
            label=f'vdW scale = {sc}'
        )
    ax.legend(loc=2, fontsize=14)
    fig.tight_layout()
    fig.savefig(
        "parity_scalecf.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def cf_verploegh2015(molecules, output_dir):
    """
    Recreate Figure 4 in verploegh2015 (D_self at 35 degC).

    D_self data hardcoded from supp info of that paper.

    """
    D_self = {
        'He': 1.61E-04,
        'H2': 1.80E-04,
        'oxygen': 9.41E-06,
        'nitrogen': 1.18E-06,
        'carbon dioxide': 2.63E-06,
        'methane': 2.79E-07,
        'SF6': 3.73E-17,
        'ethene': 4.69E-08,  # C2H4=
        'ethane': 2.36E-08,  # C2H6
        'propene': 2.58E-09,  # C3H6=
        'n-propane': 1.38E-10,  # C3H8
        '1-butene': 1.48E-10,  # 1-C4H8=
        'n-butane': 9.53E-11,  # n-C4H10
        'i-butene': 5.92E-15,  # iso-C4H8=
        'i-butane': 1.36E-16,  # iso-C4H10
    }

    fig, ax = plt.subplots()
    for name in molecules:
        out_file = (
            f"{output_dir}/{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        if name in D_self:
            DS = D_self[name]
            ax.scatter(
                mid_diam,
                DS,
                c='#FF7900',
                edgecolors='k',
                marker='o',
                alpha=1.0,
                s=120
            )

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(
        r'intermediate diameter [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_ylabel('self-diffusivity [cm$^2$s$^{-1}$]', fontsize=16)
    ax.set_xlim(2, 6)
    ax.set_ylim(1E-18, 1E-3)
    ax.set_yscale("log", nonposy='clip')
    fig.tight_layout()
    fig.savefig(
        "verploegh2015_cf.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def categorical_with_known(molecules, known_df, threshold, output_dir):
    """
    Categorical scatter plot considering experimental results.

    """
    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))
    for name in molecules:
        out_file = (
            f"{output_dir}/{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        lit_d = known_df[
            known_df['molecule'] == name
        ]['diffuse'].iloc[0]
        if lit_d == 't':
            if mid_diam <= threshold:
                C = 'b'
                M = 'o'
                D = 0.25
            else:
                C = 'b'
                M = 'o'
                D = 0.25
        elif lit_d == 'f':
            if mid_diam <= threshold:
                C = 'r'
                M = 'o'
                D = 0.75
            else:
                C = 'r'
                M = 'o'
                D = 0.75
        else:
            continue
        ax.scatter(
            D+(dx*(np.random.random() - 0.5) * 2),
            mid_diam,
            c=C,
            edgecolors='k',
            marker=M,
            alpha=1.0,
            s=80)

    ax.axhspan(ymin=3.2, ymax=threshold, facecolor='k', alpha=0.2)

    pfn.define_diff_categ_plot(
        ax,
        title='',
        xtitle='',
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xlim=(0, 1),
        ylim=(0, 10)
    )
    fig.tight_layout()
    fig.savefig("categorical.pdf", dpi=720, bbox_inches='tight')


def seed_test(seeds):
    """
    Compares the minimum diameter obtained for a set of molecules with
    different random seeds for the ETKDG algorithm.

    """

    molecules = {
        'n-hexane': 'CCCCCC',
        'n-heptane': 'CCCCCCC',
        'n-octane': 'CCCCCCCC',
        'toluene': 'CC1=CC=CC=C1',
        'p-nitrophenol': 'C1=CC(=CC=C1[N+](=O)[O-])O',
        'p-nitrophenyl butyrate': 'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]',
        'butyric acid': 'CCCC(=O)O',
    }
    colours = {
        'n-hexane': 'k',
        'n-heptane': 'r',
        'n-octane': 'b',
        'toluene': 'green',
        'p-nitrophenol': 'purple',
        'p-nitrophenyl butyrate': 'orange',
        'butyric acid': 'darkgray',
    }
    markers = {
        'n-hexane': 'o',
        'n-heptane': 'X',
        'n-octane': 'D',
        'toluene': 'P',
        'p-nitrophenol': '^',
        'p-nitrophenyl butyrate': '>',
        'butyric acid': '<',
    }

    seed_output = "seed_test.pkl"

    if os.path.exists(seed_output):
        # load results
        full_results = pickle.load(open(seed_output, "rb"))
    else:
        full_results = {}
        for t in seeds:
            full_results[t] = {}
            for name in molecules:
                full_results[t][name] = {}

        for name in molecules:
            for t in seeds:
                output_dir = f'seeds_{t}'
                out_file = (
                    f"{output_dir}/"
                    f"{name.replace(' ', '_').replace('/', '__')}"
                    '_diam_result.csv'
                )
                if os.path.exists(out_file) is False:
                    continue
                results = pd.read_csv(out_file)
                if len(results) == 0:
                    continue
                min_diam_avg = np.average(results['diam1'])
                min_diam_std = np.std(results['diam1'])
                mid_diam_avg = np.average(results['diam2'])
                mid_diam_std = np.std(results['diam2'])
                min_mid = min(results['diam2'])
                result = (
                    min_diam_avg,
                    min_diam_std,
                    mid_diam_avg,
                    mid_diam_std,
                    min_mid
                )
                full_results[t][name] = result
        # save file
        pickle.dump(full_results, open("seed_test.pkl", "wb"))

    fig, ax = plt.subplots()
    for name in molecules:
        X = []
        Y = []
        for t in seeds:
            RES = full_results[t][name]
            _, _, _, _, min_mid = RES
            X.append(int(t))
            Y.append(min_mid)
        ax.scatter(X, Y, c=colours[name], marker=markers[name],
                   label=name)
    t_lim = (0, 850000)
    t_name = 'random seed'
    pfn.define_standard_plot(
        ax,
        xtitle=t_name,
        ytitle=r'$d$ [$\mathrm{\AA}$]',
        xlim=t_lim,
        ylim=(4, 8)
    )
    # ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
    # ax.set_xticklabels([str(i) for i in seeds])
    # ax.legend(fontsize=16, ncol=3)
    fig.tight_layout()
    fig.savefig("min_of_mid_seeds.pdf", dpi=720, bbox_inches='tight')


def parameter_tests(molecules):
    """
    Function that scans over multiple parameters for a set of molecules
    and outputs figures.

    """

    markers = {
        'carbon dioxide': 'o',
        'n-butane': 'X',
        'n-hexane': 'D',
        'n-heptane': 'o',
        'n-octane': 'P',
        'ethanol': 'o',
        'para-xylene': '^',
        'meta-xylene': 'o',
        'toluene': '<',
        'napthalene': '>'
    }
    colours = {
        'carbon dioxide': 'k',
        'n-butane': 'r',
        'n-hexane': 'purple',
        'n-heptane': 'orange',
        'n-octane': 'k',
        'ethanol': 'k',
        'para-xylene': 'b',
        'meta-xylene': 'b',
        'toluene': 'green',
        'napthalene': 'darkgray'
    }
    # calculated molecule properties
    # (MW, no heavy atoms, no rotatable bonds)
    properties = {
        'carbon dioxide': (43.98982924, 3, 0),
        'n-butane': (58.078250319999995, 4, 1),
        'n-hexane': (86.109550448, 6, 3),
        'n-heptane': (100.12520051199999, 7, 4),
        'n-octane': (114.14085057599999, 8, 5),
        'ethanol': (46.041864812, 3, 0),
        'para-xylene': (106.07825032, 8, 0),
        'meta-xylene': (106.07825032, 8, 0),
        'toluene': (92.062600256, 7, 0),
        'napthalene': (128.062600256, 10, 0)
    }

    # Parameter tests.
    parameter_sets = {
        'spacing': [0.3, 0.4, 0.5, 0.6],
        'N_conformers': [10, 50, 100, 200, 300, 400, 600, 1000],
        'boxMargin': [4, 6, 8]
    }

    test_mol = [
        'n-butane', 'meta-xylene', 'n-hexane', 'n-heptane',
        'n-octane', 'toluene', 'napthalene'
    ]

    param_test_output = "param_test.pkl"

    if os.path.exists(param_test_output):
        # load results
        full_results = pickle.load(open(param_test_output, "rb"))
    else:
        full_results = {}
        for t in parameter_sets:
            full_results[t] = {}
            for name in molecules:
                if name in test_mol:
                    full_results[t][name] = {}

        for name in molecules:
            if name not in test_mol:
                continue
            for t in parameter_sets:
                for v in parameter_sets[t]:
                    output_dir = f'ptests_{t}_{v}'
                    out_file = (
                        f"{output_dir}/"
                        f"{name.replace(' ', '_').replace('/', '__')}"
                        '_diam_result.csv'
                    )
                    if os.path.exists(out_file) is False:
                        continue
                    results = pd.read_csv(out_file)
                    if len(results) == 0:
                        continue
                    min_diam_avg = np.average(results['diam1'])
                    min_diam_std = np.std(results['diam1'])
                    mid_diam_avg = np.average(results['diam2'])
                    mid_diam_std = np.std(results['diam2'])
                    min_mid = min(results['diam2'])
                    result = (
                        min_diam_avg,
                        min_diam_std,
                        mid_diam_avg,
                        mid_diam_std,
                        min_mid
                    )
                    full_results[t][name][v] = result
        # save file
        pickle.dump(full_results, open(param_test_output, "wb"))

    min_plots(
        parameter_sets=parameter_sets,
        molecules=molecules,
        test_mol=test_mol,
        full_results=full_results,
        colours=colours,
        markers=markers
    )
    mid_plots(
        parameter_sets=parameter_sets,
        molecules=molecules,
        test_mol=test_mol,
        full_results=full_results,
        colours=colours,
        markers=markers
    )
    min_of_mid_plots(
        parameter_sets=parameter_sets,
        molecules=molecules,
        test_mol=test_mol,
        full_results=full_results,
        colours=colours,
        markers=markers
    )
    target_conformer_plot(
        parameter_sets=parameter_sets,
        molecules=molecules,
        test_mol=test_mol,
        full_results=full_results,
        colours=colours,
        markers=markers,
        properties=properties
    )


def min_plots(
    parameter_sets,
    molecules,
    test_mol,
    full_results,
    colours,
    markers
):
    for t in parameter_sets:
        fig, ax = plt.subplots()
        for name in molecules:
            if name not in test_mol:
                continue
            X = []
            Y = []
            Y_err = []
            for i, v in enumerate(parameter_sets[t]):
                RES = full_results[t][name][v]
                min_diam_avg, min_diam_std, _, _, _ = RES
                avg = float(min_diam_avg)
                std = float(min_diam_std)
                # if i == 0:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #             yerr=std, fmt=markers[name], label=name)
                # else:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #             yerr=std, fmt=markers[name])
                X.append(float(v))
                Y.append(avg)
                Y_err.append(std)
            ax.plot(
                X, Y, c=colours[name], marker=markers[name], label=name
            )
        if t == 'N_conformers':
            t_lim = (0, 1100)
            t_name = '$N$'  # 'no. conformers'
        if t == 'spacing':
            t_lim = (0, 1.2)
            t_name = r'grid spacing [$\mathrm{\AA}$]'
        if t == 'vdw':
            t_lim = (0.4, 1.2)
            t_name = r'vdW scale parameter'
        if t == 'boxMargin':
            t_lim = (2, 10)
            t_name = r'box margin [$\mathrm{\AA}$]'
        pfn.define_standard_plot(
            ax,
            xtitle=t_name,
            ytitle=r'avg. minimum diameter [$\mathrm{\AA}$]',
            xlim=t_lim,
            ylim=(0, 10)
        )
        ax.legend(loc=1, fontsize=16)
        fig.tight_layout()
        fig.savefig(f"min_{t}.pdf", dpi=720, bbox_inches='tight')


def mid_plots(
    parameter_sets,
    molecules,
    test_mol,
    full_results,
    colours,
    markers
):
    for t in parameter_sets:
        fig, ax = plt.subplots()
        for name in molecules:
            if name not in test_mol:
                continue
            X = []
            Y = []
            Y_err = []
            for i, v in enumerate(parameter_sets[t]):
                RES = full_results[t][name][v]
                _, _, mid_diam_avg, mid_diam_std, _ = RES
                avg = float(mid_diam_avg)
                std = float(mid_diam_std)
                # if i == 0:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #         yerr=std, fmt=markers[name], label=name)
                # else:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #         yerr=std, fmt=markers[name])
                X.append(float(v))
                Y.append(avg)
                Y_err.append(std)
            X = np.asarray(X)
            Y = np.asarray(Y)
            Y_err = np.asarray(Y_err)
            ax.plot(
                X, Y,
                c=colours[name],
                marker=markers[name],
                label=name
            )
            ax.fill_between(
                X, Y-Y_err, Y+Y_err,
                alpha=0.2,
                facecolor=colours[name]
            )
        if t == 'N_conformers':
            t_lim = (0, 1100)
            t_name = '$N$'  # 'no. conformers'
        if t == 'spacing':
            t_lim = (0.2, 1.1)
            t_name = r'grid spacing [$\mathrm{\AA}$]'
        if t == 'vdw':
            t_lim = (0.4, 1.1)
            t_name = 'vdW scale parameter'
        if t == 'boxMargin':
            t_lim = (3, 9)
            t_name = r'box margin [$\mathrm{\AA}$]'
        pfn.define_standard_plot(
            ax,
            xtitle=t_name,
            ytitle=r'avg. intermediate diameter [$\mathrm{\AA}$]',
            xlim=t_lim,
            ylim=(3.5, 9)
        )
        # ax.legend(fontsize=16, ncol=2)
        fig.tight_layout()
        fig.savefig(f"mid_{t}.pdf", dpi=720, bbox_inches='tight')


def min_of_mid_plots(
    parameter_sets,
    molecules,
    test_mol,
    full_results,
    colours,
    markers
):
    for t in parameter_sets:
        fig, ax = plt.subplots()
        for name in molecules:
            if name not in test_mol:
                continue
            X = []
            Y = []
            Y_err = []
            for i, v in enumerate(parameter_sets[t]):
                _, _, _, _, min_mid = full_results[t][name][v]
                # if i == 0:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #         yerr=std, fmt=markers[name], label=name)
                # else:
                #     ax.errorbar(float(v), avg, c=colours[name],
                #         yerr=std, fmt=markers[name])
                X.append(float(v))
                Y.append(min_mid)
            X = np.asarray(X)
            Y = np.asarray(Y)
            Y_err = np.asarray(Y_err)
            ax.plot(
                X, Y,
                c=colours[name],
                marker=markers[name],
                label=name
            )
        if t == 'N_conformers':
            t_lim = (0, 1100)
            t_name = '$N$'  # 'no. conformers'
        if t == 'spacing':
            t_lim = (0.2, 0.7)
            t_name = r'grid spacing [$\mathrm{\AA}$]'
        if t == 'vdw':
            t_lim = (0.4, 1.1)
            t_name = 'vdW scale parameter'
        if t == 'boxMargin':
            t_lim = (3, 9)
            t_name = r'box margin [$\mathrm{\AA}$]'
        pfn.define_standard_plot(
            ax,
            xtitle=t_name,
            ytitle=r'$d$ [$\mathrm{\AA}$]',
            xlim=t_lim,
            ylim=(3.5, 8)
        )
        # ax.legend(fontsize=16, ncol=3)
        fig.tight_layout()
        fig.savefig(
            f"min_of_mid_{t}.pdf",
            dpi=720,
            bbox_inches='tight'
        )


def target_conformer_plot(
    parameter_sets,
    molecules,
    test_mol,
    full_results,
    colours,
    markers,
    properties
):
    # target no conformers
    targ_confs = [50, 200]
    # set property
    for p, PROP in enumerate(['MW', 'NHA', 'NRB']):
        if p == 0:
            PROP_lab = 'MW [g/mol]'
            p_lim = (0, 120)
        if p == 1:
            PROP_lab = 'no. heavy atoms'
            p_lim = (0, 9)
        if p == 2:
            PROP_lab = 'no. rotatable bonds'
            p_lim = (0, 6)
        for t in parameter_sets:
            if t != 'N_conformers':
                continue
            # fig = plt.figure()  # figsize=(8, 8))
            # ax = fig.add_subplot(111, projection='3d')
            fig, ax = plt.subplots()
            for name in molecules:
                if name not in test_mol:
                    continue
                X = []
                Y = []
                Z = []
                for i, v in enumerate(parameter_sets[t]):
                    RES = full_results[t][name][v]
                    _, _, _, _, min_mid = RES
                    # if i == 0:
                    #     ax.errorbar(
                    #     float(v), avg, c=colours[name],
                    #     yerr=std, fmt=markers[name], label=name)
                    # else:
                    #     ax.errorbar(float(v), avg, c=colours[name],
                    #                 yerr=std, fmt=markers[name])
                    X.append(float(v))
                    Y.append(min_mid)
                    Z.append(properties[name][p])
                X = np.asarray(X)
                Y = np.asarray(Y)
                Z = np.asarray(Z)
                for targ_conf in targ_confs:
                    Y2 = Y-Y[-1]
                    Z2 = Z[X == targ_conf]
                    Y2 = Y2[X == targ_conf]
                    # plot points
                    # ax.scatter(X, Y-Y[-1], Z, s=60,
                    #            c=colours[name], marker=markers[name])
                    ax.scatter(
                        Z2, Y2,
                        c=colours[name],
                        marker=markers[name],
                        label=name,
                        s=80
                    )

            pfn.define_standard_plot(
                ax,
                xtitle=PROP_lab,
                # ytitle=(
                #     '$d_{\mathrm{i, min}}$ - '
                #     '$d_{\mathrm{i, min}}$(1000) '
                #     '[$\mathrm{\AA}$]'
                # ),
                ytitle=r'$d-d$(1000) [$\mathrm{\AA}$]',
                xlim=p_lim,
                ylim=(-0.1, 0.5)
            )
            ax.axhline(y=0, c='k', linestyle='--')
            # ax.set_xlabel(t_name, fontsize=16)
            # ax.set_ylabel(
            #     '$d_{\mathrm{i, min}}-d_{\mathrm{i, min}}$(1000)'
            #     ' [$\mathrm{\AA}$]'
            # ),
            #               fontsize=16)
            # ax.set_zlabel(PROP_lab, fontsize=16)
            # ax.set_xlim(t_lim)
            # ax.set_ylim(-0.1, 0.5)
            # ax.set_zlim(p_lim)
            # ax.set_aspect('equal', 'box')
            # dist = 30
            # angles = 10
            # ax.view_init(dist, angles)
            # ax.legend(fontsize=14, ncol=2)
            fig.tight_layout()
            fig.savefig(
                f"min_of_mid_{t}_v_prop_{PROP}.pdf",
                bbox_inches='tight',
                dpi=720
            )


def shapes_with_known(molecules, known_df, threshold, output_dir):
    """
    Plot molecule shapes considering experimental results.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name in molecules:
        out_file = (
            f"{output_dir}/{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        if len(results) == 0:
            continue
        mid_diam = min(results['diam2'])
        lit_d = known_df[
            known_df['molecule'] == name
        ]['diffuse'].iloc[0]
        if lit_d == 't':
            if mid_diam <= threshold:
                C = 'b'
                M = 'o'
            else:
                C = 'b'
                M = 'X'
        elif lit_d == 'f':
            if mid_diam <= threshold:
                C = 'r'
                M = 'X'
            else:
                C = 'r'
                M = 'o'
        else:
            continue
        ax.scatter(
            np.average(results['ratio_1']),
            np.average(results['ratio_2']),
            c=C,
            edgecolors='k',
            marker=M,
            alpha=1.0,
            s=80
        )

    ax.plot([0, 0.5, 1, 0], [1, 0.5, 1, 1], c='k', lw=2)
    ax.text(0.75, 1.03, 'sphere', fontsize=20)
    ax.text(0.4, 0.45, 'oblate', fontsize=20)
    ax.text(-0.05, 1.03, 'prolate', fontsize=20)

    pfn.define_standard_plot(
        ax,
        xtitle='$I_1$ / $I_3$',
        ytitle='$I_2$ / $I_3$',
        xlim=(-0.1, 1.1),
        ylim=(0.4, 1.1)
    )
    fig.tight_layout()
    fig.savefig("shape.pdf", dpi=720, bbox_inches='tight')


def main():
    if (not len(sys.argv) == 5):
        print("""
    Usage: param_screening.py

        molecule_file :
            psp_source/data/test_molecules.txt

        rerun_diameter_calc :
            t or f

        param_file :
            psp_source/data/param_file.txt

        do_parity :
            t or f to make the figure or not

        """)
        sys.exit()
    else:
        molecule_file = sys.argv[1]
        rerun_diameter_calc = True if sys.argv[2] == 't' else False
        pars = utilities.read_params(sys.argv[3])
        do_parity = True if sys.argv[4] == 't' else False

    start = time.time()

    df, molecules, diameters = rdkf.read_mol_txt_file(
        molecule_file
    )

    # Ignore MW restrictions.
    pars['MW_thresh'] = 2000

    # draw 2D structures
    print('--- draw 2D structures...')
    rdkf.draw_svg_for_all_molecules(molecules)

    parameter_sets = {
        'spacing': [0.3, 0.4, 0.5, 0.6],
        'N_conformers': [10, 50, 100, 200, 300, 400, 600, 1000],
        'boxMargin': [4, 6, 8]
    }
    seeds = [
        1, 1000, 500, 50000, 2123, 345555, 542221, 679293,
        2755, 99982, 825412, 342, 54638, 1982, 77654, 8553, 4
    ]

    # calculate the size of the ellipsoid surroudning all molecules
    # using input pars
    if rerun_diameter_calc:
        print('--- calculating molecular diameters for all tests...')
        rdkf.calc_molecule_diameters(
            molecules,
            pars=pars,
            out_dir='orig_pars',
        )

        # Scale test.
        new_pars = pars.copy()
        new_pars['vdwScale'] = 1.0
        rdkf.calc_molecule_diameters(
            molecules,
            pars=new_pars,
            out_dir='scale1_test',
        )

        new_pars = pars.copy()
        new_pars['vdwScale'] = 0.9
        rdkf.calc_molecule_diameters(
            molecules,
            pars=new_pars,
            out_dir='scale09_test',
        )

        # Scale test.
        new_pars = pars.copy()
        new_pars['vdwScale'] = 0.8
        rdkf.calc_molecule_diameters(
            molecules,
            pars=new_pars,
            out_dir='scale08_test',
        )

        # Seed test.
        print('--------- seed tests! ----------------')
        for seed in seeds:
            new_pars = pars.copy()
            new_pars['seed'] = seed
            print(f'doing seed {seed}')
            print(new_pars)
            new_molecules = {
                'n-hexane': 'CCCCCC',
                'n-heptane': 'CCCCCCC',
                'n-octane': 'CCCCCCCC',
                'toluene': 'CC1=CC=CC=C1',
                'p-nitrophenol': 'C1=CC(=CC=C1[N+](=O)[O-])O',
                'p-nitrophenyl butyrate': (
                    'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                'butyric acid': 'CCCC(=O)O',
            }
            print(new_molecules)
            rdkf.calc_molecule_diameters(
                new_molecules,
                pars=new_pars,
                out_dir=f'seeds_{seed}',
            )

        # Parameter tests.
        print('--------- param tests! ----------------')
        test_mol = [
            'n-butane', 'meta-xylene', 'n-hexane', 'n-heptane',
            'n-octane', 'toluene', 'napthalene'
        ]

        new_molecules = {
            i: molecules[i] for i in molecules if i in test_mol
        }
        print(new_molecules)

        for t in parameter_sets:
            for v in parameter_sets[t]:
                print(f'doing test {t} with value {v}')
                new_pars = pars.copy()
                new_pars[t] = v

                rdkf.calc_molecule_diameters(
                    new_molecules,
                    pars=new_pars,
                    out_dir=f'ptests_{t}_{v}',
                )

    cf_verploegh2015(molecules, output_dir='orig_pars')

    print_results_cf_known(
        molecules,
        known_df=df,
        threshold=pars['size_thresh'],
        output_dir='orig_pars'
    )

    parity_with_known(molecules, diameters, output_dir='orig_pars')

    categorical_with_known(
        molecules,
        known_df=df,
        threshold=pars['size_thresh'],
        output_dir='orig_pars'
    )

    shapes_with_known(
        molecules,
        known_df=df,
        threshold=pars['size_thresh'],
        output_dir='orig_pars'
    )

    if do_parity:
        scale_info = {
            # DIR: (scale, C, M, alpha, edgecolor)
            'scale1_test': (1.0, 'k', 'o', 0.5, 'k'),
            'scale09_test': (0.9, 'b', 'D', 0.5, 'k'),
            'scale08_test': (0.8, 'r', 'P', 0.5, 'k')
        }
        parity_cf_scale_with_known(
            molecules,
            diameters,
            known_df=df,
            pars=pars,
            scale_info=scale_info
        )

    seed_test(seeds=seeds)

    parameter_tests(
        molecules
    )

    end = time.time()
    print(f'---- total time taken = {round(end-start, 2)} s')


if __name__ == "__main__":
    main()
