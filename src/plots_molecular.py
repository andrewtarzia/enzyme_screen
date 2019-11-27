#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting molecule-based functions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os

import plotting_fn as pfn


def print_results(molecules, threshold, output_dir):
    """
    Print calculated ability to diffuse for all molecules in dict.

    """
    diffuse = {}
    no_diffuse = {}
    for name in molecules:
        smile = molecules[name]
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        print(out_file)
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        min_diam = min(results['diam1'])
        mid_diam = min(results['diam2'])
        if mid_diam <= threshold:
            print(name+':')
            print('can diffuse')
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


def categorical(molecules, threshold, output_dir, plot_suffix):
    """
    Categorical scatter plot of all molecules in dictionary.

    """
    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))
    for name in molecules:
        smile = molecules[name]
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        if mid_diam <= threshold:
            C = 'b'
            M = 'o'
            E = 'k'
            D = 0.25
        else:
            C = 'r'
            M = 'o'
            E = 'k'
            D = 0.75

        ax.scatter(
            D+(dx*(np.random.random() - 0.5) * 2),
            mid_diam,
            c=C,
            edgecolors=E,
            marker=M,
            alpha=1.0,
            s=80
        )

    ax.axhline(y=threshold, c='k')

    define_diff_categ_plot(
        ax,
        title='',
        xtitle='',
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xlim=(0, 1),
        ylim=(0, 10)
    )
    fig.tight_layout()
    fig.savefig(
        f"categorical_{plot_suffix}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def categorical_moloutput(
    mol_output_file, threshold, output_dir, plot_suffix
):
    """
    Categorical scatter plot of all molecules molecule_output file.

    """
    import DB_functions
    molecule_output = DB_functions.initialize_mol_output_DF(
        mol_output_file,
        overwrite=False
    )
    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))
    for idx, row in molecule_output.iterrows():
        mid_diam = row['mid_diam']
        if mid_diam == 0:
            continue
        if mid_diam <= threshold:
            M = 'o'
            E = 'k'
            D = 0.25
            print('unique molecule that fits:', row['name'],
                  '- DB:', row['DB'], '- ID:', row['DB_ID'])
        else:
            M = 'o'
            E = 'k'
            D = 0.75

        # set colour based on role
        if row['role'] == 'reactant':
            C = 'b'
        elif row['role'] == 'product':
            C = 'r'
        elif row['role'] == 'both':
            C = 'purple'

        ax.scatter(D+(dx*(np.random.random() - 0.5) * 2),
                   mid_diam, c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=80)

    # decoy legend
    ax.scatter(-100, 100,
               c='b',
               edgecolors=E, marker='o', alpha=1.0,
               s=100,
               label='reactant')
    ax.scatter(-100, 100,
               c='r',
               edgecolors=E, marker='o', alpha=1.0,
               s=100,
               label='product')
    ax.scatter(-100, 100,
               c='purple',
               edgecolors=E, marker='o', alpha=1.0,
               s=100,
               label='either')
    ax.legend(loc=2, fontsize=12)
    ax.axhline(y=threshold, c='k')
    define_diff_categ_plot(
        ax,
        title='',
        xtitle='',
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xlim=(0, 1),
        ylim=(0, 10)
    )
    fig.tight_layout()
    fig.savefig(output_dir+"categorical_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def shapes(molecules, threshold, output_dir, plot_suffix):
    """
    Plot molecule shapes of all molecules in dictionary.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name in molecules:
        smile = molecules[name]
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        if mid_diam <= threshold:
            C = 'b'
            M = 'o'
            E = 'k'
        else:
            C = 'r'
            M = 'o'
            E = 'k'

        ax.scatter(np.average(results['ratio_1']),
                   np.average(results['ratio_2']),
                   c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=80)

    ax.plot([0, 0.5, 1, 0], [1, 0.5, 1, 1], c='k', lw=2)
    ax.text(0.75, 1.03, 'sphere', fontsize=20)
    ax.text(0.4, 0.45, 'oblate', fontsize=20)
    ax.text(-0.05, 1.03, 'prolate', fontsize=20)

    define_standard_plot(
        ax,
        title='',
        xtitle='$I_1$ / $I_3$',
        ytitle='$I_2$ / $I_3$',
        xlim=(-0.1, 1.1),
        ylim=(0.4, 1.1)
    )
    fig.tight_layout()
    fig.savefig(
        f"shape_{plot_suffix}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def mol_SA_vs_compl(output_dir, plot_suffix):
    """Plot the synthetic accessibility of a molecules VS its complexity.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    # iterate over molecules
    for m in yield_molecules(directory=output_dir):
        K_count = 0
        for R in m.rs_pkls:
            if 'KEGG' in R:
                K_count += 1
        if K_count == 0:
            continue
        if m.Synth_score is None:  # or mol.Synth_score == 0:
            continue
        if m.complexity is None:  # or mol.complexity == 0:
            continue
        M = 'o'
        E = 'k'
        ax.scatter(m.complexity,
                   m.Synth_score,
                   c='orange',
                   edgecolors=E,
                   marker=M,
                   alpha=0.8,
                   s=40)
    define_standard_plot(ax,
                         title='',
                         xtitle='complexity',
                         ytitle='SAscore',
                         xlim=(0, 5000.1),
                         ylim=(0, 10.1))
    fig.tight_layout()
    fig.savefig(output_dir+"SA_VS_compl_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def mol_SA_vs_NHA(output_dir, plot_suffix):
    """
    Plot the synthetic accessibility of a molecules VS its no. heavy
    atoms.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    # iterate over molecules
    for m in yield_molecules(directory=output_dir):
        K_count = 0
        for R in m.rs_pkls:
            if 'KEGG' in R:
                K_count += 1
        if K_count == 0:
            continue
        if m.Synth_score is None:  # or mol.Synth_score == 0:
            continue
        if m.mol is None:
            continue
        M = 'o'
        E = 'k'
        # no. heavy atoms

        ax.scatter(NHA,
                   m.Synth_score,
                   c='orange',
                   edgecolors=E,
                   marker=M,
                   alpha=0.8,
                   s=40)
    define_standard_plot(ax,
                         title='',
                         xtitle='no. heavy atoms',
                         ytitle='SAscore',
                         xlim=(0, 175.1),
                         # xlim=(0, 7.1),
                         ylim=(0, 10.1))
    fig.tight_layout()
    fig.savefig(output_dir+"SA_VS_NHA_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def mol_SA_vs_NRB(output_dir, plot_suffix):
    """
    Plot the synthetic accessibility of a molecules VS its no.
    rotatable bonds.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    # iterate over molecules
    for m in yield_molecules(directory=output_dir):
        K_count = 0
        for R in m.rs_pkls:
            if 'KEGG' in R:
                K_count += 1
        if K_count == 0:
            continue
        if m.Synth_score is None:  # or mol.Synth_score == 0:
            continue
        if m.mol is None:
            continue
        M = 'o'
        E = 'k'
        # no. rotatable bonds

        ax.scatter(NRB,
                   m.Synth_score,
                   c='purple',
                   edgecolors=E,
                   marker=M,
                   alpha=1.0,
                   s=60)
    define_standard_plot(ax,
                         title='',
                         xtitle='no. rotatable bonds',
                         ytitle='SAscore',
                         xlim=(0, 100.1),
                         ylim=(0, 10.1))
    fig.tight_layout()
    fig.savefig(output_dir+"SA_VS_NRB_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def mol_logP_vs_logS(output_dir, plot_suffix):
    """
    Plot the logP VS logS of all molecules.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    # iterate over molecules
    for m in yield_molecules(directory=output_dir):
        K_count = 0
        for R in m.rs_pkls:
            if 'KEGG' in R:
                K_count += 1
        if K_count == 0:
            continue
        if m.logS is None:  # or mol.Synth_score == 0:
            continue
        if m.logP is None:
            continue
        M = 'o'
        E = 'k'
        ax.scatter(m.logP,
                   m.logS,
                   c='orange',
                   edgecolors=E,
                   marker=M,
                   alpha=0.8,
                   s=40)
    define_standard_plot(ax,
                         title='',
                         xtitle='logP',
                         ytitle='logS',
                         xlim=(-20, 30),
                         ylim=(-30, 10))
    fig.tight_layout()
    fig.savefig(output_dir+"logS_VS_logP_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def mol_logP_vs_XlogP(output_dir, plot_suffix):
    """
    Plot the logP VS XlogP (PubChem) of all molecules.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    # iterate over molecules
    for m in yield_molecules(directory=output_dir):
        K_count = 0
        for R in m.rs_pkls:
            if 'KEGG' in R:
                K_count += 1
        if K_count == 0:
            continue
        if m.XlogP is None:
            continue
        if m.logP is None:
            continue
        M = 'o'
        E = 'k'
        ax.scatter(m.logP,
                   m.XlogP,
                   c='orange',
                   edgecolors=E,
                   marker=M,
                   alpha=0.8,
                   s=40)
    ax.plot(np.linspace(-40, 40, 2), np.linspace(-40, 40, 2), c='k', alpha=0.4)
    define_standard_plot(ax,
                         title='',
                         xtitle='logP',
                         ytitle='XlogP3-AA',
                         xlim=(-20, 30),
                         ylim=(-20, 30))
    fig.tight_layout()
    fig.savefig(output_dir+"logP_VS_XlogP_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def mol_all_dist(output_dir, plot_suffix):
    """
    Plot distributions of molecule attributes.

    """
    prop_to_plot = {'logP': [], 'logS': [], 'SAscore': []}
    for m in yield_molecules(directory=output_dir):
        if m.Synth_score == 0 or m.Synth_score is None:
            continue
        if m.logS is None:
            continue
        if m.logP is None:
            continue
        prop_to_plot['logP'].append(m.logP)
        prop_to_plot['logS'].append(m.logS)
        prop_to_plot['SAscore'].append(m.Synth_score)

    # do plots
    mol_all_logP(output_dir, prop_to_plot['logP'], plot_suffix)
    mol_all_logS(output_dir, prop_to_plot['logS'], plot_suffix)
    mol_all_SA(output_dir, prop_to_plot['SAscore'], plot_suffix)


def mol_all_logP(output_dir, data, plot_suffix):
    """
    Plot distribution of logP.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.5
    X_bins = np.arange(-20, 20, width)
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='purple',
           edgecolor='k')
    dist_plot(fig, ax, name='all_logP', xlim=(-20, 20),
              xtitle='logP', plot_suffix=plot_suffix)


def mol_all_logS(output_dir, data, plot_suffix):
    """
    Plot distribution of logS.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.5
    X_bins = np.arange(-20, 20, width)
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='purple',
           edgecolor='k')
    dist_plot(fig, ax, name='all_logS', xlim=(-20, 20),
              xtitle='logS', plot_suffix=plot_suffix)


def mol_all_SA(output_dir, data, plot_suffix):
    """
    Plot distribution of SAscore.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.1
    X_bins = np.arange(0, 10, width)
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='purple',
           edgecolor='k')
    dist_plot(fig, ax, name='all_SA', xlim=(0, 10),
              xtitle='SAscore', plot_suffix=plot_suffix)
