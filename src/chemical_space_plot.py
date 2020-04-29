#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot the chemical space of molecules.

Author: Andrew Tarzia

Date Created: 16 Feb 2020

"""

from os.path import exists
import sys
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import json
import numpy as np
import pandas as pd

from reaction import KEGG_IDs_to_ignore
import plotting_fn as pfn


def cs_purch(purch, not_purch):
    fig, ax = plt.subplots(figsize=(8, 5))
    plot_prop = {
        't': {
            'c': '#FA7268',
            'e': 'none',
            'a': 0.5,
            'm': 'o',
            's': 50,
            'label': 'purchasable'
        },
        'f': {
            'c': '#DAF7A6',
            'e': 'none',
            'a': 0.5,
            'm': 'x',
            's': 50,
            'label': 'not purchasable'
        }
    }

    # bin each of the sets of data based on X value
    for p in plot_prop:
        pp = plot_prop[p]
        if p == 't':
            data = purch
        else:
            data = not_purch
        width = 0.5
        X_bins = np.arange(0, 15.5, width)
        hist, bin_edges = np.histogram(
            a=data,
            bins=X_bins,
            density=True
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            align='edge',
            alpha=0.8,
            width=width,
            color=pp['c'],
            edgecolor='k',
            label=pp['label'],
        )

    # for X, Y, Z in zip(Xs, Ys, Zs):
    #     if Z:
    #         pp = plot_prop['t']
    #     else:
    #         pp = plot_prop['f']
    #
    #     ax.scatter(
    #         X,
    #         Y,
    #         c=pp['c'],
    #         edgecolors=pp['e'],
    #         marker=pp['m'],
    #         alpha=pp['a'],
    #         s=pp['s']
    #     )

    # Vertical lines for different materials.
    ax.axvspan(xmin=4.0, xmax=6.6, facecolor='k', alpha=0.2, hatch="/")
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    # HOF size limit:
    ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    # # Legend.
    # for p in plot_prop:
    #     pp = plot_prop[p]
    #     ax.scatter(
    #         X,
    #         Y,
    #         c=pp['c'],
    #         edgecolors=pp['e'],
    #         marker=pp['m'],
    #         alpha=pp['a'],
    #         s=pp['s'],
    #         label=pp['label']
    #     )
    ax.legend(fontsize=16)

    pfn.define_standard_plot(
        ax,
        # xtitle='number of heavy atoms',
        xtitle=r'intermediate diameter [$\mathrm{\AA}$]',
        ytitle='frequency',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_purch.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def cs_purchCT(purch, not_purch):
    fig, ax = plt.subplots(figsize=(8, 5))
    plot_prop = {
        't': {
            'c': '#FA7268',
            'e': 'none',
            'a': 0.5,
            'm': 'o',
            's': 50,
            'label': 'purchasable'
        },
        'f': {
            'c': '#DAF7A6',
            'e': 'none',
            'a': 0.5,
            'm': 'x',
            's': 50,
            'label': 'not purchasable'
        }
    }

    # bin each of the sets of data based on X value
    for p in plot_prop:
        pp = plot_prop[p]
        if p == 't':
            data = purch
        else:
            data = not_purch
        width = 50
        X_bins = np.arange(0, 2000, width)
        hist, bin_edges = np.histogram(
            a=data,
            bins=X_bins,
            density=False
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            align='edge',
            alpha=0.8,
            width=width,
            color=pp['c'],
            edgecolor='k',
            label=pp['label'],
        )

    ax.legend(fontsize=16)

    pfn.define_standard_plot(
        ax,
        # xtitle='number of heavy atoms',
        xtitle=r'BertzCT',
        ytitle='frequency',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_purchCT.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def cs_MW(Xs, Ys):
    fig, ax = plt.subplots(figsize=(8, 5))
    ylim = (0, 17)
    xlim = (0, 550)
    CS = [(1.0, 1.0, 1.0), (44/255, 62/255, 80/255)]
    cm = colors.LinearSegmentedColormap.from_list('test', CS, N=10)
    fig, ax, hist = twoD_histogram(
        X_data=Xs,
        Y_data=Ys,
        xlim=xlim,
        ylim=ylim,
        cmap=cm,
        fig=fig,
        ax=ax
    )
    cbar = fig.colorbar(hist[3], ax=ax)
    cbar.ax.set_ylabel('count', fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    # ax.scatter(
    #     Xs,
    #     Ys,
    #     c='#FF7900',
    #     edgecolors='k',
    #     marker='o',
    #     alpha=1.0,
    #     s=40
    # )

    # Horizontal lines for different materials.
    ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2)
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    # HOF size limit:
    # ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    pfn.define_standard_plot(
        ax,
        # xtitle='number of heavy atoms',
        ylim=ylim,
        xlim=xlim,
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xtitle=r'MW [g/mol]',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_MW.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def cs_NHA(Xs, Ys):
    fig, ax = plt.subplots(figsize=(8, 5))
    ylim = (0, 17)
    xlim = (0, 40)
    CS = [(1.0, 1.0, 1.0), (44/255, 62/255, 80/255)]
    cm = colors.LinearSegmentedColormap.from_list('test', CS, N=10)
    fig, ax, hist = twoD_histogram(
        X_data=Xs,
        Y_data=Ys,
        xlim=xlim,
        ylim=ylim,
        cmap=cm,
        fig=fig,
        ax=ax
    )
    cbar = fig.colorbar(hist[3], ax=ax)
    cbar.ax.set_ylabel('count', fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    #
    # ax.scatter(
    #     Xs,
    #     Ys,
    #     c='#FF7900',
    #     edgecolors='k',
    #     marker='o',
    #     alpha=1.0,
    #     s=120
    # )

    # Horizontal lines for different materials.
    ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2)
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    # HOF size limit:
    # ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    pfn.define_standard_plot(
        ax,
        # xtitle='number of heavy atoms',
        ylim=ylim,
        xlim=xlim,
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xtitle=r'no. heavy atoms',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_NHA.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def cs_NRB(Xs, Ys):
    fig, ax = plt.subplots(figsize=(8, 5))
    ylim = (0, 17)
    xlim = (-1, 30)
    CS = [(1.0, 1.0, 1.0), (44/255, 62/255, 80/255)]
    cm = colors.LinearSegmentedColormap.from_list('test', CS, N=10)
    fig, ax, hist = twoD_histogram(
        X_data=Xs,
        Y_data=Ys,
        xlim=xlim,
        ylim=ylim,
        cmap=cm,
        fig=fig,
        ax=ax
    )
    cbar = fig.colorbar(hist[3], ax=ax)
    cbar.ax.set_ylabel('count', fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    # ax.scatter(
    #     Xs,
    #     Ys,
    #     c='#FF7900',
    #     edgecolors='k',
    #     marker='o',
    #     alpha=1.0,
    #     s=120
    # )

    # Horizontal lines for different materials.
    ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2)
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    # HOF size limit:
    # ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    pfn.define_standard_plot(
        ax,
        ylim=ylim,
        xlim=xlim,
        # xtitle='number of heavy atoms',
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xtitle=r'no. rotatable bonds',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_NRB.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def cs_sol(logPs, logSs, HlogPs, HlogSs):
    fig, ax = plt.subplots(figsize=(8, 5))
    ylim = (-13, 4)
    xlim = (-9, 14)
    CS = [(1.0, 1.0, 1.0), (44/255, 62/255, 80/255)]
    cm = colors.LinearSegmentedColormap.from_list('test', CS, N=10)
    fig, ax, hist = twoD_histogram(
        X_data=logPs,
        Y_data=logSs,
        xlim=xlim,
        ylim=ylim,
        cmap=cm,
        fig=fig,
        ax=ax
    )
    cbar = fig.colorbar(hist[3], ax=ax)
    cbar.ax.set_ylabel('count', fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    ax.scatter(
        HlogPs,
        HlogSs,
        c='#E74C3C',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=80
    )

    pfn.define_standard_plot(
        ax,
        ylim=ylim,
        xlim=xlim,
        # xtitle='number of heavy atoms',
        xtitle=r'logP',
        ytitle=r'logS$_{\mathrm{w}}$',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_sol.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def cs_logPvsNHA(logPs, Xs, HlogPs, HXs):
    fig, ax = plt.subplots(figsize=(8, 5))
    xlim = (0, 40)
    ylim = (-9, 14)
    CS = [(1.0, 1.0, 1.0), (44/255, 62/255, 80/255)]
    cm = colors.LinearSegmentedColormap.from_list('test', CS, N=10)
    fig, ax, hist = twoD_histogram(
        X_data=Xs,
        Y_data=logPs,
        xlim=xlim,
        ylim=ylim,
        cmap=cm,
        fig=fig,
        ax=ax
    )
    cbar = fig.colorbar(hist[3], ax=ax)
    cbar.ax.set_ylabel('count', fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    ax.scatter(
        HXs,
        HlogPs,
        c='#E74C3C',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=80
    )

    pfn.define_standard_plot(
        ax,
        ylim=ylim,
        xlim=xlim,
        # xtitle='number of heavy atoms',
        ytitle=r'logP',
        xtitle=r'number of heavy atoms',
    )
    fig.tight_layout()
    fig.savefig(
        f'chemical_space_logPNHA.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def chemical_space_plot():
    """
    Output chemical space plot.

    """

    molecule_list = glob.glob('*_unopt.mol')

    print(f'{len(molecule_list)} molecules in DB.')

    KEGG_IDs_to_highlight = [
        'C01387', 'C00756', 'C00123', 'C00183', 'C00041',
        'C00079', 'C00407', 'C00078', 'C00073', 'C00082'
    ]

    Xs = []
    Ys = []
    MWs = []
    NRBs = []
    Zs = []
    purch = []
    not_purch = []
    purch_CT = []
    not_purch_CT = []
    logPs = []
    logSs = []
    HlogPs = []
    HXs = []
    HlogSs = []
    COUNTER = 0
    for mol in sorted(molecule_list):

        name = mol.replace('_unopt.mol', '')
        etkdg_fail = name+'_unopt.ETKDGFAILED'
        diam_file = name+'_size.csv'
        toobig_file = name+'_size.TOOBIG'
        prop_file = name+'_prop.json'

        if exists(toobig_file):
            continue
        if exists(etkdg_fail):
            continue
        if name in KEGG_IDs_to_ignore():
            continue

        if exists(prop_file) and exists(diam_file):
            # Get molecular properties from 2D structure.
            with open(prop_file, 'r') as f:
                prop_dict = json.load(f)
            # Get size and update output lists.
            results = pd.read_csv(diam_file)
            min_mid_diam = min(results['diam2'])
            Xs.append(prop_dict['NHA'])
            MWs.append(prop_dict['MW'])
            NRBs.append(prop_dict['NRB'])
            Zs.append(prop_dict['purchasability'])
            Ys.append(min_mid_diam)
            logPs.append(prop_dict['logP'])
            logSs.append(prop_dict['logS'])
            if name in KEGG_IDs_to_highlight:
                HlogPs.append(prop_dict['logP'])
                HlogSs.append(prop_dict['logS'])
                HXs.append(prop_dict['NHA'])
            if prop_dict['purchasability']:
                purch.append(min_mid_diam)
                purch_CT.append(prop_dict['bertzCT'])
            else:
                not_purch.append(min_mid_diam)
                not_purch_CT.append(prop_dict['bertzCT'])
            if prop_dict['NRB'] > 10 or min_mid_diam > 10:
                print(name, prop_dict['NRB'], min_mid_diam)
            if min_mid_diam < 6.6:
                COUNTER += 1

    print(f'{COUNTER} with intermediate diameter < 6.6 A')
    # rn [
    # '#FA7268', '#F8A72A', '#DAF7A6', '#900C3F', '#6BADB0',
    # '#DB869D', '#F6D973', 'mediumvioletred'
    cs_purch(purch, not_purch)
    cs_purchCT(purch_CT, not_purch_CT)
    cs_MW(Xs=MWs, Ys=Ys)
    cs_NHA(Xs=Xs, Ys=Ys)
    cs_NRB(Xs=NRBs, Ys=Ys)
    cs_sol(logPs, logSs, HlogPs, HlogSs)
    cs_logPvsNHA(logPs, Xs, HlogPs, HXs)


def twoD_histogram(
    X_data,
    Y_data,
    fig,
    ax,
    cmap,
    xlim,
    ylim
):
    """
    2D-histogram of pair data.

    """

    hist = ax.hist2d(
        X_data,
        Y_data,
        bins=[20, 20],
        range=[xlim, ylim],
        # no normalization
        density=False,
        # log colour map.
        norm=colors.LogNorm(),
        cmap=cmap
    )

    return fig, ax, hist


def main():
    if (not len(sys.argv) == 1):
        print("""
Usage: chemical_space_plot.py

""")
        sys.exit()
    else:
        pass

    chemical_space_plot()


if __name__ == "__main__":
    main()
