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
import json
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

import IO
import rdkit_functions as rdkf
import plots_molecular as pm
import utilities
import chemcost_IO
from reaction import KEGG_IDs_to_ignore
import plotting_fn as pfn


def chemical_space_plot():
    """
    Output chemical space plot.

    """

    molecule_list = glob.glob('*_unopt.mol')

    print(f'{len(molecule_list)} molecules in DB.')

    Xs = []
    Ys = []
    Zs = []
    purch = []
    not_purch = []
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
            Zs.append(prop_dict['purchasability'])
            Ys.append(min_mid_diam)
            if prop_dict['purchasability']:
                purch.append(min_mid_diam)
            else:
                not_purch.append(min_mid_diam)

        # rn [
        # '#FA7268', '#F8A72A', '#DAF7A6', '#900C3F', '#6BADB0',
        # '#DB869D', '#F6D973', 'mediumvioletred'

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
        f'chemical_space.pdf',
        dpi=720,
        bbox_inches='tight'
    )


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
