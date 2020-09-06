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
import glob
import matplotlib.pyplot as plt
import os
from os.path import exists
import json

import IO
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

    pfn.define_diff_categ_plot(
        ax,
        title='',
        xtitle='',
        ytitle=r'intermediate diameter [$\mathrm{\AA}$]',
        xlim=(0, 1),
        ylim=(0, 10),
        xticks=[0.25, 0.75],
        xlabels=['diffuses', 'does not diffuse']
    )
    fig.tight_layout()
    fig.savefig(
        f"categorical_{plot_suffix}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def shapes(molecules, threshold, output_dir, plot_suffix):
    """
    Plot molecule shapes of all molecules in dictionary.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name in molecules:
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

    pfn.define_standard_plot(
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


def mol_parity(propx, propy, file, xtitle, ytitle, mol_file=None):
    """
    Plot a parity of two molecular properties.

    """

    if mol_file is None:
        molecule_list = glob.glob('*_unopt.mol')
    else:
        molecule_list = IO.read_molecule_list(mol_file)

    # iterate over molecules
    Xs = []
    Ys = []
    for mol in molecule_list:
        name = mol.replace('_unopt.mol', '')
        prop_file = name+'_prop.json'

        if not exists(prop_file):
            continue

        with open(prop_file, 'r') as f:
            prop_dict = json.load(f)

        Xs.append(prop_dict[propx])
        Ys.append(prop_dict[propy])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(
        Xs,
        Ys,
        c='#FA7268',
        edgecolors='k',
        marker='o',
        alpha=1.0,
        s=80
    )
    xlim = None
    ylim = None
    if propx == 'Synth_score':
        xlim = (0, 10)
    elif propy == 'Synth_score':
        ylim = (0, 10)

    pfn.define_standard_plot(
        ax,
        xtitle=xtitle,
        ytitle=ytitle,
        xlim=xlim,
        ylim=ylim
    )
    fig.tight_layout()
    fig.savefig(
        f'parity_{file}.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def mol_categ(propx, propy, file, xtitle, ytitle, mol_file=None):
    """
    Plot a categorization of two molecular properties.

    """

    if mol_file is None:
        molecule_list = glob.glob('*_unopt.mol')
    else:
        molecule_list = IO.read_molecule_list(mol_file)

    # iterate over molecules
    Ys = {'true': [], 'false': []}
    for mol in molecule_list:
        name = mol.replace('_unopt.mol', '')
        prop_file = name+'_prop.json'
        size_file = name+'_size.csv'

        if not exists(prop_file):
            continue
        if not exists(size_file):
            continue
        if exists(size_file.replace('.csv', '.TOOBIG')):
            continue
        if exists(size_file.replace(
            'size.csv',
            'unopt.ETKDGFAILED'
        )):
            continue

        with open(prop_file, 'r') as f:
            prop_dict = json.load(f)

        if propy == 'size':
            results = pd.read_csv(size_file)
            size = min(results['diam2'])
            Y = size
        else:
            Y = prop_dict[propy]

        if propx == 'size':
            results = pd.read_csv(size_file)
            size = min(results['diam2'])
            if size > 6.6:
                Ys['true'].append(Y)
            else:
                Ys['false'].append(Y)
        elif prop_dict[propx]:
            Ys['true'].append(Y)
        else:
            Ys['false'].append(Y)

    fig, ax = plt.subplots(figsize=(4, 5))

    for keys in Ys:
        values = Ys[keys]
        if keys == 'true':
            number = 1
            C = '#2C3E50'
        elif keys == 'false':
            number = 0
            C = '#CCD1D1'
        else:
            number = int(keys)
        parts = ax.violinplot(
            values,
            [number],
            showmeans=False,
            showmedians=False,
            showextrema=False
        )
        for pc in parts['bodies']:
            pc.set_facecolor(C)
            pc.set_edgecolor('black')
            pc.set_alpha(1.0)

    xlim = (-1, 2)
    ylim = None
    if propy == 'Synth_score':
        ylim = (1, 10)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([0, 1])

    if propx == 'purchasability':
        ax.set_xticklabels(
            ['\nnot purchasable', 'purchasable']
        )
    elif propx == 'size':
        ax.set_xticklabels(
            ['can diffuse', '\ncannot diffuse']
        )
    ax.set_ylim(ylim)
    fig.tight_layout()
    fig.savefig(
        f'categ_{file}.pdf',
        dpi=720,
        bbox_inches='tight'
    )


def mol_pie(data_dict):
    """
    Plot distribution of a molecular property.

    """

    labels = ['purchasable', 'not purchasable']
    colours = [data_dict['c'], 'gray']
    sizes = [
        len([i for i in data_dict['d'] if i is True]),
        len([i for i in data_dict['d'] if i is False])
    ]
    # explode = (0.0, 0.0)

    fig, ax = plt.subplots(figsize=(5, 5))
    wedges, _, _ = ax.pie(
        sizes,
        colors=colours,
        # explode=explode,
        labels=labels,
        autopct='%1.1f%%',
        # shadow=True,
        startangle=90
    )

    for w in wedges:
        w.set_linewidth(1.5)
        w.set_edgecolor('k')
    # Equal aspect ratio ensures that pie is drawn as a circle.
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(
        f"cat_{data_dict['file']}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def mol_dist(data_dict):
    """
    Plot distribution of a molecular property.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = data_dict['width']
    X_bins = np.arange(
        data_dict['xlim'][0],
        data_dict['xlim'][1],
        width
    )
    hist, bin_edges = np.histogram(a=data_dict['d'], bins=X_bins)

    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color=data_dict['c'],
        edgecolor='k'
    )
    pfn.define_standard_plot(
        ax,
        xtitle=data_dict['xtitle'],
        ytitle='count',
        xlim=data_dict['xlim'],
        ylim=None
    )
    fig.tight_layout()
    fig.savefig(
        f"hist_{data_dict['file']}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def mol_all_dist(plot_suffix, mol_file=None):
    """
    Plot distributions of molecule attributes.

    """

    if mol_file is None:
        molecule_list = glob.glob('*_unopt.mol')
    else:
        molecule_list = IO.read_molecule_list(mol_file)

    prop_to_plot = {
        'logP': {
            'd': [],
            'width': 0.5,
            'xlim': (-20, 20),
            'xtitle': 'logP',
            'c': '#FA7268',
            'file': f'logP_{plot_suffix}'
        },
        'logS': {
            'd': [],
            'width': 0.5,
            'xlim': (-20, 20),
            'xtitle': 'logS',
            'c': '#DAF7A6',
            'file': f'logS_{plot_suffix}'
        },
        'Synth_score': {
            'd': [],
            'width': 0.25,
            'xlim': (0, 10),
            'xtitle': 'SAScore',
            'c': '#6BADB0',
            'file': f'SA_{plot_suffix}'
        },
        'purchasability': {
            'd': [],
            'width': 0.25,
            'xlim': (0, 10),
            'xtitle': 'pur',
            'c': '#5499C7',
            'file': f'purch_{plot_suffix}'
        },
        'MW': {
            'd': [],
            'width': 100,
            'xlim': (0, 10000),
            'xtitle': 'molecular weight [g/mol]',
            'c': '#5499C7',
            'file': f'MW_{plot_suffix}'
        }
    }
    for mol in molecule_list:
        name = mol.replace('_unopt.mol', '')
        prop_file = name+'_prop.json'

        if not exists(prop_file):
            continue

        with open(prop_file, 'r') as f:
            prop_dict = json.load(f)

        print('>>>>>', mol)
        prop_to_plot['logP']['d'].append(prop_dict['logP'])
        prop_to_plot['MW']['d'].append(prop_dict['MW'])
        prop_to_plot['logS']['d'].append(prop_dict['logS'])
        prop_to_plot['Synth_score']['d'].append(
            prop_dict['Synth_score']
        )
        prop_to_plot['purchasability']['d'].append(
            prop_dict['purchasability']
        )

    # do plots
    for prop in prop_to_plot:
        d = prop_to_plot[prop]
        if prop == 'purchasability':
            mol_pie(d)
        else:
            mol_dist(d)
