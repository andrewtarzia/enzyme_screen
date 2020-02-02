#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting reaction-based functions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import plotting_fn as pfn


def no_rxns_vs_size(data, params, plot_suffix):
    """
    Plot number of possible reactions as a function of size threshold.

    """

    fig, ax = plt.subplots(figsize=(8, 5))

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(0, 20.5, width)
    hist, bin_edges = np.histogram(a=data['max_mid_diam'], bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=0.9, width=width,
        color='#2980B9',
        edgecolor='k'
    )

    # cumulative plot
    cumul = np.cumsum(hist)
    ax.plot(
        bin_edges[:-1],
        cumul,
        alpha=1.0,
        label='max component < threshold',
        color='k',
        marker='o'
    )

    ax.axvspan(xmin=4.0, xmax=6.6, facecolor='k', alpha=0.2, hatch="/")
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    pfn.define_standard_plot(
        ax,
        xtitle=r'$d$ of largest component [$\mathrm{\AA}$]',
        ytitle='# reactions',
        xlim=(0, 17),
        ylim=(0, int(max(cumul)+max(cumul)*0.1))
    )
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.tight_layout()
    fig.savefig(
        f"{plot_suffix}/size_threshold_{plot_suffix}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def save_candidates(data, params, filename):
    """
    Save candidates to file.

    """

    all_fit = data.sort_values(by='max_mid_diam')
    all_fit.to_csv(filename, index=False)

    print(f'There are {len(all_fit)} candidate reactions!')
    print('---------------------------------------------------')


def stacked_dist(data, col, xtitle, xlim, width):
    """
    Plot histograms of data stacked by top level EC no.

    """

    delta_data = {'total': []}

    for i, row in data.iterrows():
        EC = row['ec']
        top_EC = EC.split('.')[0]
        if top_EC not in list(delta_data.keys()):
            delta_data[top_EC] = []
        delta_data[top_EC].append(row[col])
        delta_data['total'].append(row[col])

    fig, ax = plt.subplots(figsize=(8, 5))

    if xlim is None:
        xlim = (
            min([min(delta_data[i]) for i in delta_data])-2*width,
            max([max(delta_data[i]) for i in delta_data])+2*width
        )
    X_bins = np.arange(xlim[0], xlim[1], width)

    for keys in delta_data:
        values = delta_data[keys]
        hist, bin_edges = np.histogram(
            a=values,
            bins=X_bins,
            density=True
        )
        ax.plot(
            X_bins[:-1]+width/2,
            hist,
            c=pfn.EC_descriptions()[keys][1],
            lw='1.5',
            marker='o',
            alpha=1.0,
            label=pfn.EC_descriptions()[keys][0]
        )
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel('frequency', fontsize=16)
    ax.set_xlim(xlim)
    # legend
    ax.legend(fontsize=16)

    return fig, ax


def dist(X, xtitle, xlim, width):
    """
    Plot histograms of data.

    """
    fig, ax = plt.subplots(figsize=(8, 5))

    if xlim is None:
        xlim = (min(X)-2*width, max(X)+2*width)
    X_bins = np.arange(xlim[0], xlim[1], width)
    hist, bin_edges = np.histogram(a=X, bins=X_bins)

    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=1.0,
        width=width,
        color='#2980B9',
        edgecolor='k'
    )
    pfn.define_standard_plot(
        ax,
        xtitle=xtitle,
        ytitle='count',
        xlim=xlim,
        ylim=None
    )

    return fig, ax


def violinplot(data, col, ytitle, ylim):
    """
    Plot violin plots of data separated by top level EC no.

    """

    delta_data = {'total': []}

    for i, row in data.iterrows():
        EC = row['ec']
        top_EC = EC.split('.')[0]
        if top_EC not in list(delta_data.keys()):
            delta_data[top_EC] = []
        delta_data[top_EC].append(row[col])
        delta_data['total'].append(row[col])

    fig, ax = plt.subplots(figsize=(8, 5))

    for keys in delta_data:
        values = delta_data[keys]
        if keys == '-':
            number = 0
        elif keys == 'total':
            number = -1
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
            pc.set_facecolor(pfn.EC_descriptions()[keys][1])
            pc.set_edgecolor('black')
            pc.set_alpha(1.0)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('EC number', fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_xlim(-2, 8)
    ax.set_xticks([-1, 0, 1, 2, 3, 4, 5, 6, 7])
    ax.set_xticklabels(
        ['all', 'unknown', '1', '2', '3', '4', '5', '6', '7']
    )
    ax.set_ylim(ylim)

    return fig, ax
