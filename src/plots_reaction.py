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
import matplotlib.colors as colors
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

    ax2 = ax.twinx()
    ax2.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        alpha=0.9, width=width,
        color='#2C3E50',
        edgecolor='k'
    )

    # cumulative plot
    cumul = np.cumsum(hist)
    ax.plot(
        bin_edges[:-1],
        cumul,
        alpha=1.0,
        label='max component < threshold',
        color='r',
        marker='o'
    )

    # ax.axvspan(xmin=4.0, xmax=6.6, facecolor='k', alpha=0.2,
    #    hatch="/")
    ax.axvspan(xmin=4.0, xmax=6.6, facecolor='k', alpha=0.2)
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    # ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    pfn.define_standard_plot(
        ax,
        xtitle=r'$d$ of largest component [$\mathrm{\AA}$]',
        ytitle='cumulative # reactions',
        xlim=(0, 17),
        ylim=(0, int(max(cumul)+max(cumul)*0.1))
    )
    ax2.set_ylim(0, int(max(hist)+max(hist)*0.2))
    ax2.set_ylabel('# reactions', fontsize=16)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Change left y axis colours.
    ax.spines['left'].set_color('red')
    ax2.spines['left'].set_color('red')

    ax2.tick_params(axis='both', which='major', labelsize=16)
    fig.tight_layout()
    fig.savefig(
        f"{plot_suffix}/size_threshold_{plot_suffix}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def rxn_space(data, filename):
    """
    Plot number of possible reactions as a function of size threshold.

    """

    plot_prop = {
        1: {
            'c': '#FA7268',
            'e': 'none',
            'a': 0.5,
            'm': 'o',
            's': 50,
            'label': 'class I'
        },
        2: {
            'c': '#DAF7A6',
            'e': 'none',
            'a': 0.5,
            'm': 'x',
            's': 50,
            'label': 'class II'
        },
        3: {
            'c': '#900C3F',
            'e': 'none',
            'a': 1.0,
            'm': 'x',
            's': 50,
            'label': 'class III'
        },
        4: {
            'c': '#F6D973',
            'e': 'none',
            'a': 0.5,
            'm': 'x',
            's': 50,
            'label': 'class IV'
        }
    }

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(0, 20.5, width)
    fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    for p in plot_prop:
        if p != 3:
            continue
        pp = plot_prop[p]
        sub_data = data[data['PC_class'] == p]

        hist, bin_edges = np.histogram(
            a=sub_data['max_mid_diam'],
            bins=X_bins
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            align='edge',
            alpha=pp['a'],
            width=width,
            color=pp['c'],
            edgecolor='k',
            label=pp['label']
        )

    ax.legend(fontsize=16)
    ax.axvspan(xmin=4.0, xmax=6.6, facecolor='k', alpha=0.2, hatch="/")
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    # HOF.
    ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    pfn.define_standard_plot(
        ax,
        xtitle=r'$d$ of largest component [$\mathrm{\AA}$]',
        ytitle='# reactions',
        xlim=(0, 17),
        ylim=None
    )

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )


def rxn_value(data, filename):
    """
    Plot the value of all reactions as violin plot.

    """

    plot_prop = {
        1: {
            'c': '#900C3F',
            'e': 'none',
            'a': 0.5,
            'm': 'o',
            's': 50,
            'label': 'class I'
        },
        2: {
            'c': '#FA7268',
            'e': 'none',
            'a': 0.5,
            'm': 'x',
            's': 50,
            'label': 'class II'
        },
        3: {
            'c': '#F6D973',
            'e': 'none',
            'a': 1.0,
            'm': 'x',
            's': 50,
            'label': 'class III'
        },
        4: {
            'c': '#DAF7A6',
            'e': 'none',
            'a': 0.5,
            'm': 'x',
            's': 50,
            'label': 'class IV'
        }
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    for p in plot_prop:
        pp = plot_prop[p]
        sub_data = data[data['PC_class'] == p]

        values = sub_data['max_mid_diam']
        number = int(p)
        parts = ax.violinplot(
            values,
            [number],
            showmeans=False,
            showmedians=False,
            showextrema=False
        )
        for pc in parts['bodies']:
            pc.set_facecolor(pp['c'])
            pc.set_edgecolor('black')
            pc.set_alpha(1.0)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('purchasability class', fontsize=16)
    ax.set_ylabel(
        r'$d$ of largest component [$\mathrm{\AA}$]',
        fontsize=16
    )
    ax.set_xlim(0.5, 4.5)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=720,
        bbox_inches='tight'
    )


def rxn_complexity(data, filename):
    """
    Plot the measures of complexity of each reaction.

    """

    fig, ax = plt.subplots(figsize=(8, 5))
    ylim = (-1000, 1000)
    xlim = (-10, 10)
    # CS = [(1.0, 1.0, 1.0), (44/255, 62/255, 80/255)]
    # cm = colors.LinearSegmentedColormap.from_list('test', CS, N=10)
    # fig, ax, hist = pfn.twoD_histogram(
    #     X_data=data['deltasa'],
    #     Y_data=data['deltabct'],
    #     xlim=xlim,
    #     ylim=ylim,
    #     cmap=cm,
    #     fig=fig,
    #     ax=ax
    # )
    # cbar = fig.colorbar(hist[3], ax=ax)
    # cbar.ax.set_ylabel('count', fontsize=16)
    # cbar.ax.tick_params(labelsize=16)
    ax.scatter(
        data['deltasa'],
        data['deltabct'],
        c='#CCD1D1',
        edgecolors='none',
        marker='o',
        alpha=1.0,
        s=40,
        label='full dataset'
    )
    small_data = data[data['max_mid_diam'] < 6.6]
    ax.scatter(
        small_data['deltasa'],
        small_data['deltabct'],
        c='#2C3E50',
        edgecolors='none',
        marker='o',
        alpha=1.0,
        s=40,
        label='viable reactions'
    )

    pfn.define_standard_plot(
        ax,
        # xtitle='number of heavy atoms',
        ylim=ylim,
        xlim=xlim,
        ytitle=r'$\Delta$ BertzCT',
        xtitle=r'$\Delta$ SAscore',
    )

    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        filename,
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

    if xtitle == 'purchasability class':
        align = 'center'
    else:
        align = 'edge'

    ax.bar(
        bin_edges[:-1],
        hist,
        align=align,
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


def pie(X, xtitle, xlim, width):
    """
    Plot pie chart of categorical data.

    """
    if xtitle == 'purchasability class':
        labels = ['class I', 'class II', 'class III', 'class IV']
        colours = ['#D2B1D1', '#3498DB', '#C0392B', '#CCD1D1']
        sizes = [
            len([i for i in X if i == 1]),
            len([i for i in X if i == 2]),
            len([i for i in X if i == 3]),
            len([i for i in X if i == 4])
        ]
    else:
        raise ValueError('this type of plot is not defined.')

    # explode = (0.0, 0.0)

    fig, ax = plt.subplots(figsize=(5, 5))
    wedges, _, _ = ax.pie(
        sizes,
        colors=colours,
        # explode=explode,
        labels=labels,
        autopct='%1.1f%%',
        # shadow=True,
        startangle=90,
        textprops={'fontsize': 16}
    )

    for w in wedges:
        w.set_linewidth(1.5)
        w.set_edgecolor('k')
    # Equal aspect ratio ensures that pie is drawn as a circle.
    ax.axis('equal')

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
    if col == 'max_mid_diam':
        ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2)
    ax.set_ylim(ylim)

    return fig, ax
