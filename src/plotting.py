#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting functions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors
from ercollect.molecule import yield_molecules, molecule
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator


def EC_descriptions():
    """Dictionary of EC descriptions + colours.

    """
    top_tier = {'-': ('unknown', '#1469b5'),
                '1': ('oxidoreductases', '#FF7900'),
                '2': ('transferases', '#00B036'),
                '3': ('hydrolases', '#EB0000'),
                '4': ('lyases', '#A440BC'),
                '5': ('isomerases', '#945348'),
                '6': ('ligases', '#FA4BBE')}

    return top_tier


def define_diff_categ_plot(ax, title, ytitle, xtitle, xlim, ylim):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_ylabel(ytitle, fontsize=16)
    # ax.legend(
    #     [y, n], ['aligned', 'not aligned'], loc=4, fancybox=True
    # )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticklabels(['diffuses', 'does not diffuse'])
    ax.set_xticks([0.25, 0.75])


def dist_plot(fig, ax, name, xlim, xtitle, plot_suffix):
    """Standard plot properties for distributions.

    """
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(xlim)
    # legend
    # ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig("dist_"+name+"_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def define_standard_plot(ax, title, ytitle, xtitle, xlim, ylim):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    # ax.legend(
    #     [y, n], ['aligned', 'not aligned'], loc=4, fancybox=True
    # )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def define_3d_plot(
    ax, title, xtitle, ytitle, ztitle, xlim, ylim, zlim
):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_ylabel(ztitle, fontsize=16)
    # ax.legend(
    #     [y, n], ['aligned', 'not aligned'], loc=4, fancybox=True
    # )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_ylim(zlim)


def print_results(molecules, threshold, output_dir):
    """Print calculated ability to diffuse for all molecules in dictionary.

    """
    diffuse = {}
    no_diffuse = {}
    for name, smile in molecules.items():
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        print(out_file)
        if os.path.isfile(out_file) is False:
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
    """Categorical scatter plot of all molecules in dictionary.

    """
    dx = 0.15
    fig, ax = plt.subplots(figsize=(5, 5))
    for name, smile in molecules.items():
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        if os.path.isfile(out_file) is False:
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

        ax.scatter(D+(dx*(np.random.random() - 0.5) * 2),
                   mid_diam, c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=80)

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


def categorical_moloutput(
    mol_output_file, threshold, output_dir, plot_suffix
):
    """Categorical scatter plot of all molecules molecule_output file.

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
                        ytitle='intermediate diameter [$\mathrm{\AA}$]',
                        xlim=(0, 1),
                        ylim=(0, 10))
    fig.tight_layout()
    fig.savefig(output_dir+"categorical_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def shapes(molecules, threshold, output_dir, plot_suffix):
    """Plot molecule shapes of all molecules in dictionary.

    """
    fig, ax = plt.subplots(figsize=(5, 5))
    for name, smile in molecules.items():
        out_file = output_dir+name.replace(' ', '_')+'_diam_result.csv'
        if os.path.isfile(out_file) is False:
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
                        ylim=(0.4, 1.1))
    fig.tight_layout()
    fig.savefig(output_dir+"shape_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')


def check_rxn_unique(reaction_reported, rs):
    """Check (using KEGG identifiers (as molecule.name)) if a rxn is
    unique.

    """
    r_n = []
    for r in rs.components:
        r_n.append(r.name)
    r_n = sorted(r_n)
    if r_n in reaction_reported:
        unique = False
    else:
        reaction_reported.append(r_n)
        unique = True
    return unique, reaction_reported


def rs_number_rxns_vs_size(output_dir, size_thresh, generator, plot_suffix):
    """Plot number of possible and unique reactions as a function of size
    threshold.

    pI count commented out.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    max_sizes = []
    reaction_reported = []
    num_duplicate = 0
    data = []
    # also plot the number of new reactions with pI < thresh
    # max_sizes_pI = []
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if rs.max_comp_size < 4.5:
            print(rs.pkl)
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            num_duplicate += 1
            continue
        try:
            if rs.max_comp_size > 0:
                max_sizes.append(rs.max_comp_size)
                data.append(rs.max_comp_size)
                # if rs.seed_MOF is True:
                #     max_sizes_pI.append(rs.max_comp_size)
        except AttributeError:
            pass

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(0, 20.5, width)
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='#1469b5',
           edgecolor='k')

    max_sizes = np.asarray(max_sizes)
    # max_sizes_pI = np.asarray(max_sizes_pI)
    # counts = []
    # counts_pI = []
    # for thr in X_bins:
    #     count_above = len(max_sizes[max_sizes < thr])
    #     counts.append(count_above)
    #     count_above_pI = len(max_sizes_pI[max_sizes_pI < thr])
    #     counts_pI.append(count_above_pI)
    # cumulative plot
    cumul = np.cumsum(hist)
    ax.plot(bin_edges[:-1], cumul, alpha=1.0,
            label='max component < threshold',
            color='k', marker='o')
    # ax.bar(threshs, counts, align='center', alpha=0.5, width=0.2,
    #        label='max component < threshold',
    #        color='b', edgecolor='k')
    # ax.bar(threshs, counts_pI, align='center', alpha=0.5, width=0.2,
    #        label='+ pI < '+str(pI_thresh),
    #        color='r', edgecolor='k')

    # ax.legend(loc=2, fontsize=12)

    # ax.axvline(x=3.4, c='k')
    ax.axvspan(xmin=4.0, xmax=6.6, facecolor='k', alpha=0.2, hatch="/")
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    # plot possible region of ZIF pore limiting diameters from
    # Banerjee 2008 - 10.1126/science.1152516
    # ax.axvspan(0.0, 13, facecolor='#2ca02c', alpha=0.2)
    ax.axvline(x=13.1, c='k', lw=2, linestyle='--')

    define_standard_plot(ax,
                         title='',
                         xtitle='$d$ of largest component [$\mathrm{\AA}$]',
                         ytitle='# reactions',
                         xlim=(0, 17),
                         ylim=(0, int(max(cumul)+max(cumul)*0.1)))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.tight_layout()
    fig.savefig(output_dir+"size_threshold_"+plot_suffix+".pdf", dpi=720,
                bbox_inches='tight')
    print('number duplicates:', num_duplicate)


def print_new_rxns(output_dir, generator):
    """Print all new possible and unique reactions that fit.

    """
    reaction_reported = []
    count = 0
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            continue
        if rs.all_fit is True:
            count += 1
            print("------ New Reaction:")
            rs.print_rxn_system()

    print("There are", count, "new reactions!")


def rs_dist_delta_size(output_dir, generator, plot_suffix):
    """Plot change in maximum size of reactants to products.

    """
    delta_data = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        max_react_size = 0
        max_prod_size = 0
        for m in rs.components:
            if m.mid_diam is None:
                continue
            if m.role == 'reactant':
                max_react_size = max([max_react_size, m.mid_diam])
            elif m.role == 'product':
                max_prod_size = max([max_prod_size, m.mid_diam])
        if max_prod_size > 0 and max_react_size > 0:
            delta_size = max_prod_size - max_react_size
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta_data.keys()):
                delta_data[top_EC] = []
            delta_data[top_EC].append(delta_size)
    fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    width = 0.3
    X_bins = np.arange(-7, 7.2, width)
    for keys, values in delta_data.items():
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        # ax.bar(bin_edges[:-1],
        #        hist,
        #        align='edge',
        #        alpha=0.4, width=0.5,
        #        color=EC_descriptions()[keys][1],
        #        edgecolor='k',
        #        label=EC_descriptions()[keys][0])
        ax.plot(X_bins[:-1]+width/2, hist, c=EC_descriptions()[keys][1],
                lw='1.5', alpha=1.0)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('maximum product size $-$ maximum reactant size [$\mathrm{\AA}$]',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(-7, 7)
    # legend
    # ax.legend(fontsize=16)
    fig.tight_layout()
    filename = output_dir+"dist_delta_size_"
    # filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
    filename += plot_suffix+".pdf"
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')


def rs_dist_logP(output_dir, generator, plot_suffix, extreme):
    """Plot distribution of min/max logP of all reactions.

    """
    if extreme != 'min' and extreme != 'max':
        import sys
        sys.exit('requires extreme == max or min')

    data = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if extreme == 'min':
            Y = rs.min_logP
        else:
            Y = rs.max_logP
        if Y is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(data.keys()):
                data[top_EC] = []
            data[top_EC].append(Y)
    # fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    width = 0.25
    X_bins = np.arange(-40, 40.2, width)
    for keys, values in data.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=width,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])
        # ax.plot(X_bins[:-1]+width/2, hist, c=EC_descriptions()[keys][1],
        #         lw='1.5', alpha=1.0)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(extreme+'. logP of all components',
                      fontsize=16)
        ax.set_ylabel('count', fontsize=16)
        ax.set_xlim(-5, 15)
        # legend
        ax.legend(fontsize=16)
        fig.tight_layout()
        filename = output_dir+"dist_"+extreme+"_logP_"
        filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
        # filename += plot_suffix+".pdf"
        fig.savefig(filename,
                    dpi=720, bbox_inches='tight')


def rs_dist_logS(output_dir, generator, plot_suffix, extreme):
    """Plot distribution of min/max logS of all reactions.

    """
    if extreme != 'min' and extreme != 'max':
        import sys
        sys.exit('requires extreme == max or min')

    data = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if extreme == 'min':
            Y = rs.min_logS
        else:
            Y = rs.max_logS
        if Y is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(data.keys()):
                data[top_EC] = []
            data[top_EC].append(Y)
    # fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    width = 0.25
    X_bins = np.arange(-40, 40.2, width)
    for keys, values in data.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=width,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])
        # ax.plot(X_bins[:-1]+width/2, hist, c=EC_descriptions()[keys][1],
        #         lw='1.5', alpha=1.0)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel(extreme+'. logS of all components',
                      fontsize=16)
        ax.set_ylabel('count', fontsize=16)
        ax.set_xlim(-15, 5)
        # legend
        ax.legend(fontsize=16)
        fig.tight_layout()
        filename = output_dir+"dist_"+extreme+"_logS_"
        filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
        # filename += plot_suffix+".pdf"
        fig.savefig(filename,
                    dpi=720, bbox_inches='tight')


def rs_dist_no_reactants(output_dir, generator, plot_suffix):
    """Plot distribution of the number of reactants in all reactions.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    no_reacts = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        nr = 0
        for m in rs.components:
            if m.role == 'reactant':
                nr += 1
        if nr > 0:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(no_reacts.keys()):
                no_reacts[top_EC] = []
            no_reacts[top_EC].append(nr)

    for keys, values in no_reacts.items():
        ax.hist(values,
                facecolor=EC_descriptions()[keys][1],
                alpha=0.4,
                histtype='stepfilled',
                bins=np.arange(-10, 10.2, 0.5),
                label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('no. reactants', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 5)
    # legend
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(output_dir+"dist_no_reacts_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_no_products(output_dir, generator, plot_suffix):
    """Plot distribution of the number of products in all reactions.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    no_reacts = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        nr = 0
        for m in rs.components:
            if m.role == 'product':
                nr += 1
        if nr > 0:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(no_reacts.keys()):
                no_reacts[top_EC] = []
            no_reacts[top_EC].append(nr)

    for keys, values in no_reacts.items():
        ax.hist(values,
                facecolor=EC_descriptions()[keys][1],
                alpha=0.4,
                histtype='stepfilled',
                bins=np.arange(-10, 10.2, 0.5),
                label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('no. products', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 5)
    # legend
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(output_dir+"dist_no_prods_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_delta_SA(output_dir, generator, plot_suffix):
    """Plot distribution of the change in synthetic accesibility from react to
    products.

    """
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if rs.delta_sa is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta.keys()):
                delta[top_EC] = []
            delta[top_EC].append(rs.delta_sa)

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(-10, 10.2, width)
    for keys, values in delta.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=width,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('$\Delta$ SAscore', fontsize=16)
        ax.set_ylabel('count', fontsize=16)
        ax.set_xlim(-10, 10)
        # legend
        ax.legend(fontsize=16)
        fig.tight_layout()
        filename = output_dir+"dist_delta_SA_"
        filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
        fig.savefig(filename,
                    dpi=720, bbox_inches='tight')


def rs_dist_delta_complexity(output_dir, generator, plot_suffix):
    """Plot distribution of the change in complexity from react to
    products.

    """
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if rs.delta_comp is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta.keys()):
                delta[top_EC] = []
            delta[top_EC].append(rs.delta_comp)
    # bin each of the sets of data based on X value
    X_bins = np.arange(-500, 500, 25)
    for keys, values in delta.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=25,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('$\Delta$ complexity', fontsize=16)
        ax.set_ylabel('count', fontsize=16)
        ax.set_xlim(-500, 500)
        # legend
        ax.legend(fontsize=16)
        fig.tight_layout()
        filename = output_dir+"dist_delta_compl_"
        filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
        fig.savefig(filename,
                    dpi=720, bbox_inches='tight')


def rs_dist_max_size(output_dir, generator, plot_suffix):
    """Plot distribution of max component size.

    """
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if rs.max_comp_size is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta.keys()):
                delta[top_EC] = []
            delta[top_EC].append(rs.max_comp_size)

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(0, 20.5, width)
    for keys, values in delta.items():
        fig, ax = plt.subplots(figsize=(8, 5))
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=width,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('$d$ of largest component [$\mathrm{\AA}$]', fontsize=16)
        ax.set_ylabel('count', fontsize=16)
        ax.set_xlim(0, 20)
        # legend
        ax.legend(fontsize=16)
        fig.tight_layout()
        filename = output_dir+"dist_max_size_"
        filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
        fig.savefig(filename,
                    dpi=720, bbox_inches='tight')


def rs_dist_max_size_noEC(output_dir, generator, plot_suffix):
    """Plot distribution of max component size.

    """
    data = []
    num_duplicate = 0
    reaction_reported = []
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            num_duplicate += 1
            continue
        if rs.max_comp_size is not None:
            data.append(rs.max_comp_size)

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(0, 20.5, width)
    fig, ax = plt.subplots(figsize=(8, 5))
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='#1469b5',
           edgecolor='k')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('$d$ of largest component [$\mathrm{\AA}$]', fontsize=16)
    # ax.set_ylabel('count', fontsize=16)
    ax.set_ylabel('# reactions', fontsize=16)
    ax.set_xlim(0, 15)
    # legend
    # ax.legend(fontsize=16)
    fig.tight_layout()
    filename = output_dir+"dist_max_size_noEC.pdf"
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')
    print('number duplicates:', num_duplicate)


def rs_dist_logP_noEC(output_dir, generator, plot_suffix, extreme):
    """Plot distribution of min/max logP of all reactions.

    """
    if extreme != 'min' and extreme != 'max':
        import sys
        sys.exit('requires extreme == max or min')

    data = []
    num_duplicate = 0
    reaction_reported = []
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            num_duplicate += 1
            continue
        if extreme == 'min':
            Y = rs.min_logP
        else:
            Y = rs.max_logP
        if Y is not None:
            data.append(Y)
    # bin each of the sets of data based on X value
    width = 0.25
    X_bins = np.arange(-40, 40.2, width)
    fig, ax = plt.subplots(figsize=(8, 5))
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='#1469b5',
           edgecolor='k')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(extreme+'. logP of all components',
                  fontsize=16)
    # ax.set_ylabel('count', fontsize=16)
    ax.set_ylabel('# reactions', fontsize=16)
    ax.set_xlim(-5, 15)
    fig.tight_layout()
    filename = output_dir+"dist_"+extreme+"_logP_noEC.pdf"
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')
    print('number duplicates:', num_duplicate)


def rs_dist_delta_SA_noEC(output_dir, generator, plot_suffix):
    """Plot distribution of the change in synthetic accesibility from react to
    products.

    """
    data = []
    num_duplicate = 0
    reaction_reported = []
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            num_duplicate += 1
            continue
        if rs.delta_sa is not None:
            data.append(rs.delta_sa)

    # bin each of the sets of data based on X value
    width = 0.5
    X_bins = np.arange(-10, 10.5, width)
    fig, ax = plt.subplots(figsize=(8, 5))
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='#1469b5',
           edgecolor='k')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('$\Delta$ SAscore', fontsize=16)
    # ax.set_ylabel('count', fontsize=16)
    ax.set_ylabel('# reactions', fontsize=16)
    ax.set_xlim(-10, 10)
    fig.tight_layout()
    filename = output_dir+"dist_delta_SA_noEC.pdf"
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')
    print('number duplicates:', num_duplicate)


def rs_dist_logS_noEC(output_dir, generator, plot_suffix, extreme):
    """Plot distribution of min/max logS of all reactions.

    """
    if extreme != 'min' and extreme != 'max':
        import sys
        sys.exit('requires extreme == max or min')

    data = []
    num_duplicate = 0
    reaction_reported = []
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            num_duplicate += 1
            continue
        if extreme == 'min':
            Y = rs.min_logS
        else:
            Y = rs.max_logS
        if Y is not None:
            data.append(Y)
    # fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    width = 0.25
    X_bins = np.arange(-40, 40.2, width)
    fig, ax = plt.subplots(figsize=(8, 5))
    hist, bin_edges = np.histogram(a=data, bins=X_bins)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color='#1469b5',
           edgecolor='k')
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(extreme+'. logS of all components',
                  fontsize=16)
    # ax.set_ylabel('count', fontsize=16)
    ax.set_ylabel('# reactions', fontsize=16)
    ax.set_xlim(-15, 5)
    fig.tight_layout()
    filename = output_dir+"dist_"+extreme+"_logS_noEC.pdf"
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')
    print('number duplicates:', num_duplicate)


def rs_dist_GRAVY(output_dir, generator, plot_suffix):
    """Plot distribution of protein GRAVY for reactions with known sequences.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        try:
            if rs.GRAVY is not None:
                top_EC = rs.EC.split('.')[0]
                if top_EC not in list(delta.keys()):
                    delta[top_EC] = []
                delta[top_EC].append(rs.GRAVY)
        except AttributeError:
            pass

    # bin each of the sets of data based on X value
    X_bins = np.arange(-2, 2.2, 0.1)
    for keys, values in delta.items():
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=0.1,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])

    # GRAVY specific visuals
    # ax.text(-1.45, 40, 'hydrophilic', fontsize=16)
    # get ylim
    ylim = ax.get_ylim()
    ax.text(0.55, max(ylim)/2 + 0.05*(max(ylim)/2), 'hydrophobic', fontsize=16)
    ax.arrow(0.5, max(ylim)/2, 0.7, 0,
             head_width=0.05*(max(ylim)/2), head_length=0.1, fc='k', ec='k')
    avg_GRAVY = -0.4
    ax.axvline(x=avg_GRAVY, c='grey', alpha=1.0, linestyle='--')
    catalase_GRAVY = -0.605
    ax.axvline(x=catalase_GRAVY, c='r', alpha=1.0)
    urease_GRAVY = -0.1524
    ax.axvline(x=urease_GRAVY, c='b', alpha=1.0)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('GRAVY', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(-1.5, 1.5)
    # legend
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(output_dir+"dist_GRAVY_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_A_index(output_dir, generator, plot_suffix):
    """Plot distribution of protein aliphatic indec for reactions with known
    sequences.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        try:
            if rs.A_index is not None:
                top_EC = rs.EC.split('.')[0]
                if top_EC not in list(delta.keys()):
                    delta[top_EC] = []
                delta[top_EC].append(rs.A_index)
        except AttributeError:
            pass

    # bin each of the sets of data based on X value
    X_bins = np.arange(0, 150, 5)
    for keys, values in delta.items():
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=5,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])

    # AI specific visuals
    ylim = ax.get_ylim()
    ax.text(10, max(ylim)/2 + 0.05*(max(ylim)/2), 'more stable', fontsize=16)
    ax.arrow(10, max(ylim)/2, 40, 0,
             head_width=0.05*(max(ylim)/2), head_length=5, fc='k', ec='k')
    catalase_AI = 68
    ax.axvline(x=catalase_AI, c='r', alpha=1.0)
    urease_AI = 90.476
    ax.axvline(x=urease_AI, c='b', alpha=1.0)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('aliphatic index', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 150)
    # legend
    ax.legend(fontsize=16, loc=2)
    fig.tight_layout()
    fig.savefig(output_dir+"dist_A_index_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_pI(output_dir, generator, plot_suffix):
    """Plot distribution of protein pI (no modifications) for reactions with
    known sequences.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        try:
            if rs.pI is not None:
                top_EC = rs.EC.split('.')[0]
                if top_EC not in list(delta.keys()):
                    delta[top_EC] = []
                delta[top_EC].append(rs.pI)
        except AttributeError:
            pass

    # bin each of the sets of data based on X value
    X_bins = np.arange(0, 14.1, 0.5)
    for keys, values in delta.items():
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=0.5,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('pI', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 14)
    # legend
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(output_dir+"dist_pI_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_I_index(output_dir, generator, plot_suffix):
    """Plot distribution of protein aliphatic indec for reactions with known
    sequences.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        try:
            if rs.I_index is not None:
                top_EC = rs.EC.split('.')[0]
                if top_EC not in list(delta.keys()):
                    delta[top_EC] = []
                delta[top_EC].append(rs.I_index)
        except AttributeError:
            pass

    # bin each of the sets of data based on X value
    X_bins = np.arange(0, 150, 5)
    for keys, values in delta.items():
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=5,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])

    # instability specific visuals
    # get ylim
    ylim = ax.get_ylim()
    ax.text(41, max(ylim)/2 + 0.05*(max(ylim)/2), 'unstable', fontsize=16)
    ax.arrow(40, max(ylim)/2, 30, 0,
             head_width=0.05*(max(ylim)/2), head_length=5, fc='k', ec='k')
    II_cutoff = 40
    ax.axvline(x=II_cutoff, c='grey', alpha=1.0, linestyle='--')
    catalase_II = 27.010
    ax.axvline(x=catalase_II, c='r', alpha=1.0)
    urease_II = 31.75
    ax.axvline(x=urease_II, c='b', alpha=1.0)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('instability index', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 100)
    # legend
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(output_dir+"dist_I_index_"+plot_suffix+".pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_TM_index(output_dir, generator, plot_suffix):
    """Plot distribution of protein aliphatic indec for reactions with known
    sequences.

    """
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        try:
            if rs.TM_index is not None:
                top_EC = rs.EC.split('.')[0]
                if top_EC not in list(delta.keys()):
                    delta[top_EC] = []
                delta[top_EC].append(rs.TM_index)
        except AttributeError:
            pass
    fig, ax = plt.subplots(figsize=(8, 5))
    # bin each of the sets of data based on X value
    X_bins = np.arange(-5, 5, 0.5)
    for keys, values in delta.items():
        hist, bin_edges = np.histogram(a=values, bins=X_bins)
        ax.bar(bin_edges[:-1],
               hist,
               align='edge',
               alpha=0.4, width=0.5,
               color=EC_descriptions()[keys][1],
               edgecolor='k',
               label=EC_descriptions()[keys][0])
    # melting temperature index specific visuals
    TM_cutoff = (0, 1)
    ax.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
               alpha=0.2)
    catalase_TMI = 1.22
    ax.axvline(x=catalase_TMI, c='r', alpha=1.0)
    urease_TMI = 0.62
    ax.axvline(x=urease_TMI, c='b', alpha=1.0)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('thermostability index', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    # ax.set_xlim(-5, 5)
    # legend
    ax.legend(fontsize=16)
    fig.tight_layout()
    filename = output_dir+"dist_TM_index_"
    # filename += EC_descriptions()[keys][0]+"_"+plot_suffix+".pdf"
    filename += plot_suffix+".pdf"
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')


# def rs_dist_TM_index_1fig(output_dir, generator, plot_suffix):
#     """Plot distribution of protein aliphatic indec for reactions with known
#     sequences.
#
#     """
#     fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(8, 10))
#     # Remove horizontal space between axes
#     fig.subplots_adjust(hspace=0)
#     delta = {}
#     # iterate over reaction system files
#     for rs in generator:
#         if rs.skip_rxn is True:
#             continue
#         try:
#             if rs.TM_index is not None:
#                 top_EC = rs.EC.split('.')[0]
#                 if top_EC not in list(delta.keys()):
#                     delta[top_EC] = []
#                 delta[top_EC].append(rs.TM_index)
#         except AttributeError:
#             pass
#
#     # bin each of the sets of data based on X value
#     X_bins = np.arange(-5, 5, 0.5)
#     max3 = 0
#     for keys, values in delta.items():
#         if keys != list(delta.keys())[0]:
#             continue
#         hist, bin_edges = np.histogram(a=values, bins=X_bins)
#         max3 = max([max3, max(hist)])
#         ax3.bar(bin_edges[:-1],
#                 hist,
#                 align='edge',
#                 alpha=0.4, width=0.5,
#                 color=EC_descriptions()[keys][1],
#                 edgecolor='k',
#                 label=EC_descriptions()[keys][0])
#     max2 = 0
#     for keys, values in delta.items():
#         if keys != list(delta.keys())[1]:
#             continue
#         hist, bin_edges = np.histogram(a=values, bins=X_bins)
#         max2 = max([max2, max(hist)])
#         ax2.bar(bin_edges[:-1],
#                 hist,
#                 align='edge',
#                 alpha=0.4, width=0.5,
#                 color=EC_descriptions()[keys][1],
#                 edgecolor='k',
#                 label=EC_descriptions()[keys][0])
#     max1 = 0
#     for keys, values in delta.items():
#         if keys != list(delta.keys())[2]:
#             continue
#         hist, bin_edges = np.histogram(a=values, bins=X_bins)
#         max1 = max([max1, max(hist)])
#         ax1.bar(bin_edges[:-1],
#                 hist,
#                 align='edge',
#                 alpha=0.4, width=0.5,
#                 color=EC_descriptions()[keys][1],
#                 edgecolor='k',
#                 label=EC_descriptions()[keys][0])
#
#     ax1.tick_params(axis='y', which='major', labelsize=16)
#     ax2.tick_params(axis='y', which='major', labelsize=16)
#     ax3.tick_params(axis='both', which='major', labelsize=16)
#     # melting temperature index specific visuals
#     TM_cutoff = (0, 1)
#     urease_TMI = 0.62
#     catalase_TMI = 1.22
#     ax1.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
#                 alpha=0.2)
#     ax1.axvline(x=catalase_TMI, c='r', alpha=1.0)
#     ax1.axvline(x=urease_TMI, c='b', alpha=1.0)
#     ax2.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
#                 alpha=0.2)
#     ax2.axvline(x=catalase_TMI, c='r', alpha=1.0)
#     ax2.axvline(x=urease_TMI, c='b', alpha=1.0)
#     ax3.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
#                 alpha=0.2)
#     ax3.axvline(x=catalase_TMI, c='r', alpha=1.0)
#     ax3.axvline(x=urease_TMI, c='b', alpha=1.0)
#     ax3.set_xlabel('thermostability index', fontsize=16)
#     ax1.set_ylabel('', fontsize=16)
#     ax2.set_ylabel('count', fontsize=16)
#     ax3.set_ylabel('', fontsize=16)
#     ax1.set_xlim(-2, 3)
#     ax2.set_xlim(-2, 3)
#     ax3.set_xlim(-2, 3)
#     ax1.set_ylim(0, max1+15)
#     ax2.set_ylim(0, max2+15)
#     ax3.set_ylim(0, max3+15)
#     start, end = ax1.get_ylim()
#     ax1.set_yticks(np.arange(0, end, int(end/4 + 1)))
#     # ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
#     start, end = ax2.get_ylim()
#     ax2.set_yticks(np.arange(0, end, int(end/4 + 1)))
#     # ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
#     start, end = ax3.get_ylim()
#     ax3.set_yticks(np.arange(0, end, int(end/4 + 1)))
#     # legend
#     ax1.legend(fontsize=16)
#     ax2.legend(fontsize=16)
#     ax3.legend(fontsize=16)
#     fig.tight_layout()
#     filename = output_dir+"dist_TM_index_1fig_"
#     filename += plot_suffix+".pdf"
#     fig.savefig(filename,
#                 dpi=720, bbox_inches='tight')


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
    """Plot the synthetic accessibility of a molecules VS its no. heavy atoms.

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
        NHA = m.mol.GetNumHeavyAtoms()
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
    """Plot the synthetic accessibility of a molecules VS its no. rotatable
    bonds.

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
        NRB = Descriptors.rdMolDescriptors.CalcNumRotatableBonds(m.mol)
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
    """Plot the logP VS logS of all molecules.

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
    """Plot the logP VS XlogP (PubChem) of all molecules.

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
    """Plot distributions of molecule attributes.

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
    """Plot distribution of logP.

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
    """Plot distribution of logS.

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
    """Plot distribution of SAscore.

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


def violin_max_size(output_dir, generator, plot_suffix):
    """Do violin plots of all properties for all EC output files.

    """
    delta = {}
    # iterate over reaction system files
    for rs in generator:
        if rs.skip_rxn is True:
            continue
        if rs.max_comp_size is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta.keys()):
                delta[top_EC] = []
            delta[top_EC].append(rs.max_comp_size)

    # bin each of the sets of data based on X value
    fig, ax = plt.subplots(figsize=(8, 5))
    for keys, values in delta.items():
        if keys == '-':
            number = 0
        else:
            number = int(keys)
        parts = ax.violinplot(values, [number],
                              showmeans=False,
                              showmedians=False,
                              showextrema=False,)
        for pc in parts['bodies']:
            pc.set_facecolor(EC_descriptions()[keys][1])
            pc.set_edgecolor('black')
            pc.set_alpha(0.6)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('EC number', fontsize=16)
    ax.set_ylabel('$d$ of largest component [$\mathrm{\AA}$]', fontsize=16)
    ax.set_xlim(-1, 7)
    ax.set_ylim(2.5, 15)
    ax.set_xticks([0, 1, 2, 3, 4, 5, 6])
    ax.set_xticklabels(['unknown', '1', '2', '3', '4', '5', '6'])
    fig.tight_layout()
    fig.savefig("violin_max_size.pdf",
                dpi=720, bbox_inches='tight')


if __name__ == "__main__":
    import sys
    from ercollect.rxn_syst import yield_rxn_syst

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

    print('--------------------------------------------------------------')
    print('Plotting all the plots')
    print('--------------------------------------------------------------')
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
        rs_number_rxns_vs_size(output_dir=search_output_dir,
                               size_thresh=size_thresh,
                               generator=yield_rxn_syst(search_output_dir),
                               plot_suffix=plot_suffix)
    if input('do dist_logPs? (t/f)') == 't':
        print('doing....')
        rs_dist_logP(output_dir=search_output_dir,
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
        rs_dist_logS(output_dir=search_output_dir,
                     generator=yield_rxn_syst(search_output_dir),
                     plot_suffix=plot_suffix,
                     extreme='min')
    if input('do dist_logPs -- no EC? (t/f)') == 't':
        print('doing....')
        rs_dist_logP_noEC(output_dir=search_output_dir,
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
        rs_dist_logS_noEC(output_dir=search_output_dir,
                          generator=yield_rxn_syst(search_output_dir),
                          plot_suffix=plot_suffix,
                          extreme='min')
    # rs_dist_delta_complexity_vs_size(
    #                 output_dir=search_output_dir,
    #                 generator=yield_rxn_syst(search_output_dir),
    #                 plot_suffix=plot_suffix)
    if DB_switch == 1:
        # print new reactions
        print_new_rxns(output_dir=search_output_dir,
                       generator=yield_rxn_syst(search_output_dir))
    if input('do dist_no prod and reacts? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the number of reactnts in each reaction system
        rs_dist_no_reactants(output_dir=search_output_dir,
                             generator=yield_rxn_syst(search_output_dir),
                             plot_suffix=plot_suffix)
        # plot a distribution of the number of products in each reaction system
        rs_dist_no_products(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)
    if input('do dist max SIZE? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        rs_dist_max_size(output_dir=search_output_dir,
                         generator=yield_rxn_syst(search_output_dir),
                         plot_suffix=plot_suffix)
    if input('do dist max SIZE -- no EC? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        rs_dist_max_size_noEC(output_dir=search_output_dir,
                              generator=yield_rxn_syst(search_output_dir),
                              plot_suffix=plot_suffix)
    if input('do violin max SIZE? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        violin_max_size(output_dir=search_output_dir,
                        generator=yield_rxn_syst(search_output_dir),
                        plot_suffix=plot_suffix)
    if input('do dist_delta SIZE? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in molecule size due to reaction
        rs_dist_delta_size(output_dir=search_output_dir,
                           generator=yield_rxn_syst(search_output_dir),
                           plot_suffix=plot_suffix)
    if input('do dist_deltaSA? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in synthetic accesibility
        rs_dist_delta_SA(output_dir=search_output_dir,
                         generator=yield_rxn_syst(search_output_dir),
                         plot_suffix=plot_suffix)
    if input('do dist_deltaSA? -- no EC (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in synthetic accesibility
        rs_dist_delta_SA_noEC(output_dir=search_output_dir,
                              generator=yield_rxn_syst(search_output_dir),
                              plot_suffix=plot_suffix)

    if input('do dist_deltacomplexity? (t/f)') == 't':
        print('doing....')
        # plot a distribution of the change in synthetic accesibility
        rs_dist_delta_complexity(output_dir=search_output_dir,
                                 generator=yield_rxn_syst(search_output_dir),
                                 plot_suffix=plot_suffix)

    # plot distributions of protein sequence properties
    if DB_switch != 3:
        if input('do dist_GRAVY? (t/f)') == 't':
            print('doing....')
            rs_dist_GRAVY(output_dir=search_output_dir,
                          generator=yield_rxn_syst(search_output_dir),
                          plot_suffix=plot_suffix)
        if input('do dist_I index? (t/f)') == 't':
            print('doing....')
            rs_dist_I_index(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)
        if input('do dist_A index? (t/f)') == 't':
            print('doing....')
            rs_dist_A_index(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)

        if input('do dist_TM index? (t/f)') == 't':
            print('doing....')
            rs_dist_TM_index(output_dir=search_output_dir,
                             generator=yield_rxn_syst(search_output_dir),
                             plot_suffix=plot_suffix)
        if input('do dist_pI? (t/f)') == 't':
            print('doing....')
            rs_dist_pI(output_dir=search_output_dir,
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
