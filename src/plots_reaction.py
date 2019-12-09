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
    for keys in delta_data:
        values = delta_data[keys]
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
    for keys in data:
        values = data[keys]
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
    for keys in data:
        values = data[keys]
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

    for keys in no_reacts:
        values = no_reacts[keys]
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

    for keys in no_reacts:
        values = no_reacts[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
    for keys in delta:
        values = delta[keys]
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
