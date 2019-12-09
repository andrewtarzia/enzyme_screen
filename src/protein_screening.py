#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script for plotting data that is used for screening.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import os
import sys

from reaction import yield_rxn_syst
import plots_reaction as pr
import utilities




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


def main():
    if (not len(sys.argv) == 3):
        print('Usage: screening.py plot_suffix target_EC\n')
        print("""
        plot_suffix: string to put at the end of plot file names.
        target_EC: file containing a target EC number (f if all EC
            should be used)
        param_file:

        """)
        sys.exit()
    else:
        plot_suffix = sys.argv[1]
        target_EC_file = sys.argv[2] if sys.argv[2] != 'f' else None
        pars = utilities.read_params(sys.argv[3])

    if target_EC_file is None:
        # Set options and get all EC numbers.
        settings = {}
        search_EC_file = pars['EC_file']
        search_ECs = utilities.get_ECs_from_file(
            EC_file=search_EC_file
        )
    else:
        settings = {}
        # Read ECs from file.
        search_ECs = utilities.get_ECs_from_file(
            EC_file=target_EC_file
        )

    print(settings)
    print(search_ECs)

    # plot distributions of protein sequence properties
    if DB_switch != 3:
        if input('do dist_GRAVY? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_GRAVY(output_dir=search_output_dir,
                          generator=yield_rxn_syst(search_output_dir),
                          plot_suffix=plot_suffix)
        if input('do dist_I index? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_I_index(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)
        if input('do dist_A index? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_A_index(output_dir=search_output_dir,
                            generator=yield_rxn_syst(search_output_dir),
                            plot_suffix=plot_suffix)

        if input('do dist_TM index? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_TM_index(output_dir=search_output_dir,
                             generator=yield_rxn_syst(search_output_dir),
                             plot_suffix=plot_suffix)
        if input('do dist_pI? (t/f)') == 't':
            print('doing....')
            pfn.rs_dist_pI(output_dir=search_output_dir,
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


if __name__ == "__main__":
    main()
