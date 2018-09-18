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
import matplotlib.pyplot as plt
import os
from rxn_syst import yield_rxn_syst
from rdkit.Chem import Descriptors


def define_diff_categ_plot(ax, title, ytitle, xtitle, xlim, ylim):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_ylabel(ytitle, fontsize=16)
    # ax.legend([y, n], ['aligned', 'not aligned'], loc=4, fancybox=True)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticklabels(['diffuses', 'does not diffuse'])
    ax.set_xticks([0.25, 0.75])


def define_standard_plot(ax, title, ytitle, xtitle, xlim, ylim):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    # ax.legend([y, n], ['aligned', 'not aligned'], loc=4, fancybox=True)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


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


def categorical(molecules, threshold, output_dir):
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
                        ytitle='intermediate diameter [$\mathrm{\AA}$]',
                        xlim=(0, 1),
                        ylim=(0, 10))
    fig.tight_layout()
    fig.savefig(output_dir+"categorical.pdf", dpi=720,
                bbox_inches='tight')


def categorical_moloutput(mol_output_file, threshold, output_dir):
    """Categorical scatter plot of all molecules molecule_output file.

    """
    import DB_functions
    molecule_output = DB_functions.initialize_mol_output_DF(mol_output_file,
                                                            overwrite=False)
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
    fig.savefig(output_dir+"categorical.pdf", dpi=720,
                bbox_inches='tight')


def shapes(molecules, threshold, output_dir):
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
    fig.savefig(output_dir+"shape.pdf", dpi=720,
                bbox_inches='tight')


def rs_pI_distribution(output_dir, cutoff_pI):
    """Plot distribution of pIs for all reaction systems.

    """
    import Uniprot_IO
    import pi_fn
    fig, ax = plt.subplots(figsize=(8, 5))
    native_pi = []
    succ_pi = []
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        # collect pIs of all sequences even if reaction is skipped elsewhere
        if rs.skip_rxn is True and rs.UniprotID != '' and rs.UniprotID is not None:
            try:
                pI = rs.pI
            except AttributeError:
                print('calculating pI...')
                IDs = rs.UniprotID.split(" ")
                total_sequence = ''
                for i in IDs:
                    sequence = Uniprot_IO.get_sequence(i)
                    total_sequence += sequence
                if len(total_sequence) > 0:
                    rs = pi_fn.calculate_rxn_syst_pI(total_sequence, rs,
                                                     cutoff_pi=cutoff_pI)
                else:
                    rs.pI = None
                rs.save_object(output_dir+rs.pkl)

        if rs.UniprotID != '' and rs.UniprotID is not None and rs.pI is not None:
            if rs.req_mod is None:
                native_pi.append(rs.pI)
            else:
                succ_pi.append(rs.pI)

    ax.hist(native_pi,
            facecolor='k',
            alpha=0.5,
            histtype='stepfilled',
            bins=np.arange(0, 14 + 0.2, 0.5),
            label='unmodified')

    ax.hist(succ_pi,
            facecolor='firebrick',
            alpha=0.5,
            histtype='stepfilled',
            bins=np.arange(0, 14 + 0.2, 0.5),
            label='succinylated')

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('calculated pI', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 14)
    # plot pI cut-off
    ax.axvline(x=cutoff_pI, c='k', lw='2', linestyle='--')
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"pI_dist.pdf",
                dpi=720, bbox_inches='tight')


def rs_size_vs_pI(output_dir, cutoff_pI, size_thresh):
    """Plot maximum component size of a reaction vs. pI.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        if rs.seed_MOF is None:
            continue
        if rs.all_fit is True and rs.seed_MOF is True:
            M = 'o'
            if rs.req_mod is not None:
                C = 'orange'
            else:
                C = 'b'
            E = 'k'
        else:
            M = 'o'
            C = 'r'
            E = 'k'

        ax.scatter(rs.pI,
                   rs.max_comp_size, c=C,
                   edgecolors=E, marker=M, alpha=1.0,
                   s=100)

    # decoy legend
    ax.scatter(-100, 100,
               c='b',
               edgecolors=E,
               marker='o',
               alpha=1.0,
               s=100,
               label='candidate - native')
    ax.scatter(-100, 100,
               c='orange',
               edgecolors=E,
               marker='o',
               alpha=1,
               s=100,
               label='candidate - modified')
    ax.scatter(-100, 100,
               c='r',
               edgecolors=E,
               marker='o',
               alpha=1,
               s=100,
               label='non-candidate')

    ax.legend(loc=1, fontsize=12)

    ax.axhline(y=size_thresh, c='k')
    ax.axvline(x=cutoff_pI, c='k')
    define_standard_plot(ax,
                         title='',
                         xtitle='pI',
                         ytitle='diameter of largest component [$\mathrm{\AA}$]',
                         xlim=(0, 14),
                         ylim=(0, 10))
    fig.tight_layout()
    fig.savefig(output_dir+"size_vs_pI.pdf", dpi=720,
                bbox_inches='tight')


def check_rxn_unique(reaction_reported, rs):
    """Check (using the sorted list of component molecule weights) if a rxn is
    unique.

    """
    # get list of SMILES of all components
    r_smiles = []
    r_MW = []
    for r in rs.components:
        r_smiles.append(r.SMILES)
        r_MW.append(Descriptors.MolWt(r.mol))
    r_smiles = [x for _, x in sorted(zip(r_MW, r_smiles))]
    if r_smiles in reaction_reported:
        unique = False
    else:
        reaction_reported.append(r_smiles)
        unique = True

    return unique, reaction_reported


def number_rxns_vs_size(output_dir, size_thresh, pI_thresh):
    """Plot number of possible and unique reactions as a function of size
    threshold.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    max_sizes = []
    reaction_reported = []
    max_sizes_pI = []  # also plot the number of new reactions with pI < thresh
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        unique, reaction_reported = check_rxn_unique(reaction_reported, rs)
        if unique is False:
            continue
        try:
            if rs.max_comp_size > 0:
                max_sizes.append(rs.max_comp_size)
                if rs.seed_MOF is True:
                    max_sizes_pI.append(rs.max_comp_size)
        except AttributeError:
            pass

    max_sizes = np.asarray(max_sizes)
    max_sizes_pI = np.asarray(max_sizes_pI)
    counts = []
    counts_pI = []
    threshs = np.arange(0.1, 21, 0.5)
    for thr in threshs:
        count_above = len(max_sizes[max_sizes < thr])
        counts.append(count_above)
        count_above_pI = len(max_sizes_pI[max_sizes_pI < thr])
        counts_pI.append(count_above_pI)

    ax.bar(threshs, counts, align='center', alpha=0.5, width=0.2,
           label='max component < threshold',
           color='b', edgecolor='k')
    ax.bar(threshs, counts_pI, align='center', alpha=0.5, width=0.2,
           label='+ pI < '+str(pI_thresh),
           color='r', edgecolor='k')

    ax.legend(loc=2, fontsize=12)

    ax.axvline(x=size_thresh, c='k')

    define_standard_plot(ax,
                         title='',
                         xtitle='diffusion threshold [$\mathrm{\AA}$]',
                         ytitle='# reactions',
                         xlim=(0, 11),
                         ylim=(0, max(counts)+10))
    fig.tight_layout()
    fig.savefig(output_dir+"size_vs_threshold.pdf", dpi=720,
                bbox_inches='tight')


def print_new_rxns(output_dir):
    """Print all new possible and unique reactions that fit.

    """
    reaction_reported = []
    count = 0
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
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


def rs_delta_size(output_dir):
    """Plot change in maximum size of reactants to products.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta_data = []
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
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
            # ax.scatter(rs.pI,
            #            rs.max_comp_size, c=C,
            #            edgecolors=E, marker=M, alpha=1.0,
            #            s=100)
            delta_data.append(delta_size)

    ax.hist(delta_data,
            facecolor='purple',
            alpha=1.0,
            histtype='stepfilled',
            bins=np.arange(-10, 10.2, 0.5))

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('maximum product size $-$ maximum reactant size [$\mathrm{\AA}$]',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(-10, 10)
    # legend
    # ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"delta_size_dist.pdf",
                dpi=720, bbox_inches='tight')
