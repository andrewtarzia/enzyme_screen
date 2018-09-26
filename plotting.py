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
import matplotlib.colors as colors
import matplotlib.cm as cm


def EC_descriptions():
    """Dictionary of EC descriptions + colours.

    """
    top_tier = {'1': ('oxidoreductases', 'k'),
                '2': ('transferases', 'b'),
                '3': ('hydrolases', 'orange'),
                '4': ('lyases', 'g'),
                '5': ('isomerases', 'purple'),
                '6': ('ligases', 'r')}

    return top_tier


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


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero


    From Stack Exchange:
        https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


def define_plot_cmap(fig, ax, mid_point, cmap, ticks, labels, cmap_label):
    """Define cmap shifted to midpoint and plot colourbar

    """
    new_cmap = shiftedColorMap(cmap, start=0, midpoint=mid_point,
                               stop=1, name='shifted')

    X = np.linspace(0, 1, 256)

    cax = ax.scatter(-X-100, -X-100, c=X, cmap=new_cmap)

    cbar = fig.colorbar(cax, ticks=ticks,
                        spacing='proportional')
    cbar.ax.set_yticklabels(labels,
                            fontsize=16)
    cbar.set_label(cmap_label, fontsize=16)

    return new_cmap


def rs_size_vs_SA_vs_logP(output_dir, size_thresh):
    """Plot maximum component size of a reaction vs. logP with synth access as
    3rd variable.

    Plot max_logP => most hydroXX component.

    """
    HRP_logP = 0.017
    fig, ax = plt.subplots(figsize=(8, 5))
    new_cmap = define_plot_cmap(fig, ax, mid_point=0.5, cmap=cm.RdBu,
                                ticks=[0, 0.25, 0.5, 0.75, 1],
                                labels=['-10', '-5', '0', '5', '10'],
                                cmap_label='$\Delta$ synthetic accessibility')
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        M = 'o'
        E = 'k'

        ax.scatter(rs.max_logP,
                   rs.max_comp_size,
                   c=new_cmap(abs((-10-rs.delta_sa))/20),
                   edgecolors=E,
                   marker=M,
                   alpha=1.0,
                   s=100)

    ax.axhline(y=size_thresh, c='gray', linestyle='--')
    ax.axvline(x=HRP_logP, c='gray', linestyle='--')
    define_standard_plot(ax,
                         title='',
                         xtitle='logP of most hydrophobic component',
                         ytitle='diameter of largest component [$\mathrm{\AA}$]',
                         xlim=(-6, 6),
                         ylim=(0, 15))
    fig.tight_layout()
    fig.savefig(output_dir+"size_vs_SA_vs_logP.pdf", dpi=720,
                bbox_inches='tight')


def rs_size_vs_complexity_vs_XlogP(output_dir, size_thresh):
    """Plot maximum component size of a reaction vs. XlogP with complexity as
    3rd variable.

    Plot max_logP => most hydrophobic component.

    """
    HRP_XlogP = -0.9
    fig, ax = plt.subplots(figsize=(8, 5))
    new_cmap = define_plot_cmap(fig, ax, mid_point=0.5, cmap=cm.RdBu,
                                ticks=[0, 0.5, 1],
                                labels=['-300', '0', '300'],
                                cmap_label='$\Delta$ complexity')
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        M = 'o'
        E = 'k'

        ax.scatter(rs.max_XlogP,
                   rs.max_comp_size,
                   c=new_cmap(abs((-300-rs.delta_comp))/600),
                   edgecolors=E,
                   marker=M,
                   alpha=1.0,
                   s=100)

    ax.axhline(y=size_thresh, c='gray', linestyle='--')
    ax.axvline(x=HRP_XlogP, c='gray', linestyle='--')
    define_standard_plot(ax,
                         title='',
                         xtitle='XlogP of most hydrophobic component',
                         ytitle='diameter of largest component [$\mathrm{\AA}$]',
                         xlim=(-6, 6),
                         ylim=(0, 15))
    fig.tight_layout()
    fig.savefig(output_dir+"size_vs_complexity_vs_XlogP.pdf", dpi=720,
                bbox_inches='tight')


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


def rs_number_rxns_vs_size(output_dir, size_thresh):
    """Plot number of possible and unique reactions as a function of size
    threshold.

    pI count commented out.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    max_sizes = []
    reaction_reported = []
    # also plot the number of new reactions with pI < thresh
    # max_sizes_pI = []
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
                # if rs.seed_MOF is True:
                #     max_sizes_pI.append(rs.max_comp_size)
        except AttributeError:
            pass

    max_sizes = np.asarray(max_sizes)
    # max_sizes_pI = np.asarray(max_sizes_pI)
    counts = []
    # counts_pI = []
    threshs = np.arange(0.1, 21, 0.5)
    for thr in threshs:
        count_above = len(max_sizes[max_sizes < thr])
        counts.append(count_above)
        # count_above_pI = len(max_sizes_pI[max_sizes_pI < thr])
        # counts_pI.append(count_above_pI)

    ax.plot(threshs, counts, alpha=1.0,
            label='max component < threshold',
            color='k', marker='o')
    # ax.bar(threshs, counts, align='center', alpha=0.5, width=0.2,
    #        label='max component < threshold',
    #        color='b', edgecolor='k')
    # ax.bar(threshs, counts_pI, align='center', alpha=0.5, width=0.2,
    #        label='+ pI < '+str(pI_thresh),
    #        color='r', edgecolor='k')

    ax.legend(loc=4, fontsize=12)

    ax.axvline(x=size_thresh, c='gray', linestyle='--')

    # plot possible region of ZIF pore limiting diameters from
    # Materials Project
    ax.axvspan(3.7, 16, facecolor='#2ca02c', alpha=0.2)

    define_standard_plot(ax,
                         title='',
                         xtitle='diffusion threshold [$\mathrm{\AA}$]',
                         ytitle='# reactions',
                         xlim=(0, 20),
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
    delta_data = {}
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
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta_data.keys()):
                delta_data[top_EC] = []
            delta_data[top_EC].append(delta_size)

    for keys, values in delta_data.items():
        ax.hist(values,
                facecolor=EC_descriptions()[keys][1],
                alpha=0.4,
                histtype='stepfilled',
                bins=np.arange(-10, 10.2, 0.5),
                label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('maximum product size $-$ maximum reactant size [$\mathrm{\AA}$]',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(-10, 10)
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"delta_size_dist.pdf",
                dpi=720, bbox_inches='tight')


def rs_no_reactants(output_dir):
    """Plot distribution of the number of reactants in all reactions.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    no_reacts = {}
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
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
    ax.set_xlabel('no. reactants',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 5)
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"no_reacts_dist.pdf",
                dpi=720, bbox_inches='tight')


def rs_no_products(output_dir):
    """Plot distribution of the number of products in all reactions.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    no_reacts = {}
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
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
    ax.set_xlabel('no. products',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 5)
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"no_prods_dist.pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_deltaSA(output_dir):
    """Plot distribution of the change in synthetic accesibility from react to
    products.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta = {}
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        if rs.delta_sa is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta.keys()):
                delta[top_EC] = []
            delta[top_EC].append(rs.delta_sa)

    for keys, values in delta.items():
        ax.hist(values,
                facecolor=EC_descriptions()[keys][1],
                alpha=0.4,
                histtype='stepfilled',
                bins=np.arange(-10, 10.2, 0.5),
                label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('$\Delta$ synthetic accessibility',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(-10, 10)
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"delta_SA_dist.pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_complexity(output_dir):
    """Plot distribution of molecule complexity in all rxns.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    dists = {}
    dists_mols = {}
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        for m in rs.components:
            if m.complexity is not None:
                top_EC = rs.EC.split('.')[0]
                if top_EC not in list(dists.keys()):
                    dists[top_EC] = []
                    dists_mols[top_EC] = []
                if m.SMILES not in dists_mols[top_EC]:
                    dists[top_EC].append(m.complexity)
                    dists_mols[top_EC].append(m.SMILES)

    for keys, values in dists.items():
        ax.hist(values,
                facecolor=EC_descriptions()[keys][1],
                alpha=0.4,
                histtype='stepfilled',
                bins=np.arange(0, 500, 10),
                label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('complexity',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 500)
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"complexity_dist.pdf",
                dpi=720, bbox_inches='tight')


def rs_dist_deltacomplexity(output_dir):
    """Plot distribution of the change in complexity from react to
    products.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    delta = {}
    # iterate over reaction system files
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is True:
            continue
        if rs.delta_comp is not None:
            top_EC = rs.EC.split('.')[0]
            if top_EC not in list(delta.keys()):
                delta[top_EC] = []
            delta[top_EC].append(rs.delta_comp)

    for keys, values in delta.items():
        ax.hist(values,
                facecolor=EC_descriptions()[keys][1],
                alpha=0.4,
                histtype='stepfilled',
                bins=np.arange(-500, 500, 10),
                label=EC_descriptions()[keys][0])

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('$\Delta$ complexity',
                  fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(-300, 300)
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+"delta_complexity_dist.pdf",
                dpi=720, bbox_inches='tight')

 # "fig1, ax1 = plt.subplots(figsize=(8,5))\n",
 # "fig2, ax2 = plt.subplots(figsize=(8,5))\n",
 # "\n",
 # "for idx, row in molecule_output.iterrows():\n",
 # "    name = row['name']\n",
 # "    smiles = row['SMILE']\n",
 # "    if row['mid_diam'] == 0:\n",
 # "        continue\n",
 # "    mol = Chem.MolFromSmiles(smiles)\n",
 # "    mol = Chem.AddHs(mol)\n",
 # "    #print('-------------------------------------')\n",
 # "    #print('molecule:', name)\n",
 # "    logP = Descriptors.MolLogP(mol, includeHs=True)\n",
 # "    #print('logP =', logP)\n",
 # "    SA = molecule.get_SynthA_score(mol)\n",
 # "    #print('synthetic accesibility =', SA)\n",
 # "    # print(name,'__',logP,'__',SA)\n",
 # "    #print('-------------------------------------')\n",
 # "    if logP <= 0:\n",
 # "        C = 'b'\n",
 # "    else:\n",
 # "        C = 'r'\n",
 # "    ax1.scatter(logP,\n",
 # "                row['mid_diam'], c=C, \n",
 # "                edgecolors='k', marker='o', alpha=1.0,\n",
 # "                s=100)\n",
 # "    ax2.scatter(SA,\n",
 # "                row['mid_diam'], c='purple', \n",
 # "                edgecolors='k', marker='o', alpha=1.0,\n",
 # "                s=100)\n",
 # "    \n",
 # "\n",
 # "ax1.tick_params(axis='both', which='major', labelsize=16)\n",
 # "ax1.set_xlabel('logP', fontsize=16)\n",
 # "ax1.set_ylabel('intermediate diameter [$\\mathrm{\\AA}$]', fontsize=16)\n",
 # "ax1.set_xlim(-10, 10)\n",
 # "ax1.set_ylim(0, 15)\n",
 # "    \n",
 # "ax2.tick_params(axis='both', which='major', labelsize=16)\n",
 # "ax2.set_xlabel('synthetic accessibility', fontsize=16)\n",
 # "ax2.set_ylabel('intermediate diameter [$\\mathrm{\\AA}$]', fontsize=16)\n",
 # "ax2.set_xlim(0, 10)\n",
 # "ax2.set_ylim(0, 15)\n",
 # "    \n",
 # "fig1.tight_layout()\n",
 # "fig1.savefig(output_dir+\"logP_scatter.pdf\",\n",
 # "             dpi=720, bbox_inches='tight')\n",
 # "fig2.tight_layout()\n",
 # "fig2.savefig(output_dir+\"SA_scatter.pdf\",\n",
 # "             dpi=720, bbox_inches='tight')"
