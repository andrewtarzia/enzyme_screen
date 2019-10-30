#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse FASTA sequence files from BRENDA in bulk.

Author: Andrew Tarzia

Date Created: 06 Dec 2018

"""
from os.path import isfile
import pandas as pd
from numpy import arange, histogram
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import IUPAC
from ercollect.rs_protein_analysis import calculate_seq_aliphatic_index
from ercollect.tm_predictor import calculate_TM_index
from ercollect.plotting import EC_descriptions


def specific_EC_descriptions():
    """Dictionary of EC descriptions + colours for target reactions.

    """
    top_tier = {'3.5.5.1': ('nitrilase', '#1469b5'),
                '3.5.5.4': ('cyanoalanine nitrilase', '#FF7900'),
                '4.2.1.65': ('3-cyanoalanine hydratase', '#00B036'),
                }

    return top_tier


def fix_fasta(FASTA_file):
    """
    Fix FASTA files to be in BioPYTHON format.

    Arguments:
        FASTA_file (str) - FASTA file name to fix

    """
    file_mod = FASTA_file.replace(".fasta", "_mod.fasta")
    if isfile(file_mod) is False:
        with open(FASTA_file, 'r', encoding='mac-roman') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if '|' in line and ">" not in line:
                # we replace spaces in header line with "__"
                # so I can manipulate that later as biopython
                # doesn't like "__"
                new_line = ">"+line
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        with open(file_mod, 'w') as f:
            for line in new_lines:
                f.write(line)
    return file_mod


def read_seq_output(output_file):
    """Read sequence information output file. Returns associated PANDAS
    dataframe.

    """
    output = pd.read_table(output_file, delimiter='@', skiprows=[0],
                           names=['acc_code', 'organism', 'EC_code',
                                  'species', 'note', 'pI', 'GRAVY',
                                  'I_index', 'A_index', 'TM_index'],
                           engine='python')
    return output


def update_seq_output(output_file, ROW):
    """Update sequence information output file with ROW.
    Returns associated PANDAS dataframe.

    """
    # output = output.append(ROW, ignore_index=True)
    # output.to_csv(output_file, index=False, sep='@')
    with open(output_file, 'a') as f:
        string = ROW.acc_code.iloc[0]+'@'
        string += ROW.organism.iloc[0]+'@'
        string += ROW.EC_code.iloc[0]+'@'
        string += ROW.species.iloc[0]+'@'
        string += ROW.note.iloc[0]+'@'
        string += str(ROW.pI.iloc[0])+'@'
        string += str(ROW.GRAVY.iloc[0])+'@'
        string += str(ROW.I_index.iloc[0])+'@'
        string += str(ROW.A_index.iloc[0])+'@'
        string += str(ROW.TM_index.iloc[0])
        string += '\n'
        f.write(string)


def write_seq_output(output_file):
    """Write new sequence information output file. Returns associated
    PANDAS dataframe.

    """
    output = pd.DataFrame(columns=['acc_code', 'organism', 'EC_code',
                                   'species', 'note', 'pI', 'GRAVY',
                                   'I_index', 'A_index', 'TM_index'])
    output.to_csv(output_file, index=False, sep='@')
    return output


def check_sequence(sequence_string):
    """Check sequence string for unknown or non-natural amino acids.

    Returns True if only natural AA is found.

    """
    nat_AA = [
        "G", "P", "V", "A", "L", "I", "M", "C", "F", "Y", "W", "H",
        "R", "K", "Q", "T", "D", "N", "S", "E"
    ]
    for AA in sequence_string:
        if AA not in nat_AA:
            return False
    return True


def get_no_seq(FASTA_file):
    """Get number of sequences in a FASTA file from number of '>'.

    """
    seq = 0
    with open(FASTA_file, 'r') as f:
        for line in f.readlines():
            if '>' in line:
                seq += 1
    return seq


def get_fasta_sequence_properties(output_file, fasta_file):
    """Get sequence properties for all reaction systems with an associated
    protein sequence.

    Currently applied to only SABIO DB.

    Properties:
        - pI: we do not consider the possibility of modifications here.
            (Biopython: http://biopython.org/DIST/docs/api/
                Bio.SeqUtils.ProtParam-pysrc.html)
        - instability index:
            (Biopython: http://biopython.org/DIST/docs/api/
                Bio.SeqUtils.ProtParam-pysrc.html)
        - aliphatic index:
            (code from: https://github.com/ddofer/ProFET/blob
                /master/ProFET/feat_extract/ProtFeat.py)
            Under GNU GPL
        - GRAVY:
            (Biopython: http://biopython.org/DIST/docs/api/
                Bio.SeqUtils.ProtParam-pysrc.html)

    Keywords:
        output_dir (str) - directory to output reaction system files

    """
    if input('load existing data? (t/f)') == 't':
        # load existing data from this FASTA file
        if isfile(output_file) is True:
            output = read_seq_output(output_file)
        else:
            write_seq_output(output_file)
            output = read_seq_output(output_file)
    else:
        # overwrite output file
        write_seq_output(output_file)
        output = read_seq_output(output_file)
    print('-------------------------------------------------------')
    print('doing calculations...')
    # need to fix the FASTA output format so BIOPYTHON can read it
    file_mod = fix_fasta(FASTA_file=fasta_file)
    total_start_time = time.time()
    total_seq = get_no_seq(FASTA_file=file_mod)
    print_opt = arange(0, total_seq, 1000)
    total_seq_done = 0
    # iterate through sequences in FASTA file
    done = list(output.acc_code)
    del output
    with open(file_mod, "r") as handle:
        generator = SeqIO.parse(
            handle,
            "fasta",
            alphabet=IUPAC.protein
        )
        for i, record in enumerate(generator):
            total_seq_done += 1
            record_list = record.description.split("|")
            # collect information on sequence and sequence
            # should be a unique descriptor
            acc_code = record_list[0].lstrip().rstrip()
            organism = record_list[1].lstrip().rstrip()
            EC_code = record_list[2].lstrip().rstrip()
            species = record_list[3].lstrip().rstrip()
            note = record_list[4]
            if acc_code in done:
                continue
            seq = record.seq
            sequence_string = str(seq)
            # check sequence string for uknown or nonnatural amino acid
            natural = check_sequence(sequence_string=sequence_string)
            seq_obj = ProteinAnalysis(''.join(seq))
            if natural is False:
                continue
            # do calculations
            pI = seq_obj.isoelectric_point()
            GRAVY = seq_obj.gravy()
            I_index = seq_obj.instability_index()
            A_index = calculate_seq_aliphatic_index(sequence_string)
            TM_index = calculate_TM_index(seq_string=sequence_string)
            ROW = pd.DataFrame({
                'acc_code': acc_code,
                'organism': organism,
                'EC_code': EC_code,
                'species': species,
                'note': note,
                'pI': pI,
                'GRAVY': GRAVY,
                'I_index': I_index,
                'A_index': A_index,
                'TM_index': TM_index
            }, index=[0])
            # save to output file
            update_seq_output(output_file, ROW)
            if i in print_opt:
                print(
                    i+1, 'done of', total_seq,
                    'in %s seconds' % ('{0:.2f}'.format(
                        time.time() - total_start_time)
                    )
                )
    print(
        '--- finished %s sequences in %s seconds ---'
        % (
            total_seq_done,
            '{0:.2f}'.format(time.time() - total_start_time)
        )
    )


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


def dist_Aindex(output, plot_suffix, EC):
    """
    Plot distribution of protein Aindex for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 5
    X_bins = arange(0, 150, width)
    hist, bin_edges = histogram(a=list(output.A_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color=EC_descriptions()[str(EC)][1],
           edgecolor='k',
           label=EC_descriptions()[str(EC)][0])

    # AI specific visuals
    ylim = ax.get_ylim()
    ax.text(
        10, max(ylim)/2 + 0.05*(max(ylim)/2), 'more stable',
        fontsize=16
    )
    ax.arrow(
        10, max(ylim)/2, 40, 0,
        head_width=0.05*(max(ylim)/2), head_length=5, fc='k', ec='k'
    )
    # catalase_AI = 68
    # ax.axvline(x=catalase_AI, c='r', alpha=1.0)
    # urease_AI = 90.476
    # ax.axvline(x=urease_AI, c='b', alpha=1.0)
    dist_plot(fig, ax, name='Aindex', xlim=(0, 150),
              xtitle='aliphatic index', plot_suffix=plot_suffix)


def dist_Iindex(output, plot_suffix, EC):
    """
    Plot distribution of protein I index for all sequences in
    FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 5
    X_bins = arange(0, 150, width)
    hist, bin_edges = histogram(a=list(output.I_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color=EC_descriptions()[str(EC)][1],
           edgecolor='k',
           label=EC_descriptions()[str(EC)][0])

    # instability specific visuals
    # get ylim
    ylim = ax.get_ylim()
    ax.text(
        51, max(ylim)/2 + 0.05*(max(ylim)/2), 'unstable', fontsize=16
    )
    ax.arrow(
        50, max(ylim)/2, 30, 0,
        head_width=0.05*(max(ylim)/2), head_length=4, fc='k', ec='k'
    )
    II_cutoff = 40
    ax.axvline(x=II_cutoff, c='k', alpha=1.0, linestyle='--', lw=2)
    # catalase_II = 27.010
    # ax.axvline(x=catalase_II, c='r', alpha=1.0)
    # urease_II = 31.75
    # ax.axvline(x=urease_II, c='b', alpha=1.0)
    dist_plot(fig, ax, name='Iindex', xlim=(0, 100),
              xtitle='instability index', plot_suffix=plot_suffix)


def dist_TMindex(output, plot_suffix, EC):
    """
    Plot distribution of protein TM index for all sequences in
    FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.2
    X_bins = arange(-5, 5.1, width)
    hist, bin_edges = histogram(a=list(output.TM_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color=EC_descriptions()[str(EC)][1],
           edgecolor='k',
           label=EC_descriptions()[str(EC)][0])

    # melting temperature index specific visuals
    TM_cutoff = (0, 1)
    ax.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
               alpha=0.2)
    # catalase_TMI = 1.22
    # ax.axvline(x=catalase_TMI, c='r', alpha=1.0)
    # urease_TMI = 0.62
    # ax.axvline(x=urease_TMI, c='b', alpha=1.0)
    dist_plot(fig, ax, name='TMindex', xlim=(-5, 5),
              xtitle='thermostability index', plot_suffix=plot_suffix)


def dist_pI(output, plot_suffix, EC):
    """
    Plot distribution of protein pI for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.5
    X_bins = arange(0, 14.1, width)
    hist, bin_edges = histogram(a=list(output.pI), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color=EC_descriptions()[str(EC)][1],
           edgecolor='k',
           label=EC_descriptions()[str(EC)][0])
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    dist_plot(fig, ax, name='pI', xlim=(0, 14),
              xtitle='pI', plot_suffix=plot_suffix)


def dist_GRAVY(output, plot_suffix, EC):
    """
    Plot distribution of protein GRAVY for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.05
    X_bins = arange(-2, 2.2, width)
    hist, bin_edges = histogram(a=list(output.GRAVY), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color=EC_descriptions()[str(EC)][1],
           edgecolor='k',
           label=EC_descriptions()[str(EC)][0])
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')

    # GRAVY specific visuals
    # ax.text(-1.45, 40, 'hydrophilic', fontsize=16)
    # get ylim
    ylim = ax.get_ylim()
    ax.text(
        0.55, max(ylim)/2 + 0.05*(max(ylim)/2),
        'hydrophobic', fontsize=16
    )
    ax.arrow(
        0.5, max(ylim)/2, 0.7, 0,
        head_width=0.05*(max(ylim)/2), head_length=0.1, fc='k', ec='k'
    )
    # avg_GRAVY = -0.4
    # ax.axvline(x=avg_GRAVY, c='grey', alpha=1.0, linestyle='--')
    # catalase_GRAVY = -0.605
    # ax.axvline(x=catalase_GRAVY, c='r', alpha=1.0)
    # urease_GRAVY = -0.1524
    # ax.axvline(x=urease_GRAVY, c='b', alpha=1.0)
    dist_plot(fig, ax, name='GRAVY', xlim=(-1.5, 1.5),
              xtitle='GRAVY', plot_suffix=plot_suffix)


def all_EC_violin_plot():
    """Do violin plots of all properties for all EC output files.

    """
    properties = ['I_index', 'A_index', 'TM_index', 'pI', 'GRAVY']
    prop_label = ['instability index', 'aliphatic index',
                  'TM index', 'pI', 'GRAVY']
    prop_lim = [(0, 100), (0, 150), (-5, 5), (0, 14), (-1.5, 1.5)]
    ECs = ['1', '2', '3', '4', '5', '6']
    output_files = [i+'__BRENDA_sequences_output.csv' for i in ECs]

    for i, prop in enumerate(properties):
        print('doing', prop, '....')
        fig, ax = plt.subplots(figsize=(8, 5))
        for out_file in output_files:
            print(out_file)
            EC = out_file[0]
            print(EC)
            output = read_seq_output(out_file)
            parts = ax.violinplot(output[prop], [int(EC)],
                                  showmeans=False,
                                  showmedians=False,
                                  showextrema=False,)
            for pc in parts['bodies']:
                pc.set_facecolor(EC_descriptions()[EC][1])
                pc.set_edgecolor('black')
                pc.set_alpha(0.6)
        if prop == 'TM_index':
            # melting temperature index specific visuals
            TM_cutoff = (0, 1)
            ax.axhspan(ymin=TM_cutoff[0], ymax=TM_cutoff[1],
                       facecolor='grey',
                       alpha=0.2)
        if prop == 'I_index':
            II_cutoff = 40
            ax.axhline(
                y=II_cutoff, c='k', alpha=1.0, linestyle='--', lw=2
            )
        if prop == 'A_index':
            ax.text(0.21, 60,
                    'more stable', fontsize=16, ha='left', va='bottom',
                    rotation=90)
            ax.arrow(0.5, 40, 0, 80,
                     head_width=0.2, head_length=10, fc='k',
                     ec='k')
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('EC number', fontsize=16)
        ax.set_ylabel(prop_label[i], fontsize=16)
        ax.set_xlim(0, 7)
        ax.set_ylim(prop_lim[i])
        ax.set_xticks([1, 2, 3, 4, 5, 6])
        ax.set_xticklabels(['1', '2', '3', '4', '5', '6'])
        fig.tight_layout()
        fig.savefig("violin_"+prop+".pdf",
                    dpi=720, bbox_inches='tight')


def fasta_plotting(output_file, plot_suffix, EC):
    """
    Plot data for bulk sequence analysis. output_file defines
    file names, titles.

    """
    output = read_seq_output(output_file)
    print('doing GRAVY...')
    dist_GRAVY(output, plot_suffix, EC)
    print('doing I index...')
    dist_Iindex(output, plot_suffix, EC)
    print('doing A index...')
    dist_Aindex(output, plot_suffix, EC)
    print('doing TM index...')
    dist_TMindex(output, plot_suffix, EC)
    print('doing pI...')
    dist_pI(output, plot_suffix, EC)


def dist_TMindex_specific(output, plot_suffix, EC):
    """
    Plot distribution of protein TM index for all sequences in
    FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.2
    X_bins = arange(-5, 5.1, width)
    hist, bin_edges = histogram(a=list(output.TM_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    # ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.4, width=width,
           color=specific_EC_descriptions()[str(EC)][1],
           edgecolor='k',
           label=specific_EC_descriptions()[str(EC)][0])

    # melting temperature index specific visuals
    TM_cutoff = (0, 1)
    ax.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
               alpha=0.2)
    # catalase_TMI = 1.22
    # ax.axvline(x=catalase_TMI, c='r', alpha=1.0)
    # urease_TMI = 0.62
    # ax.axvline(x=urease_TMI, c='b', alpha=1.0)
    ax.set_title('EC '+EC+': '+specific_EC_descriptions()[str(EC)][0],
                 fontsize=16)
    dist_plot(fig, ax, name='TMindex', xlim=(-5, 5),
              xtitle='thermostability index', plot_suffix=plot_suffix)


def main_analysis(plot_suffix, fasta_file, output_file, EC):
    """Analyse all sequences in FASTA file from BRENDA.

    """
    print('---------------------------------------------------------')
    print('Analyse properties of all sequences in FASTA:', fasta_file)
    print('---------------------------------------------------------')
    # percent_w_sequence(output_dir=search_output_dir)
    temp_time = time.time()
    if input('run calculations? (t/f)') == 't':
        get_fasta_sequence_properties(output_file=output_file,
                                      fasta_file=fasta_file)
    # do plotting + analysis -- plots
    if EC in list(EC_descriptions().keys()):
        print('------------------------------------------------------')
        print('doing analysis...')
        # load existing data from this FASTA file
        fasta_plotting(
            output_file=output_file,
            plot_suffix=plot_suffix,
            EC=EC
        )
        print(
            '--- time taken =',
            '{0:.2f}'.format(time.time()-temp_time),
            's'
        )
    else:
        print('------------------------------------------------------')
        print('doing specific sequence analysis...')
        output = read_seq_output(output_file)
        dist_TMindex_specific(output, plot_suffix, EC)


if __name__ == "__main__":
    import sys
    import time

    if (not len(sys.argv) == 4):
        print('Usage: bulk_protein_analysis.py FASTA plot_suffix\n')
        print('   FASTA: FASTA file containing sequences to analyse.')
        print(
            '   plot_suffix: file extension.'
        )
        print('   EC: EC number.')
        sys.exit()
    else:
        FASTA = sys.argv[1]
        plot_suffix = sys.argv[2]
        EC = sys.argv[3]

    if input('do all EC plot? (t/f)') == 't':
        all_EC_violin_plot()
    if input('do single EC analysis? (t/f)') == 't':
        output_file = FASTA.replace('.fasta', '_output.csv')
        print(output_file)
        main_analysis(plot_suffix=plot_suffix, fasta_file=FASTA,
                      output_file=output_file, EC=EC)
    sys.exit()
