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
                # so I can manipulate that later as biopython doesn't like "__"
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


def update_seq_output(output_file, output, ROW):
    """Update sequence information output file with ROW.
    Returns associated PANDAS dataframe.

    """
    output = output.append(ROW, ignore_index=True)
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
    return output


def write_seq_output(output_file):
    """Write new sequence information output file. Returns associated PANDAS
    dataframe.

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
    nat_AA = ["G", "P", "V", "A", "L", "I", "M", "C", "F", "Y", "W", "H",
              "R", "K", "Q", "T", "D", "N", "S", "E"]
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
            (Biopython: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-pysrc.html)
        - instability index:
            (Biopython: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-pysrc.html)
        - aliphatic index:
            (code from: https://github.com/ddofer/ProFET/blob/master/ProFET/feat_extract/ProtFeat.py)
            Under GNU GPL
        - GRAVY:
            (Biopython: http://biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-pysrc.html)

    Keywords:
        output_dir (str) - directory to output reaction system files

    """
    if input('load existing data? (t/f)') == 't':
        # load existing data from this FASTA file
        if isfile(output_file) is True:
            output = read_seq_output(output_file)
        else:
            output = write_seq_output(output_file)
    else:
        # overwrite output file
        output = write_seq_output(output_file)
    print('-----------------------------------------------------------')
    print('doing calculations...')
    # need to fix the FASTA output format so BIOPYTHON can read it
    file_mod = fix_fasta(FASTA_file=fasta_file)
    total_start_time = time.time()
    total_seq = get_no_seq(FASTA_file=file_mod)
    print_opt = arange(0, total_seq, 1000)
    total_seq_done = 0
    # iterate through sequences in FASTA file
    done = list(output.acc_code)
    with open(file_mod, "r") as handle:
        generator = SeqIO.parse(handle, "fasta", alphabet=IUPAC.protein)
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
            # check sequence string for uknown or non-natural amino acid
            natural = check_sequence(sequence_string=sequence_string)
            seq_obj = ProteinAnalysis(''.join(seq))
            if natural is False:
                done.append(acc_code)
                continue
            # do calculations
            pI = seq_obj.isoelectric_point()
            GRAVY = seq_obj.gravy()
            I_index = seq_obj.instability_index()
            A_index = calculate_seq_aliphatic_index(sequence_string)
            TM_index = calculate_TM_index(seq_string=sequence_string)
            done.append(acc_code)
            ROW = pd.DataFrame({'acc_code': acc_code, 'organism': organism,
                                'EC_code': EC_code,  'species': species,
                                'note': note, 'pI': pI, 'GRAVY': GRAVY,
                                'I_index': I_index, 'A_index': A_index,
                                'TM_index': TM_index}, index=[0])
            # save to output file
            output = update_seq_output(output_file, output, ROW)
            if i in print_opt:
                print(i+1, 'done of', total_seq)
    print('--- finished %s sequences in %s seconds ---'
          % (total_seq_done, '{0:.2f}'.format(time.time() - total_start_time)))


def dist_plot(fig, ax, name, xlim, xtitle):
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


def dist_Aindex(output_file, plot_suffix):
    """Plot distribution of protein Aindex for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    output = read_seq_output(output_file)
    width = 5
    X_bins = arange(0, 150, width)
    hist, bin_edges = histogram(a=list(output.A_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')

    # AI specific visuals
    ylim = ax.get_ylim()
    ax.text(10, max(ylim)/2 + 0.05*(max(ylim)/2), 'more stable', fontsize=16)
    ax.arrow(10, max(ylim)/2, 40, 0,
             head_width=0.05*(max(ylim)/2), head_length=5, fc='k', ec='k')
    catalase_AI = 68
    ax.axvline(x=catalase_AI, c='r', alpha=1.0)
    urease_AI = 90.476
    ax.axvline(x=urease_AI, c='b', alpha=1.0)
    dist_plot(fig, ax, name='Aindex', xlim=(0, 150),
              xtitle='aliphatic index')


def dist_Iindex(output_file, plot_suffix):
    """Plot distribution of protein Aindex for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    output = read_seq_output(output_file)
    width = 5
    X_bins = arange(0, 150, width)
    hist, bin_edges = histogram(a=list(output.I_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')

    # instability specific visuals
    # get ylim
    ylim = ax.get_ylim()
    ax.text(41, max(ylim)/2 + 0.05*(max(ylim)/2), 'unstable', fontsize=16)
    ax.arrow(40, max(ylim)/2, 30, 0,
             head_width=0.05*(max(ylim)/2), head_length=4, fc='k', ec='k')
    II_cutoff = 40
    ax.axvline(x=II_cutoff, c='grey', alpha=1.0, linestyle='--')
    catalase_II = 27.010
    ax.axvline(x=catalase_II, c='r', alpha=1.0)
    urease_II = 31.75
    ax.axvline(x=urease_II, c='b', alpha=1.0)
    dist_plot(fig, ax, name='Iindex', xlim=(0, 100),
              xtitle='instability index')


def dist_TMindex(output_file, plot_suffix):
    """Plot distribution of protein Aindex for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    output = read_seq_output(output_file)
    width = 0.5
    X_bins = arange(-5, 5.1, width)
    hist, bin_edges = histogram(a=list(output.TM_index), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')

    # melting temperature index specific visuals
    TM_cutoff = (0, 1)
    ax.axvspan(xmin=TM_cutoff[0], xmax=TM_cutoff[1], facecolor='grey',
               alpha=0.2)
    catalase_TMI = 1.22
    ax.axvline(x=catalase_TMI, c='r', alpha=1.0)
    urease_TMI = 0.62
    ax.axvline(x=urease_TMI, c='b', alpha=1.0)
    dist_plot(fig, ax, name='TMindex', xlim=(-5, 5),
              xtitle='thermostability index')


def dist_pI(output_file, plot_suffix):
    """Plot distribution of protein Aindex for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    output = read_seq_output(output_file)
    width = 0.5
    X_bins = arange(0, 14.1, width)
    hist, bin_edges = histogram(a=list(output.pI), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')
    dist_plot(fig, ax, name='pI', xlim=(0, 14),
              xtitle='pI')


def dist_GRAVY(output_file, plot_suffix):
    """Plot distribution of protein GRAVY for all sequences in FASTA file.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    output = read_seq_output(output_file)
    width = 0.05
    X_bins = arange(-2, 2.2, width)
    hist, bin_edges = histogram(a=list(output.GRAVY), bins=X_bins)
    # output.GRAVY.plot.hist(bins=50,
    #                        color='#607c8e')
    ax.plot(X_bins[:-1]+width/2, hist, c='k', lw='2')

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
    dist_plot(fig, ax, name='GRAVY', xlim=(-1.5, 1.5),
              xtitle='GRAVY')


def fasta_plotting(output_file, plot_suffix):
    """Plot data for bulk sequence analysis. output_file defines file names,
    titles.

    """
    dist_GRAVY(output_file, plot_suffix)
    dist_Iindex(output_file, plot_suffix)
    dist_Aindex(output_file, plot_suffix)
    dist_TMindex(output_file, plot_suffix)
    dist_pI(output_file, plot_suffix)


def main_analysis(plot_suffix, fasta_file, output_file):
    """Analyse all sequences in FASTA file from BRENDA.

    """
    print('--------------------------------------------------------------')
    print('Analyse properties of all sequences in FASTA:', fasta_file)
    print('--------------------------------------------------------------')
    # percent_w_sequence(output_dir=search_output_dir)
    temp_time = time.time()
    if input('run calculations? (t/f)') == 't':
        get_fasta_sequence_properties(output_file=output_file,
                                      fasta_file=fasta_file)
    # do plotting + analysis -- plots
    print('-----------------------------------------------------------')
    print('doing analysis...')
    # load existing data from this FASTA file
    fasta_plotting(output_file=output_file, plot_suffix=plot_suffix)
    print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')


if __name__ == "__main__":
    import sys
    import time

    if (not len(sys.argv) == 3):
        print('Usage: bulk_protein_analysis.py FASTA plot_suffix\n')
        print('   FASTA: FASTA file containing sequences to analyse.')
        print('   plot_suffix: string to put at the end of plot file names.')
        sys.exit()
    else:
        FASTA = sys.argv[1]
        plot_suffix = sys.argv[2]

    output_file = FASTA.replace('.fasta', '_output.csv')
    print(output_file)
    main_analysis(plot_suffix=plot_suffix, fasta_file=FASTA,
                  output_file=output_file)
    sys.exit()
