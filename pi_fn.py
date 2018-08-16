#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for pI calculation program.

Author: Andrew Tarzia

Date Created: 24 Apr 2018

License:


"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import IUPAC


def fix_fasta(database_names):
    """
    Fix FASTA files to be in BioPYTHON format.

    Arguments:
        database_names (list) - list of FASTA file names

    """
    for file in database_names:
        file_mod = file.replace(".fasta", "_mod.fasta")
        with open(file, 'r') as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if '|' in line and ">" not in line:
                # we replace spaces in header line with "__"
                # so I can manipulate that later as biopython doesn't like "__"
                new_line = ">"+line.replace(" ", "__")
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        with open(file_mod, 'w') as f:
            for line in new_lines:
                f.write(line)


def convert_to_one_letter_code_sing(seq):
    """
    Converts a single amino acid one letter code to the 3 letter code.

    Arguments:
        -seq (str) - one letter AA code.

    """
    conversion = {"GLY": "G", "PRO": "P", "VAL": "V", "ALA": "A", "LEU": "L",
                  "ILE": "I", "MET": "M", "CYS": "C", "PHE": "F", "TYR": "Y",
                  "TRP": "W", "HIS": "H", "ARG": "R", "LYS": "K", "GLN": "Q",
                  "THR": "T", "ASP": "D", "ASN": "N", "SER": "S", "GLU": "E"}
    n_seq = conversion[seq]
    return n_seq


def plot_pi_data(param_dict):
    """
    Plot and save histogram of all pI data calculated in this run.

    Arguments:
        param_dict (dict) - dictionary of parameters.

    """
    # read in CSV
    pi_data = pd.read_csv(param_dict['out_CSV_pi'], index_col=False)

    plt.rcParams['font.size'] = 12

    fig, ax = plt.subplots(figsize=(8, 5))
    # unmodifed
    mod_dict = param_dict['modifications']['0']
    data = pi_data[pi_data['modification'] == 0]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # modification 1 - succinylation
    mod_dict = param_dict['modifications']['1']
    data = pi_data[pi_data['modification'] == 1]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    ax.set_xlabel('calculated pI')
    ax.set_ylabel('count [arb. units]')
    ax.set_yticklabels([])

    # plot pI cut-off
    ax.axvline(x=param_dict['cutoff_pi'], c='k', lw='2', linestyle='--')

    ax.set_xlim(0, 14)
    ax.legend()
    fig.savefig('pI_histogram.pdf', dpi=720, bbox_inches='tight')


def plot_EC_pI_dist(EC_pi_data, param_dict, filename, title):
    """
    Plot and save histogram of all pI data calculated in this run.

    Arguments:
        param_dict (dict) - dictionary of parameters.

    """
    fig, ax = plt.subplots()

    # unmodifed
    mod_dict = param_dict['modifications']['0']
    data = EC_pi_data[EC_pi_data['modification'] == 0]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # modification 1 - succinylation
    mod_dict = param_dict['modifications']['1']
    data = EC_pi_data[EC_pi_data['modification'] == 1]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('calculated pI', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 14)
    # plot pI cut-off
    ax.axvline(x=param_dict['cutoff_pi'], c='k', lw='2', linestyle='--')
    # legend
    ax.legend(fontsize=16)
    # title
    ax.set_title('EC = '+title, fontsize=16)

    fig.tight_layout()
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')


def calculate_pI_from_file(file, param_dict, output_dir):
    """Calculate the pI of all sequences in FASTA file.

    """
    modifications = param_dict['modifications']
    count_sequences_done = 0
    total_start_time = time.time()
    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta", alphabet=IUPAC.protein):
            record_list = record.description.split("|")
            # get meta data
            acc_code, organism, EC_code, species, note = get_record_meta(record_list)
            # get unmodified pI
            seq_obj = ProteinAnalysis(''.join(record.seq))
            pi = seq_obj.isoelectric_point()
            count_sequences_done += 1
            modifier = '0'
            if pi < param_dict['cutoff_pi']:
                category = '0'
            else:
                category = '1'
            # output to CSV
            output_pI_row(output_dir, param_dict, file,
                          acc_code, organism, EC_code,
                          species, note,
                          pi, modifier, category)

            # if the category is 1 - i.e. pi > cutoff
            # then we test modification
            if category == '1':
                modifier = '1'
                # get modified pI
                seq = record.seq
                # replace target amino acid residue
                # with replacement amino acid residue
                # one letter codes
                targ = convert_to_one_letter_code_sing(modifications[modifier]['target_res'])
                replacement = convert_to_one_letter_code_sing(modifications[modifier]['replace_res'])
                mod_seq = ''.join(seq).replace(targ, replacement)
                seq_obj = ProteinAnalysis(mod_seq)
                pi = seq_obj.isoelectric_point()
                count_sequences_done += 1
                if pi < param_dict['cutoff_pi']:
                    category = '0'
                else:
                    category = '1'
                # output to CSV
                output_pI_row(output_dir, param_dict, file,
                              acc_code, organism, EC_code,
                              species, note,
                              pi, modifier, category)
            # break
    print('--- finished %s sequences in %s seconds ---'
          % (count_sequences_done, '{0:.2f}'.format(time.time() - total_start_time)))


def output_pI_row(output_dir, param_dict, file,
                  acc_code, organism, EC_code, species, note,
                  pi, modifier, category):
    """Output results of pI calculation to CSV.

    """
    with open(output_dir+param_dict['out_CSV_pi'], 'a') as f:
        string = file+','
        string += acc_code+','
        string += organism.replace(',', '_')+','
        string += EC_code+','
        string += species.replace(',', '_')+','
        string += note.replace(',', '_')+','
        string += '{0:.2f}'.format(pi)+','
        string += modifier+','
        string += category+'\n'
        f.write(string)


def get_record_meta(record_list):
    """Get meta data of FASTA record.

    """
    acc_code = record_list[0]
    organism = record_list[1]
    EC_code = record_list[2].replace("__", " ")
    species = record_list[3].replace("__", " ")
    note = record_list[4]
    return acc_code, organism, EC_code, species, note


def prepare_out_csv(output_dir, param_dict):
    """Prepare headers for output CSV file. Overwrites existing files.

    """
    string = ''
    for i in param_dict['out_columns_pi']:
        if i == param_dict['out_columns_pi'][-1]:
            string += i
        else:
            string += i+','
    string += '\n'
    with open(output_dir+param_dict['out_CSV_pi'], 'w') as f:
        f.write(string)
