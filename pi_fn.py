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
