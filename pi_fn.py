#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for pI calculation program.

Author: Andrew Tarzia

Date Created: 24 Apr 2018

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import time
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import IUPAC
from plotting import EC_descriptions


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


def plot_pI_dist(pi_data, filename, output_dir, cutoff_pi):
    """Plot and save histogram of all pI data calculated in this run.

    Arguments:

    """
    fig, ax = plt.subplots(figsize=(8, 5))

    modifications = define_seq_modifications()

    # unmodifed
    mod_dict = modifications['0']
    data = pi_data[pi_data['modification'] == 0]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               histtype='stepfilled',
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # modification 1 - succinylation
    mod_dict = modifications['1']
    data = pi_data[pi_data['modification'] == 1]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               histtype='stepfilled',
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.set_xlabel('calculated pI', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 14)
    # plot pI cut-off
    ax.axvline(x=cutoff_pi, c='k', lw='2', linestyle='--')
    # legend
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(output_dir+filename,
                dpi=720, bbox_inches='tight')


def plot_EC_pI_dist(EC_pi_data, filename, title, cutoff_pi):
    """Plot and save histogram of all pI data calculated in this run separated
    by EC class.

    Arguments:

    """
    fig, ax = plt.subplots()

    modifications = define_seq_modifications()

    # unmodifed
    mod_dict = modifications['0']
    data = EC_pi_data[EC_pi_data['modification'] == 0]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               histtype='stepfilled',
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # modification 1 - succinylation
    mod_dict = modifications['1']
    data = EC_pi_data[EC_pi_data['modification'] == 1]
    n, bins, patches = ax.hist(data['pi'],
                               facecolor=mod_dict['colour'],
                               alpha=0.5,
                               histtype='stepfilled',
                               bins=np.arange(0, 14 + 0.2, 0.5),
                               label=mod_dict['name'])

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('calculated pI', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(0, 14)
    # plot pI cut-off
    ax.axvline(x=cutoff_pi, c='k', lw='2', linestyle='--')
    # legend
    ax.legend(fontsize=16)
    # title
    ax.set_title(title, fontsize=16)

    fig.tight_layout()
    fig.savefig(filename,
                dpi=720, bbox_inches='tight')


def define_seq_modifications():
    """Define possible sequence modifications.

    # modification types + colours
    # just implement succinylation for now
    # why would you want to do acetylation if you can succinylate??
    # currently succinylation == LYS swapped with GLU

    """
    modifications = {'0': {
            'colour': 'k',
            'name': 'unmodified',
        },
                     '1': {
            'colour': 'firebrick',
            'name': 'succinylated',
            'target_res': 'LYS',
            'replace_res': 'GLU',
        }
        }

    return modifications


def calculate_pI_from_file(file, output_dir, cutoff_pi, out_CSV_pi):
    """Calculate the pI of all sequences in FASTA file.

    """
    modifications = define_seq_modifications()
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
            if pi < cutoff_pi:
                category = '0'
            else:
                category = '1'
            # output to CSV
            output_pI_row(output_dir, out_CSV_pi, file,
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
                if pi < cutoff_pi:
                    category = '0'
                else:
                    category = '1'
                # output to CSV
                output_pI_row(output_dir, out_CSV_pi, file,
                              acc_code, organism, EC_code,
                              species, note,
                              pi, modifier, category)
            # break
    print('--- finished %s sequences in %s seconds ---'
          % (count_sequences_done, '{0:.2f}'.format(time.time() - total_start_time)))


def calculate_rxn_syst_pI(sequence, rxn_syst, cutoff_pi):
    """Calculate the pI of a sequence associated with a reaction system.

    """
    modifications = define_seq_modifications()
    seq_obj = ProteinAnalysis(sequence)
    pi = seq_obj.isoelectric_point()
    modifier = '0'
    if pi < cutoff_pi:
        category = '0'
    else:
        category = '1'

    if category == '0':
        rxn_syst.seed_MOF = True
        rxn_syst.pI = pi

    # if the category is 1 - i.e. pi > cutoff
    # then we test modification
    elif category == '1':
        # report unmodified pI if modification isn't successful
        rxn_syst.pI = pi
        modifier = '1'
        # get modified pI
        seq = sequence
        # replace target amino acid residue
        # with replacement amino acid residue
        # one letter codes
        targ = convert_to_one_letter_code_sing(modifications[modifier]['target_res'])
        replacement = convert_to_one_letter_code_sing(modifications[modifier]['replace_res'])
        mod_seq = ''.join(seq).replace(targ, replacement)
        seq_obj = ProteinAnalysis(mod_seq)
        pi = seq_obj.isoelectric_point()
        if pi < cutoff_pi:
            category = '0'
        else:
            category = '1'

        if category == '0':
            rxn_syst.seed_MOF = True
            rxn_syst.req_mod = modifier
            rxn_syst.pI = pi
        else:
            rxn_syst.seed_MOF = False

    return rxn_syst


def output_pI_row(output_dir, out_CSV_pi, file,
                  acc_code, organism, EC_code, species, note,
                  pi, modifier, category):
    """Output results of pI calculation to CSV.

    """
    with open(output_dir+out_CSV_pi, 'a') as f:
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


def prepare_out_csv(output_dir, filename):
    """Prepare headers for output CSV file. Overwrites existing files.

    """
    out_columns_pi = ['fasta_file', 'acc.code',
                      'organism', 'EC.code', 'species',
                      'note', 'pi', 'modification', 'category']
    string = ''
    for i in out_columns_pi:
        if i == out_columns_pi[-1]:
            string += i
        else:
            string += i+','
    string += '\n'
    with open(output_dir+filename, 'w') as f:
        f.write(string)


def prepare_pI_calc(database_directory, redo_pi, output_dir, csv):
    """Prepare for pI screening.

    """
    # get input FASTA file names
    database_names = []
    for i in glob.glob(database_directory+"*fasta"):
        if "_mod" not in i:
            database_names.append(i)
    database_names = sorted(database_names)
    print('databases:')
    for i in database_names:
        print('--', i.replace(database_directory, ''))

    # prepare output CSV file

    if redo_pi is True:
        prepare_out_csv(output_dir, csv)
        # fix formatting of FASTA files to match BIOPYTHON readable
        fix_fasta(database_names)

    return database_names


def screen_pIs(database_names, redo_pI, redo_pI_plots, pI_csv, pI_output_dir,
               cutoff_pi, descriptors):
    """Screen the pI of all sequences with chosen EC numbers.

    """
    if descriptors is None:
        descriptors = {}
    for EC_file in database_names:
        EC = EC_file.replace(pI_output_dir, '')
        EC = EC.replace('__BRENDA_sequences.fasta', '').replace('_', '.')
        top_EC = EC.split('.')[0]
        # read the file but to avoid memory issues # we will calculate the pI
        # on the fly using the bio python module
        print('doing:', EC_file)
        file_mod = EC_file.replace(".fasta", "_mod.fasta")
        if redo_pI is True:
            calculate_pI_from_file(file_mod, pI_output_dir,
                                   cutoff_pi, pI_csv)
        if redo_pI_plots is True:
            print('plot distribution of pIs')
            pi_data = pd.read_csv(pI_output_dir+pI_csv, index_col=False)
            EC_pi_data = pi_data[pi_data['fasta_file'] == file_mod]
            plot_EC_pI_dist(EC_pi_data,
                            filename=file_mod.replace('.fasta', '.pdf'),
                            title=EC_descriptions()[top_EC][0],
                            cutoff_pi=cutoff_pi)
        print('done')
    if redo_pI_plots is True:
        print('plot full distribution of pIs')
        pi_data = pd.read_csv(pI_output_dir+pI_csv, index_col=False)
        plot_pI_dist(pi_data,
                     filename='full_pI_dist.pdf',
                     output_dir=pI_output_dir,
                     cutoff_pi=cutoff_pi)
