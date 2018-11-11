#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Program that calculates the pI of an input sequence string or sequences in an
input FASTA file.

Author: Andrew Tarzia

Date Created: 24 Apr 2018

TODO: get another method for calculating exposure without structure?
TODO: generalise code to multiple target residues
("targ = convert_to_one_letter_code_sing('target_res')" needs to act on a list
"""

# %%
# IMPORTS
import glob
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import IUPAC
import time
# my imports
from ercollect import parameter_file as params
from ercollect import pi_functions as pi_f

# %%
print("=====================================================================")
print('                            pI calculator                            ')
print("=====================================================================")
print("This program calculates the pI of all sequences in all FASTA files in")
print("the current directory using the BIOPYTHON SeqIO and ProteinAnalysis ")
print("modules.")
print("=====================================================================")


# get parameters
param_dict = params.get_parameters()

# get input FASTA file names
database_names = []
for i in glob.glob("*fasta"):
    if "_mod" not in i:
        database_names.append(i)

# prepare output CSV file
redo_pi = True
if redo_pi is True:
    string = ''
    for i in param_dict['out_columns_pi']:
        if i == param_dict['out_columns_pi'][-1]:
            string += i
        else:
            string += i+','
    string += '\n'
    with open(param_dict['out_CSV_pi'], 'w') as f:
        f.write(string)

# need to fix the FASTA output format so BIOPYTHON can read it
# for each FASTA file
pi_f.fix_fasta(database_names=database_names)


# %%
# run program
modifications = param_dict['modifications']
# for each FASTA file
for file in database_names:
    count_sequences_done = 0
    total_start_time = time.time()
    # read the file
    # but to avoid memory issues
    # we will calculate the pI on the fly
    # using the bio python module
    file_mod = file.replace(".fasta", "_mod.fasta")
    with open(file_mod, "r") as handle:
        print("--- Doing:", file_mod, "---")
        for record in SeqIO.parse(handle, "fasta", alphabet=IUPAC.protein):
            start_time = time.time()
            record_list = record.description.split("|")
            # get meta data
            acc_code = record_list[0]
            organism = record_list[1]
            EC_code = record_list[2].replace("__", " ")
            species = record_list[3].replace("__", " ")
            note = record_list[4]
            # get unmodified pI
            seq = record.seq
            seq_obj = ProteinAnalysis(''.join(seq))
            pi = seq_obj.isoelectric_point()
            count_sequences_done += 1
            modifier = '0'
            if pi < param_dict['cutoff_pi']:
                category = '0'
            else:
                category = '1'
            # output to CSV
            with open(param_dict['out_CSV_pi'], 'a') as f:
                string = file+','
                string += acc_code+','
                string += organism+','
                string += EC_code+','
                string += species+','
                string += note+','
                string += '{0:.2f}'.format(pi)+','
                string += modifier+','
                string += category+',\n'
                f.write(string)

            # if the category is 1 - i.e. pi > cutoff
            # then we test modification - else we do not
            if category == '1':
                modifier = '1'
                # get modified pI
                seq = record.seq
                # replace target amino acid residue
                # with replacement amino acid residue
                # one letter codes
                targ = pi_f.convert_to_one_letter_code_sing(modifications[modifier]['target_res'])
                replacement = pi_f.convert_to_one_letter_code_sing(modifications[modifier]['replace_res'])
                mod_seq = ''.join(seq).replace(targ, replacement)
                seq_obj = ProteinAnalysis(mod_seq)
                pi = seq_obj.isoelectric_point()
                count_sequences_done += 1
                if pi < param_dict['cutoff_pi']:
                    category = '0'
                else:
                    category = '1'
                # output to CSV
                with open(param_dict['out_CSV_pi'], 'a') as f:
                    string = file+','
                    string += acc_code+','
                    string += organism+','
                    string += EC_code+','
                    string += species+','
                    string += note+','
                    string += '{0:.2f}'.format(pi)+','
                    string += modifier+','
                    string += category+'\n'
                    f.write(string)

    print('--- finished %s sequences in %s seconds ---'
          % (count_sequences_done,
             '{0:.2f}'.format(time.time() - total_start_time)))
print('=====================================================================')
print('                              Complete!                              ')
print('=====================================================================')

# %%
# plot the data
pi_f.plot_pi_data(param_dict)
