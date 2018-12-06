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
        with open(FASTA_file, 'r') as f:
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
    output.to_csv(output_file, index=False, sep='@')
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
        if input('run calculations? (t/f)') == 't':
            run = True
        else:
            run = False
    else:
        # overwrite output file
        output = write_seq_output(output_file)
    print('-----------------------------------------------------------')
    print('doing calculations...')
    # need to fix the FASTA output format so BIOPYTHON can read it
    file_mod = fix_fasta(FASTA_file=fasta_file)
    total_start_time = time.time()
    total_seq = 0
    if run is True:
        # iterate through sequences in FASTA file
        done = list(output.acc_code)
        with open(file_mod, "r") as handle:
            for record in SeqIO.parse(handle, "fasta", alphabet=IUPAC.protein):
                total_seq += 1
                record_list = record.description.split("|")
                # collect information on sequence and sequence
                # should be a unique descriptor
                acc_code = record_list[0].lstrip().rstrip()
                organism = record_list[1].lstrip().rstrip()
                EC_code = record_list[2].replace("__", " ").lstrip().rstrip()
                species = record_list[3].replace("__", " ").lstrip().rstrip()
                note = record_list[4]
                if acc_code in done:
                    continue
                seq = record.seq
                sequence_string = str(seq)
                seq_obj = ProteinAnalysis(''.join(seq))
                # do calculations
                pI = seq_obj.isoelectric_point()
                GRAVY = seq_obj.gravy()
                I_index = seq_obj.instability_index()
                A_index = calculate_seq_aliphatic_index(sequence_string)
                TM_index = calculate_TM_index(seq_string=sequence_string)
                print(seq)
                print(pI, GRAVY, I_index, A_index, TM_index)
                print(record_list)
                input('ok?')
                done.append(acc_code)
                ROW = pd.DataFrame({'acc_code': acc_code, 'organism': organism,
                                    'EC_code': EC_code,  'species': species,
                                    'note': note, 'pI': pI, 'GRAVY': GRAVY,
                                    'I_index': I_index, 'A_index': A_index,
                                    'TM_index': TM_index}, index=[0])
                print(ROW)
                input('pl?')
                # save to output file
                output = update_seq_output(output_file, output, ROW)
    print('--- finished %s sequences in %s seconds ---'
          % (total_seq, '{0:.2f}'.format(time.time() - total_start_time)))


def fasta_plotting(output_file, plot_suffix):
    """Plot data for bulk sequence analysis. output_file defines file names,
    titles.

    """
    print('ti')


def main_analysis(plot_suffix, fasta_file, output_file):
    """Analyse all sequences in FASTA file from BRENDA.

    """
    print('--------------------------------------------------------------')
    print('Analyse properties of all sequences in FASTA:', fasta_file)
    print('--------------------------------------------------------------')
    # percent_w_sequence(output_dir=search_output_dir)
    temp_time = time.time()
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

    output_file = FASTA.replace('_sequences.fasta', '_output.csv')
    main_analysis(plot_suffix=plot_suffix, fasta_file=FASTA,
                  output_file=output_file)
    sys.exit()
