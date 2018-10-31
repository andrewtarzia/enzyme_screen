#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the reaction system class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import glob
import os
import Uniprot_IO
import pi_fn
from rxn_syst import yield_rxn_syst, yield_rxn_syst_filelist
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
from tm_predictor import calculate_TM_index


def percent_w_sequence(output_dir):
    """Print the percent of all reaction systems with a protein sequence.

    """
    # what percentage of reaction systems have skip_rxn = False
    count = 0
    react_syst_files = sorted(glob.glob(output_dir+'sRS-*.gpkl'))
    for rs in yield_rxn_syst(output_dir):
        if rs.UniprotID != '':
            if rs.UniprotID is not None:
                count += 1

    print('-----------------------------------')
    print(count, 'reaction systems of', len(react_syst_files),
          'had a sequence.')
    print('=>', round(count/len(react_syst_files), 4)*100, 'percent')
    print('-----------------------------------')


def get_pI_of_UNIRPOT_IDs(UniProtIDs, rs, pI_thresh, output_dir):
    """Get the pI of all UniProtIDs summed together. Checks for precalculated
    pI.

    """
    print('move to UNIPROT')
    # iterate over all UniProtIDs
    # collate all sequences into one pI
    total_sequence = ''
    for i in UniProtIDs:
        sequence = Uniprot_IO.get_sequence(i)
        total_sequence += sequence
    if len(total_sequence) > 0:
        rs = pi_fn.calculate_rxn_syst_pI(total_sequence, rs,
                                         cutoff_pi=pI_thresh)
        print('seed MOF?', rs.seed_MOF)
    else:
        rs.pI = None
        rs.seed_MOF = 'unclear'
    rs.save_object(output_dir+rs.pkl)


def check_all_seedMOF(output_dir, pI_thresh):
    """Check all reaction systems an associated protein sequence and whether it
    seeds MOF growth.

    Keywords:
        output_dir (str) - directory to output molecule files
        pI_thresh (float) - thrshold pI for MOF growth

    """
    # iterate over reaction system files
    react_syst_files = sorted(glob.glob(output_dir+'sRS-*.gpkl'))
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            continue
        print('change this based on other IDs used')
        # this is only possible for reaction systems with an ID for a sequence
        # SABIO uses Uniprot
        if rs.UniprotID is None or rs.UniprotID == '':
            continue
        # KEGG AND ATLAS?
        print("#####################")
        # pI already checked?
        if rs.seed_MOF is not None:
            continue
        print('checking rxn', count, 'of', len(react_syst_files))
        # DBs for which protein sequences were possible
        if rs.DB == 'SABIO':
            # split UniprotID for the cases where multiple subunits exist
            IDs = rs.UniprotID.split(" ")
            print('Uniprot IDs:', IDs)
            if len(IDs) > 0:
                get_pI_of_UNIRPOT_IDs(UniProtIDs=IDs, rs=rs,
                                      pI_thresh=pI_thresh,
                                      output_dir=output_dir)
            else:
                rs.save_object(output_dir+rs.pkl)
                print('seed MOF?', rs.seed_MOF)
        if rs.DB == 'KEGG':
            continue
        if rs.DB == 'ATLAS':
            continue


def IDs2sequence(UniProtIDs):
    """Obtain Uniprot sequences from IDs and add together.

    """
    print('move to UNIPROT')
    # iterate over all UniProtIDs
    # collate all sequences into one pI
    total_sequence = ''
    for i in UniProtIDs:
        sequence = Uniprot_IO.get_sequence(i)
        total_sequence += sequence
    if len(total_sequence) > 0:
        return total_sequence
    else:
        return None


def calculate_seq_aliphatic_index(seq_string):
    """Calculate aliphatic index of an amino acid sequence (str).

    Defined in this paper: https://www.jstage.jst.go.jp/article/biochemistry1922/88/6/88_6_1895/_article
    Code inspired by: https://github.com/ddofer/ProFET/blob/master/ProFET/feat_extract/ProtFeat.py

    """
    # set parameters
    a = 2.9
    b = 3.9
    length = float(len(seq_string))
    AA_Counts = Counter(seq_string)
    alanine_per = (AA_Counts['A'] / length)
    valine_per = (AA_Counts['V'] / length)
    isoleucine_per = (AA_Counts['I'] / length)
    leucine_per = (AA_Counts['L'] / length)
    # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
    aliphatic_index = (100 * (alanine_per + a * valine_per + b * (isoleucine_per + leucine_per )))
    return aliphatic_index


def get_RS_sequence_properties(output_dir, filelist):
    """Get sequence properties for all reaction systems with an associated
    protein sequence.

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
    # iterate over reaction system files
    react_syst_files = sorted(glob.glob(output_dir+'sRS-*.gpkl'))
    count = 0
    for rs in yield_rxn_syst_filelist(output_dir, filelist, verbose=True):
        count += 1
        if rs.skip_rxn is True:
            continue
        print('change this based on other IDs used')
        # this is only possible for reaction systems with an ID for a sequence
        # SABIO uses Uniprot
        if rs.UniprotID is None or rs.UniprotID == '':
            continue
        # KEGG AND ATLAS?
        print("#####################")
        # sequence properties checked?
        try:
            if rs.sequence is not None:
                continue
        except AttributeError:
            rs.sequence = None
            rs.save_object(output_dir+rs.pkl)
            pass
        print('checking rxn', count, 'of', len(react_syst_files))
        # DBs for which protein sequences were possible
        if rs.DB == 'SABIO':
            # split UniprotID for the cases where multiple subunits exist
            IDs = rs.UniprotID.split(" ")
            print('Uniprot IDs:', IDs)
            if len(IDs) > 0:
                if rs.sequence is None:
                    sequence_string = IDs2sequence(UniProtIDs=IDs)
                    if sequence_string is None:
                        rs.sequence = None
                        continue
                    # set RS sequence
                    rs.sequence = sequence_string
                else:
                    sequence_string = rs.sequence
                # convert to BioPYTHON object
                seq_obj = ProteinAnalysis(sequence_string)
                # get sequence properties
                # this does not include modified pI
                try:
                    print('calculating protein properties')
                    rs.pI = seq_obj.isoelectric_point()
                    rs.GRAVY = seq_obj.gravy()
                    rs.I_index = seq_obj.instability_index()
                    rs.A_index = calculate_seq_aliphatic_index(sequence_string)
                    rs.TM_index = calculate_TM_index(
                        seq_string=sequence_string)
                except KeyError:
                    print('sequence has non-natural amino acid.')
                    rs.pI = None
                    rs.GRAVY = None
                    rs.I_index = None
                    rs.A_index = None
                    rs.TM_index = None
                rs.save_object(output_dir+rs.pkl)
            else:
                rs.save_object(output_dir+rs.pkl)
        if rs.DB == 'KEGG':
            rs.save_object(output_dir+rs.pkl)
            continue
        if rs.DB == 'ATLAS':
            rs.save_object(output_dir+rs.pkl)
            continue


def wipe_reaction_properties(rs, output_dir):
    """Set attributes of rxn system to None.

    """
    print('wiping: sequence, pI, GRVAY, I_index, A_index.')
    rs.sequence = None
    rs.pI = None
    rs.GRAVY = None
    rs.I_index = None
    rs.A_index = None
    rs.save_object(output_dir+rs.pkl)


def main_wipe():
    """Wipe reaction system properties to rerun analysis

    """
    print('--------------------------------------------------------------')
    print('Wipe reaction properties')
    print('--------------------------------------------------------------')
    search_output_dir = os.getcwd()+'/'
    react_syst_files = sorted(glob.glob(search_output_dir+'sRS-*.gpkl'))
    count = 0
    for rs in yield_rxn_syst(search_output_dir):
        print('wiping', count, 'of', len(react_syst_files))
        count += 1
        wipe_reaction_properties(rs, search_output_dir)


def main_analysis(plot_suffix, file):
    """Analyse all reaction systems.

    """
    print('--------------------------------------------------------------')
    print('Collect reaction properties for reactions with known sequences')
    print('--------------------------------------------------------------')
    NP = 1  # number of processes
    pI_thresh = 6
    print('settings:')
    print('    Number of processes:', NP)
    print('    pI threshold:', pI_thresh)
    inp = input('happy with these? (T/F)')
    if inp == 'F':
        sys.exit('change them in the source code')
    elif inp != 'T':
        sys.exit('I dont understand, T or F?')
    search_output_dir = os.getcwd()+'/'
    # percent_w_sequence(output_dir=search_output_dir)
    temp_time = time.time()
    print('get sequence properties for reaction systems...')
    get_RS_sequence_properties(output_dir=search_output_dir, filelist=file)
    print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
    # print('get pI of sequences associated with reaction systems...')
    # check_all_seedMOF(search_output_dir, pI_thresh)
    # print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
    # temp_time = time.time()
    # print('get hydropathy of sequences associated with reaction systems...')
    # get_RS_sequence_hydropathy(search_output_dir)
    # print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')


if __name__ == "__main__":
    import sys
    import time
    from rxn_syst import reaction

    if (not len(sys.argv) == 5):
        print('Usage: rs_protein_analysis.py properties wipe plot_suffix RS_file\n')
        print('   properties: T to get properties of reaction systems in cwd.')
        print('   wipe: T to wipe properties of reaction systems in cwd.')
        print('   plot_suffix: string to put at the end of plot file names.')
        print('   RS_file: name of file containing list of RS.')
        sys.exit()
    else:
        properties = sys.argv[1]
        wipe = sys.argv[2]
        plot_suffix = sys.argv[3]
        RS_file = sys.argv[4]

    if wipe == 'T':
        main_wipe()
    if properties == 'T':
        main_analysis(plot_suffix=plot_suffix, file=RS_file)

    sys.exit()

    ### TESTING ATLAS AND KEGG
    # search_output_dir = '/home/atarzia/psp/screening_results/new_reactions/'
    # for rs in yield_rxn_syst(search_output_dir):
    #     if rs.skip_rxn:
    #         continue
    #     if rs.DB == 'KEGG':
    #         break
    #
    # rs.print_rxn_system()
    #
    # import requests
    # URL = 'http://rest.kegg.jp/get/reaction:R00963'
    # request = requests.post(URL)
    # request.raise_for_status()
    #
    # request.text
    # ortho_codes = []
    # if 'ORTHOLOGY' in request.text:
    #     first_part = request.text.split('ORTHOLOGY   ')[1].split('\n')
    #     for fir in first_part:
    #         print(fir.lstrip().rstrip())
    #         line = fir.lstrip().rstrip()
    #         if len(line) > 0 and line[0] == 'K':  # is orthology code
    #             print(line)
    #             ortho_code = line.split(' ')[0]
    #             print(ortho_code)
    #             ortho_codes.append(ortho_code)
    #
    # ortho_codes
    #
    # # get all genes associated with orthologies
    # for orth in ortho_codes:
    #     print(orth)
    #     URL = 'http://rest.kegg.jp/get/orthology:'+orth
    #     request = requests.post(URL)
    #     request.raise_for_status()
    #
    #     request.text
    #     break
