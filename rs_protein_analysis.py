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
from rxn_syst import yield_rxn_syst, load_molecule
import plotting


def percent_w_sequence(output_dir):
    """Print the percent of all reaction systems with a protein sequence.

    """
    # what percentage of reaction systems have skip_rxn = False
    count = 0
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
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
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
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


def wipe_reaction_properties(rs, output_dir):
    """Set attributes of rxn system to None.

    """
    rs.save_object(output_dir+rs.pkl)


def main_wipe():
    """Wipe reaction system properties to rerun analysis

    """
    print('NEED TO DETERMINE WHAT IS TO BE WIPED')
    sys.exit()
    print('--------------------------------------------------------------')
    print('Wipe reaction properties')
    print('--------------------------------------------------------------')
    search_output_dir = os.getcwd()+'/'
    react_syst_files = glob.glob(search_output_dir+'sRS-*.pkl')
    count = 0
    for rs in yield_rxn_syst(search_output_dir):
        print('wiping', count, 'of', len(react_syst_files))
        count += 1
        wipe_reaction_properties(rs, search_output_dir)


def main_analysis():
    """Analyse all reaction systems.

    """
    print('--------------------------------------------------------------')
    print('Collect reaction properties')
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
    percent_w_sequence(output_dir=search_output_dir)
    temp_time = time.time()
    print('get disitribution of pIs of all known protein sequences...')
    plotting.rs_pI_distribution(output_dir=search_output_dir,
                                cutoff_pI=pI_thresh,
                                generator=yield_rxn_syst(search_output_dir))
    print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
    temp_time = time.time()
    print('get subset of reactions with known protein sequences...')
    check_all_seedMOF(search_output_dir, pI_thresh)
    print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
    temp_time = time.time()

    print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')


if __name__ == "__main__":
    import sys
    import time
    from multiprocessing import Pool
    from molecule import molecule
    from rxn_syst import reaction

    if (not len(sys.argv) == 3):
        print('Usage: rs_protein_analysis.py properties wipe\n')
        print('   properties: T to get properties of reaction systems in cwd.')
        print('   wipe: T to wipe properties of reaction systems in cwd.')
        sys.exit()
    else:
        properties = sys.argv[1]
        wipe = sys.argv[2]

    if wipe == 'T':
        main_wipe()
    if properties == 'T':
        main_analysis()
