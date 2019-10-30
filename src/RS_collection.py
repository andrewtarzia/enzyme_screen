#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to collect all RS.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import time
from multiprocessing import Pool
import pickle
import gzip
import glob
from os.path import isfile
from os import getcwd
import pandas as pd
import sys
from ercollect.molecule import (
    read_molecule_lookup_file,
    write_lookup_files
)


def get_reaction_systems(EC, DB, output_dir, molecule_dataset,
                         clean_system=False,
                         verbose=False):
    """Get reaction system from SABIO reaction ID (rID).

    Keywords:
        EC (str) - Enzyme commision number (X.X.X.X)
        DB (str) - name of Database
        output_dir (str) - directory where all data should be saved
        molecule_dataset (Pandas DataFrame) -
        look up for known molecules
        clean_system (bool) - wipe the data in reaction systems
        for fresh start
            default = False
        verbose (bool) - print update
            default = False

    """
    if DB == 'SABIO':
        print('SABIO DB cannot be used at this current time.')
        print(' - KEGG OR ATLAS only')
        sys.exit('exitting.')
        from ercollect.SABIO_IO import get_rxn_systems
    elif DB == 'KEGG':
        print('searching for EC:', EC)
        from ercollect.KEGG_IO import get_rxn_systems
    elif DB == 'BKMS':
        print('BKMS DB cannot be used at this current time.')
        print(' - KEGG OR ATLAS only')
        sys.exit('exitting.')
        from ercollect.BKMS_IO import get_rxn_systems
    elif DB == 'BRENDA':
        print('BRENDA DB cannot be used at this current time.')
        print(' - KEGG OR ATLAS only')
        sys.exit('exitting.')
        from ercollect.BRENDA_IO import get_rxn_systems
    elif DB == 'ATLAS':
        print('searching for EC:', EC)
        from ercollect.ATLAS_IO import get_rxn_systems
    get_rxn_systems(EC, output_dir, molecule_dataset=molecule_dataset,
                    clean_system=clean_system,
                    verbose=verbose)


def get_RS(filename, output_dir, verbose=False):
    """Read in reaction system from filename.

    """
    _rsf = filename.replace(output_dir+'sRS-', '').replace('.gpkl', '')
    EC_, DB, DB_ID = _rsf.split('-')
    EC = EC_.replace("_", ".").replace('XX', '-')
    rs = reaction(EC, DB, DB_ID)
    if isfile(output_dir+rs.pkl) is False:
        print('you have not collected all reaction systems.')
        print('Exitting.')
        import sys
        sys.exit()
    # load in rxn system
    if verbose:
        print('loading:', rs.pkl)
    rs = rs.load_object(output_dir+rs.pkl, verbose=False)
    return rs


def percent_skipped(output_dir):
    """Print the percent of all reaction systems that will NOT be skipped.

    """
    # what percentage of reaction systems have skip_rxn = False
    count = 0
    count_atlas = 0
    count_brenda = 0
    count_bkms = 0
    count_kegg = 0
    count_sabio = 0
    react_syst_files = glob.glob(output_dir+'sRS-*.gpkl')
    rsf_atlas = glob.glob(output_dir+'sRS-*ATLAS*.gpkl')
    rsf_brenda = glob.glob(output_dir+'sRS-*BRENDA*.gpkl')
    rsf_bkms = glob.glob(output_dir+'sRS-*BKMS*.gpkl')
    rsf_kegg = glob.glob(output_dir+'sRS-*KEGG*.gpkl')
    rsf_sabio = glob.glob(output_dir+'sRS-*SABIO*.gpkl')
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is False:
            count += 1
            if 'ATLAS' in rs.pkl:
                count_atlas += 1
            if 'BRENDA' in rs.pkl:
                count_brenda += 1
            if 'BKMS' in rs.pkl:
                count_bkms += 1
            if 'KEGG' in rs.pkl:
                count_kegg += 1
            if 'SABIO' in rs.pkl:
                count_sabio += 1

    print('-----------------------------------')
    print(count, 'reaction systems of', len(react_syst_files),
          'are NOT skipped.')
    print('=>', round(count/len(react_syst_files), 4)*100, 'percent')
    print('-----------------------------------')
    print(count_atlas, 'reaction systems of', len(rsf_atlas),
          'are NOT skipped in the ATLAS data set.')
    print('-----------------------------------')
    print(count_brenda, 'reaction systems of', len(rsf_brenda),
          'are NOT skipped in the BRENDA data set.')
    print('-----------------------------------')
    print(count_bkms, 'reaction systems of', len(rsf_bkms),
          'are NOT skipped in the BKMS data set.')
    print('-----------------------------------')
    print(count_kegg, 'reaction systems of', len(rsf_kegg),
          'are NOT skipped in the KEGG data set.')
    print('-----------------------------------')
    print(count_sabio, 'reaction systems of', len(rsf_sabio),
          'are NOT skipped in the SABIO data set.')
    print('-----------------------------------')


def get_ECs_from_file(EC_file):
    """Read in ECs to search from a file.

    """
    # get search EC numbers from file:
    # set EC numbers of interest
    # get from a data file - manually made from
    # https://enzyme.expasy.org/enzyme-byclass.html
    EC_DF = pd.read_table(
        EC_file,
        delimiter='__',
        names=['EC_no', 'description'],
        engine='python'
    )
    search_ECs = list(EC_DF['EC_no'])

    # remove all spaces within EC numbers
    search_ECs = [i.replace(' ', '') for i in search_ECs]

    # add check for '1' from '1.-.-.-'
    new_search_ECs = []
    for EC in search_ECs:
        if '-' in EC:
            new_search_ECs.append(EC.replace('.-', ''))
            new_search_ECs.append(EC)
        else:
            new_search_ECs.append(EC)

    print(len(search_ECs), 'EC numbers to test')
    print('first EC:', search_ECs[0], '---- last EC:', search_ECs[-1])
    print('collect all reaction systems (ONLINE)...')
    return new_search_ECs


def main_run(redo):
    """Run reaction system collection.

    """
    if redo == 'T':
        redo = True
    else:
        redo = False
    print('--------------------------------------------------------')
    print('Screen new reactions')
    print('--------------------------------------------------------')
    temp_time = time.time()
    DB_switch = input(
        'biomin (1) or new (2) or SABIO (3) or KEGG (4) or ATLAS (5)?'
    )
    if DB_switch == '1':
        search_DBs = ['BRENDA', 'SABIO', 'KEGG', 'BKMS', ]
    elif DB_switch == '2':
        search_DBs = ['SABIO', 'ATLAS', 'KEGG', 'BRENDA', 'BKMS']
    elif DB_switch == '3':
        search_DBs = ['SABIO']
    elif DB_switch == '4':
        search_DBs = ['KEGG']
    elif DB_switch == '5':
        search_DBs = ['ATLAS', ]
    else:
        print('answer correctly...')
        sys.exit()
    NP = 1  # number of processes
    search_EC_file = 'desired_EC.txt'
    lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
    translator = '/home/atarzia/psp/molecule_DBs/KEGG/translator.txt'
    molecule_DB_directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
    # write molecule look up files based on molecule DB
    write_lookup_files(lookup_file, translator, molecule_DB_directory)
    molecule_dataset = read_molecule_lookup_file(
        lookup_file=lookup_file
    )
    print('settings:')
    print('    EC file:', search_EC_file)
    print('    Number of processes:', NP)
    print('    DBs to search:', search_DBs)
    print('    Molecule DB lookup file:', lookup_file)
    # inp = input('happy with these? (T/F)')
    # if inp == 'F':
    #     sys.exit('change them in the source code')
    # elif inp != 'T':
    #     sys.exit('I dont understand, T or F?')
    print('collect all reaction systems (ONLINE)...')
    search_ECs = get_ECs_from_file(EC_file=search_EC_file)
    search_output_dir = getcwd()+'/'
    for DB in search_DBs:
        # iterate over EC numbers of interest
        if NP > 1:
            # Create a multiprocessing Pool
            with Pool(NP) as pool:
                # process data_inputs iterable with pool
                # func(EC, DB, search_output_dir, mol dataset,
                # search_redo,
                #      verbose)
                args = [
                    (
                        EC,
                        DB,
                        search_output_dir,
                        molecule_dataset,
                        redo,
                        True
                    )
                    for EC in search_ECs
                ]
                pool.starmap(get_reaction_systems, args)
        # in serial
        else:
            for EC in search_ECs:
                get_reaction_systems(EC=EC, DB=DB,
                                     output_dir=search_output_dir,
                                     molecule_dataset=molecule_dataset,
                                     clean_system=redo, verbose=True)
    percent_skipped(search_output_dir)
    print('---- time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')


def main():
    if (not len(sys.argv) == 4):
        print("""
Usage: RS_collection.py run redo skipped
    run: T to run search for new rxn systems into current dir.
    redo: T to overwrite all rxn systems.
    skipped: T to see the number of skipped rxns in cwd.
""")
        sys.exit()
    else:
        run = sys.argv[1]
        redo = sys.argv[2]
        skipped = sys.argv[4]

    if run == 'T':
        main_run(redo)
    if skipped == 'T':
        search_output_dir = getcwd()+'/'
        percent_skipped(search_output_dir)

    print('All done!')


if __name__ == "__main__":
    main()
