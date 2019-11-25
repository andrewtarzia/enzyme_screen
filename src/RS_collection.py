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
import glob
from os.path import join
from os import getcwd
import pandas as pd
import sys

import IO
import utilities
import reaction
import KEGG_IO


def get_reaction_systems(
    EC,
    DB,
    params,
    output_dir,
    molecule_dataset,
    clean_system=False,
    verbose=False
):
    """
    Get reaction system from SABIO reaction ID (rID).

    Keywords:
        EC (str) :
            Enzyme commision number (X.X.X.X)
        DB (str) :
            name of Database
        output_dir (str) :
            directory where all data should be saved
        molecule_dataset (pd.DataFrame) :
            dataset of existing molecules to use
        clean_system (bool) :
            wipe the data in reaction systems for fresh start
            (default = False)
        verbose (bool) :
            print update to terminal
            (default = False)

    """

    if DB == 'SABIO':
        raise NotImplementedError(
            'SABIO DB cannot be used at this current time.'
        )
        # rxn_fn = SABIO_IO.get_rxn_systems
    elif DB == 'KEGG':
        print('searching for EC:', EC)
        rxn_fn = KEGG_IO.get_rxn_systems
    elif DB == 'BKMS':
        print('BKMS DB cannot be used at this current time.')
        print(' - KEGG OR ATLAS only')
        raise NotImplementedError(
            'BKMS DB cannot be used at this current time.'
        )
        # rxn_fn = BKMS_IO.get_rxn_systems
    elif DB == 'BRENDA':
        raise NotImplementedError(
            'BRENDA DB cannot be used at this current time.'
        )
        # rxn_fn = BRENDA_IO.get_rxn_systems
    elif DB == 'ATLAS':
        raise NotImplementedError(
            'ATLAS DB cannot be used at this current time.'
        )
        print('searching for EC:', EC)
        # rxn_fn = ATLAS_IO.get_rxn_systems

    rxn_fn(
        EC,
        output_dir,
        params=params,
        molecule_dataset=molecule_dataset,
        clean_system=clean_system,
        verbose=verbose
    )


def percent_skipped(output_dir):
    """
    Print the percent of all reaction systems that will NOT be skipped.

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
    for rs in reaction.yield_rxn_syst(output_dir):
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
    print(
        f'{count} reaction systems of '
        f'{len(react_syst_files)} are NOT skipped.'
    )
    print('=>', round(count/len(react_syst_files), 4)*100, 'percent')
    print('-----------------------------------')
    print(
        f'{count_atlas} reaction systems of '
        f'{len(rsf_atlas)} are NOT skipped in the ATLAS data set.'
    )
    print('-----------------------------------')
    print(
        f'{count_brenda} reaction systems of '
        f'{len(rsf_brenda)} are NOT skipped in the BRENDA data set.'
    )
    print('-----------------------------------')
    print(
        f'{count_bkms} reaction systems of '
        f'{len(rsf_bkms)} are NOT skipped in the BKMS data set.'
    )
    print('-----------------------------------')
    print(
        f'{count_kegg} reaction systems of '
        f'{len(rsf_kegg)} are NOT skipped in the KEGG data set.'
    )
    print('-----------------------------------')
    print(
        f'{count_sabio} reaction systems of '
        f'{len(rsf_sabio)} are NOT skipped in the SABIO data set.'
    )
    print('-----------------------------------')


def get_ECs_from_file(EC_file):
    """
    Read in ECs to search from a file.

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


def main_run(redo, pars):
    """
    Run reaction system collection.

    """

    print('--------------------------------------------------------')
    print('Screen new reactions')
    print('--------------------------------------------------------')
    temp_time = time.time()
    search_DBs = pars['DBs'].split('_')

    NP = 1  # number of processes
    search_EC_file = pars['EC_file']
    translator = pars['translator_kegg']
    molecule_DB_directory = pars['molec_dir']
    lookup_file = join(molecule_DB_directory, 'lookup.txt')
    # write molecule look up files based on molecule DB
    IO.write_lookup_files(
        lookup_file,
        translator,
        molecule_DB_directory
    )
    molecule_dataset = IO.read_molecule_lookup_file(
        lookup_file=lookup_file
    )
    print('--- settings:')
    print('    EC file:', search_EC_file)
    print('    Number of processes:', NP)
    print('    DBs to search:', search_DBs)
    print('    Molecule DB lookup file:', lookup_file)

    print('--- collect all reaction systems (ONLINE)...')
    search_ECs = get_ECs_from_file(EC_file=search_EC_file)
    search_output_dir = getcwd()
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
                        pars,
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
                get_reaction_systems(
                    EC=EC,
                    DB=DB,
                    params=pars,
                    output_dir=search_output_dir,
                    molecule_dataset=molecule_dataset,
                    clean_system=redo,
                    verbose=True
                )
    percent_skipped(search_output_dir)
    print(
        '---- time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )


def main():
    if (not len(sys.argv) == 5):
        print("""
Usage: RS_collection.py run redo skipped
    run: t to run search for new rxn systems into current dir.
    redo: t to overwrite all rxn systems.
    param_file:
    skipped: t to see the number of skipped rxns in cwd.
""")
        sys.exit()
    else:
        run = True if sys.argv[1] == 't' else False
        redo = True if sys.argv[2] == 't' else False
        pars = utilities.read_params(sys.argv[3])
        skipped = True if sys.argv[4] == 't' else False

    if run:
        main_run(redo, pars)
    if skipped:
        search_output_dir = getcwd()+'/'
        percent_skipped(search_output_dir)

    print('----- All done! ------')


if __name__ == "__main__":
    main()
