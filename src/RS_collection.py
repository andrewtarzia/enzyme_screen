#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to collect all RS.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import time
from os import getcwd
import sys

import utilities
import reaction
import KEGG_IO


def get_reaction_systems(
    EC,
    DB,
    params,
    output_dir,
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
        # rxn_fn = ATLAS_IO.get_rxn_systems

    rxn_fn(
        EC,
        output_dir,
        params=params,
        clean_system=clean_system,
        verbose=verbose
    )


def percent_skipped(output_dir, params):
    """
    Print the percent of all reaction systems that will NOT be skipped.

    """
    # what percentage of reaction systems have skip_rxn = False
    count_failed = 0
    total = 0
    for count, rs in reaction.yield_rxn_syst(output_dir, pars=params):
        print(f'{rs.DB_ID} ------ {rs.skip_rxn}')
        if rs.skip_rxn is False:
            count_failed += 1
        total += 1

    print('-----------------------------------')
    print(f'{count_failed} reaction systems of {total} pass!')
    print('=>', round(count_failed/total, 4)*100, 'percent')
    print('-----------------------------------')


def main_run(redo, pars):
    """
    Run reaction system collection.

    """

    print('--------------------------------------------------------')
    print('Screen new reactions')
    print('--------------------------------------------------------')
    temp_time = time.time()
    search_DBs = pars['DBs'].split('_')

    search_EC_file = pars['EC_file']

    print('--- settings:')
    print('    EC file:', search_EC_file)
    print('    DBs to search:', search_DBs)

    print('--- collect all reaction systems (ONLINE)...')
    search_ECs = utilities.get_ECs_from_file(EC_file=search_EC_file)
    search_output_dir = getcwd()
    for DB in search_DBs:
        for EC in search_ECs:
            print(f'------- searching for EC: {EC}')
            get_reaction_systems(
                EC=EC,
                DB=DB,
                params=pars,
                output_dir=search_output_dir,
                clean_system=redo,
                verbose=True
            )
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
        percent_skipped(search_output_dir, pars)

    print('----- All done! ------')


if __name__ == "__main__":
    main()
