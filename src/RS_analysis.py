#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse all RS properties.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import time
from os import getcwd
from os.path import join, exists
import sys

import reaction
import utilities


def main_analysis(prop_redo, pars):
    """
    Analyse all reaction systems.

    """

    print('------------------------------------------------------')
    print('collect component properties and analyse reaction systems:')
    print('    - diffusion of components')
    print('    - solubility (logP) of components')
    print('    - change in synthetic accessibility of components')
    print('------------------------------------------------------')

    print('settings:')
    print('    Diffusion threshold:', pars['size_thresh'], 'Angstrom')
    print('    Molecule database:', pars['molec_dir'])
    search_output_dir = getcwd()

    prop_output_file = join(search_output_dir, 'rs_properties.csv')

    if prop_redo or not exists(prop_output_file):
        with open(prop_output_file, 'w') as f:
            f.write(
                f'db_id,ec,max_mid_diam,minlogs,maxlogs,minlogp,'
                f'maxlogp,deltasa,rmaxsa,pmaxsa,nr,np,PC_class,'
                f'deltabct,rmaxbct,pmaxbct\n'
            )
        output_data = []
    else:
        output_data = []
        with open(prop_output_file, 'r') as f:
            for line in f.readlines():
                output_data.append(line)

    # in serial
    generator = reaction.yield_rxn_syst(search_output_dir, pars)
    for i, (count, rs) in enumerate(generator):
        if 'KEGG' not in rs.pkl:
            continue
        if rs.DB_ID in [i.rstrip().split(',')[0] for i in output_data]:
            continue
        if rs.skip_rxn:
            continue
        if rs.components is None:
            continue
        print('---------------------------------------------')
        print(f'checking RS {rs.pkl}, which is {i} of {count}')

        # Get all component properties from prop and diam files.
        rs.get_component_properties()
        if rs.max_min_mid_diam == 0:
            rs.fail_size()

        if rs.min_logP == 1E10 or rs.max_logP == -1E10:
            rs.fail_properties()

        if rs.min_logS == 1E10 or rs.max_logS == -1E10:
            rs.fail_properties()

        if rs.r_max_SA == 0 or rs.p_max_SA == 0:
            rs.fail_properties()

        # Save reaction system, but do not save data is failed.
        rs.save_reaction(filename=join(search_output_dir, rs.pkl))
        if rs.skip_rxn:
            continue

        no_react = len([
            i for i in rs.components if i.role == 'reactant'
        ])
        no_prods = len([
            i for i in rs.components if i.role == 'product'
        ])

        line_data = (
            f'{rs.DB_ID},{rs.EC},{rs.max_min_mid_diam},'
            f'{rs.min_logS},{rs.max_logS},{rs.min_logP},{rs.max_logP},'
            f'{rs.delta_SA},{rs.r_max_SA},{rs.p_max_SA},'
            f'{no_react},{no_prods},{rs.p_class},'
            f'{rs.delta_bCT},{rs.r_max_bCT},{rs.p_max_bCT}\n'
        )
        output_data.append(line_data)
        with open(prop_output_file, 'a') as f:
            f.write(line_data)


def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: RS_analysis.py prop_redo param_file
    prop_redo: t for rerun, else to read from prop_done.txt.
    param_file:
""")
        sys.exit()
    else:
        prop_redo = True if sys.argv[1] == 't' else False
        pars = utilities.read_params(sys.argv[2])

    temp_time = time.time()
    main_analysis(prop_redo=prop_redo, pars=pars)
    print(
        '--- time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )


if __name__ == "__main__":
    main()
