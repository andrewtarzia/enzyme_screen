#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to wipe all RS properties.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import glob
from os import getcwd
import sys


def wipe_reaction_properties(rs, output_dir):
    """Set attributes of rxn system to None.

    """
    print('wiping: r_max_sa, p_max_sa, delta_sa.')
    # rs.skip_rxn = False
    # rs.all_fit = None  # do all the components fit?
    # rs.max_comp_size = None
    rs.r_max_sa = None
    rs.p_max_sa = None
    rs.delta_sa = None
    rs.save_object(output_dir+rs.pkl)


def main_wipe():
    """Wipe reaction system properties to rerun analysis

    """
    print('------------------------------------------------')
    print('Wipe reaction properties')
    print('------------------------------------------------')
    inp = input('are you sure? (T/F)')
    if inp == 'F':
        sys.exit('')
    elif inp != 'T':
        sys.exit('I dont understand, T or F?')
    search_output_dir = getcwd()+'/'
    react_syst_files = glob.glob(search_output_dir+'sRS-*.gpkl')
    count = 0
    for rs in yield_rxn_syst(search_output_dir):
        print('wiping', count, 'of', len(react_syst_files))
        count += 1
        wipe_reaction_properties(rs, search_output_dir)


def main():
    main_wipe()

    print('All done!')


if __name__ == "__main__":
    main()
