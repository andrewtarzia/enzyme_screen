#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to populate the properties of all molecules in database.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os import getcwd
from ercollect import rxn_syst
from plotting import (
    mol_SA_vs_compl,
    mol_SA_vs_NRB,
    mol_SA_vs_NHA,
    mol_logP_vs_logS,
    mol_all_dist,
    mol_logP_vs_XlogP
)
import sys


def main():
    if (not len(sys.argv) == 5):
        print("""
Usage: molecule.py get_mol pop_mol mol_file
    get_mol: T for overwrite and collection of molecules from RS in
        current dir
        (F for skip) ---- this function is useful if you update the
        base attributes of the molecule class.
    scratch: T for restart collection of molecules from RS.
    pop_mol: T to run population of molecule properties
        (does not overwrite).
    mol_file: file name of list of molecules, F if not specified.
""")
        sys.exit()
    else:
        get_mol = sys.argv[1]
        scratch = sys.argv[2]
        pop_mol = sys.argv[3]
        mol_file = sys.argv[4]

    if get_mol == 'T':
        print(
            'extract all molecules from reaction systems in '
            'current dir...'
        )
        curr_dir = getcwd()
        done_file = curr_dir+'/done_RS.txt'
        lookup_file = (
            '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
        )
        molecule_dataset = read_molecule_lookup_file(
            lookup_file=lookup_file
        )
        update_molecule_DB(
            rxn_syst.yield_rxn_syst(curr_dir+'/'),
            from_scratch=scratch,
            dataset=molecule_dataset,
            done_file=done_file
        )

    if pop_mol == 'T':
        vdwScale = 0.8
        boxMargin = 4.0
        spacing = 0.5
        N_conformers = 100
        MW_thresh = 500
        print('settings:')
        print('    VDW scale:', vdwScale)
        print('    Box Margin:', boxMargin, 'Angstrom')
        print('    Grid spacing:', spacing, 'Angstrom')
        print('    No Conformers:', N_conformers)
        print('    MW threshold:', MW_thresh, 'g/mol')
        print('    Molecule file:', mol_file)
        # inp = input('happy with these? (T/F)')
        # if inp == 'F':
        #     sys.exit('change them in the source code')
        # elif inp != 'T':
        #     sys.exit('I dont understand, T or F?')
        print(
            'populate the properties attributes for all '
            'molecules in DB...'
        )
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        if mol_file == 'F':
            mol_file = False
        populate_all_molecules(directory=directory,
                               vdwScale=vdwScale,
                               boxMargin=boxMargin,
                               spacing=spacing,
                               N_conformers=N_conformers,
                               MW_thresh=MW_thresh,
                               mol_file=mol_file)

    if input('do plots? (t/f)') == 't':
        directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
        #######
        # molecule distributions
        #######
        # plot synthetic accessibility VS no heavy atoms
        mol_logP_vs_logS(output_dir=directory, plot_suffix='mol_cf')
        mol_logP_vs_XlogP(output_dir=directory, plot_suffix='mol_cf')
        # plot synthetic accessibility VS complexity
        mol_SA_vs_compl(output_dir=directory, plot_suffix='mol_cf')
        # plot synthetic accessibility VS no rotatable bonds
        mol_SA_vs_NRB(output_dir=directory, plot_suffix='mol_cf')
        # # plot synthetic accessibility VS no heavy atoms
        mol_SA_vs_NHA(output_dir=directory, plot_suffix='mol_cf')
        # # plot distributions of attributes
        mol_all_dist(output_dir=directory, plot_suffix='mol_cf')
    sys.exit()


if __name__ == "__main__":
    main()
