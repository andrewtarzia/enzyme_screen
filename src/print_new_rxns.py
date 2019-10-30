#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to print all new reactions in a data set.

Author: Andrew Tarzia

Date Created: 14 Oct 2018

"""
import os
from rdkit.Chem import Descriptors


def check_rxn_unique(reaction_reported, rs):
    """Check (using the sorted list of component molecule weights)
    if a rxn is unique.

    """
    # get list of SMILES of all components
    r_smiles = []
    r_MW = []
    for r in rs.components:
        r_smiles.append(r.SMILES)
        r_MW.append(Descriptors.MolWt(r.mol))
    r_smiles = [x for _, x in sorted(zip(r_MW, r_smiles))]
    if r_smiles in reaction_reported:
        unique = False
    else:
        reaction_reported.append(r_smiles)
        unique = True

    return unique, reaction_reported


def print_new_rxns(output_dir, generator, output_file):
    """Print all new possible and unique reactions that fit.

    """
    reaction_reported = []
    count = 0
    # iterate over reaction system files
    with open(output_file, 'w') as f:
        f.write('rs_pkl,EC,all_fit\n')
        for rs in generator:
            if rs.skip_rxn is True:
                continue
            unique, reaction_reported = check_rxn_unique(
                reaction_reported,
                rs
            )
            if unique is False:
                continue
            if rs.all_fit is True:
                count += 1
                print("------ New Reaction:")
                rs.print_rxn_system()
                f.write(rs.pkl+','+rs.EC+','+str(rs.all_fit)+'\n')

    print("There are", count, "new reactions!")


if __name__ == "__main__":
    import sys
    from rxn_syst import yield_rxn_syst

    if (not len(sys.argv) == 2):
        print('Usage: print_new_rxns.py out_file\n')
        print('   out_file: file name to write output to.')
        sys.exit()
    else:
        out_file = sys.argv[1]

    search_output_dir = os.getcwd()+'/'
    size_thresh = 4.2  # angstroms
    print('settings:')
    print('    Diffusion threshold:', size_thresh, 'Angstrom')
    inp = input('happy with these? (T/F)')
    if inp == 'F':
        sys.exit('change them in the source code')
    elif inp != 'T':
        sys.exit('I dont understand, T or F?')

    # # print new reactions
    print_new_rxns(output_dir=search_output_dir,
                   generator=yield_rxn_syst(search_output_dir),
                   output_file=out_file)

    sys.exit()
