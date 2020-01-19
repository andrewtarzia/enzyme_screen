#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to visualise the properties of a reaction system.

Author: Andrew Tarzia

Date Created: 19 Jan 2020

"""

import sys
import os

from reaction import get_RS
import utilities


def main():
    if (not len(sys.argv) == 3):
        print("""
Usage: visualise_reaction_system.py param_file file

    param_file : (str)

    file : (str)
        pkl file containing reaction system
""")
        sys.exit()
    else:
        params = utilities.read_params(sys.argv[1])
        file = sys.argv[2]

    file = os.path.join(os.getcwd(), file)
    print(file)

    # Read in reaction system.
    rs = get_RS(
        filename=file,
        output_dir=os.getcwd(),
        pars=params,
        verbose=True
    )
    print(rs)
    if rs.skip_rxn:
        print(f'>>> {rs.skip_reason}')

    print(f'max d = {rs.max_min_mid_diam}\n')

    # Output molecular components and their properties.
    for rsc in rs.components:
        print(rsc)
        print(f'SMILEs = {rsc.SMILES}')
        print(f'logP = {rsc.logP}')
        print(f'logS = {rsc.logS}')
        print(f'SA = {rsc.Synth_score}')


if __name__ == "__main__":
    main()
