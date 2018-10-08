#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to print reaction system information.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import pickle
import bz2
import sys

if (not len(sys.argv) == 2):
    print('Usage: print_rxn_syst.py file\n')
    print('    file: compressed (.bpkl) pickled reaction system\n')
    sys.exit()
else:
    file = sys.argv[1]

with bz2.BZ2File(file, 'rb') as input:
    rs = pickle.load(input)

rs.print_rxn_system()

try:
    print('all fit?', rs.all_fit)
    print('largest component size?', rs.max_comp_size)
except AttributeError:
    pass
