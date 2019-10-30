#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to print reaction system information.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import pickle
import gzip
import sys

if (not len(sys.argv) == 2):
    print('Usage: print_rxn_syst.py file\n')
    print('    file: compressed (.gpkl) pickled reaction system\n')
    sys.exit()
else:
    file = sys.argv[1]

with gzip.GzipFile(file, 'rb') as input:
    rs = pickle.load(input)

rs.print_rxn_system()

try:
    print('all fit?', rs.all_fit)
    print('largest component size?', rs.max_comp_size)
except AttributeError:
    pass
