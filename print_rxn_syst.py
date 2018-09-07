#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to print reaction system information.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
# ensure cpickle usage
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import sys
import rxn_syst

if (not len(sys.argv) == 2):
    print('Usage: print_rxn_syst.py file\n')
    print('    file: pickle file containing reaction system\n')
    sys.exit()
else:
    file = sys.argv[1]

with open(file, 'rb') as input:
    rs = pickle.load(input)

rs.print_rxn_system()
