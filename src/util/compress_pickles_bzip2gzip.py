#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to compress pickle files using bzip2.

Author: Andrew Tarzia

Date Created: 08 Oct 2018

"""
import pickle
import bz2
import gzip
import glob
import os
from rxn_syst import reaction
from molecule import molecule

files = glob.glob("*.bpkl")

for file in files:
    if os.path.exists(file.replace('.bpkl', '.gpkl')):
        continue
    print(file)
    # load pickle
    with bz2.BZ2File(file, 'rb') as input:
        obj = pickle.load(input)
    # resave with bzip
    with gzip.GzipFile(file.replace('.bpkl', '.gpkl'), 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    # # reload with bzip
    # with bz2.BZ2File(file.replace('.pkl', '.bpkl'), 'rb') as output:
    #     obj2 = pickle.load(output)
    # print(obj.__dict__)
    # print(obj2.__dict__)
    # print(obj.__dict__ == obj2.__dict__)
    # break
