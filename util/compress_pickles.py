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
import glob

files = glob.glob("*.pkl")

for file in files:
    print(file)
    # load pickle
    with open(file, 'rb') as input:
        obj = pickle.load(input)
    # resave with bzip
    with bz2.BZ2File(file.replace('.pkl', '.bpkl'), 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    # # reload with bzip
    # with bz2.BZ2File(file.replace('.pkl', '.bpkl'), 'rb') as output:
    #     obj2 = pickle.load(output)
    # print(obj.__dict__)
    # print(obj2.__dict__)
    # print(obj.__dict__ == obj2.__dict__)
    # break
