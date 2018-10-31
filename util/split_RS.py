#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to split moleucle pkl files into NP lists for trivial parallelisation.

Author: Andrew Tarzia

Date Created: 06 Oct 2018

"""
import glob
import os

NP = 2

curr_dir = os.getcwd()+'/'

files = sorted(glob.glob("*.gpkl"))

print(len(files), 'files --', len(files)/NP, 'per process')

out_files = ['file_list_'+str(i+1)+'.txt' for i in range(NP)]
print(out_files)

files_sep = {}
# separate all files into different processors
for n in range(NP):
    files_sep[str(n+1)] = []

n = 1
for m in files:
    files_sep[str(n)].append(curr_dir+m)
    if n == NP:
        n = 1
    else:
        n += 1

# write files list to file
for n in range(NP):
    print(len(files_sep[str(n+1)]))
    with open(out_files[n], 'w') as f:
        for i in files_sep[str(n+1)]:
            f.write(i+'\n')
