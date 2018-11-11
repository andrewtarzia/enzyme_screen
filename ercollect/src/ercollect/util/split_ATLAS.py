#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to split ATLAS FULL CSV into a CSV for each top tier EC number.
Must be run in directory with ATLAS DB.

Author: Andrew Tarzia

Date Created: 02 Oct 2018

"""
import pandas as pd

atlas = pd.read_csv('ATLAS-FULL.csv')
atlas_1 = pd.DataFrame(columns=atlas.columns)
atlas_2 = pd.DataFrame(columns=atlas.columns)
atlas_3 = pd.DataFrame(columns=atlas.columns)
atlas_4 = pd.DataFrame(columns=atlas.columns)
atlas_5 = pd.DataFrame(columns=atlas.columns)
atlas_6 = pd.DataFrame(columns=atlas.columns)
for idx, row in atlas.iterrows():
    ECs = row['REACTIONRULE'].split('|')
#    print(ECs)
    top_tiers = [i.split('.')[0].lstrip().rstrip() for i in ECs]
#    print(top_tiers)
    top_tiers = list(set(top_tiers))
    print(top_tiers)
    for i in top_tiers:
        if i == '1':
            atlas_1 = atlas_1.append(row)
        if i == '2':
            atlas_2 = atlas_2.append(row)
        if i == '3':
            atlas_3 = atlas_3.append(row)
        if i == '4':
            atlas_4 = atlas_4.append(row)
        if i == '5':
            atlas_5 = atlas_5.append(row)
        if i == '6':
            atlas_6 = atlas_6.append(row)
print(len(atlas_1), len(atlas_2), len(atlas_3), len(atlas_4), len(atlas_5),
      len(atlas_6))
atlas_1.to_csv('ATLAS-1.csv', index=False)
atlas_2.to_csv('ATLAS-2.csv', index=False)
atlas_3.to_csv('ATLAS-3.csv', index=False)
atlas_4.to_csv('ATLAS-4.csv', index=False)
atlas_5.to_csv('ATLAS-5.csv', index=False)
atlas_6.to_csv('ATLAS-6.csv', index=False)
