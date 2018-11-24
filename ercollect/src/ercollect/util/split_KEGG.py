#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to convert KEGG JSON file into a JSON file with EC numbers at the top
level.

Author: Andrew Tarzia

Date Created: 24 Nov 2018

"""
import json
rxn_DB_file = 'br08201.json'
with open(rxn_DB_file, 'r') as data_file:
        rxn_DB = json.load(data_file)

all_EC_dict = {}
for i in rxn_DB['children']:
    for j in i['children']:
        for k in j['children']:
            for m in k['children']:
                # print(m['name'])
                try:
                    # print(m['children'])
                    all_EC_dict[m['name']] = m['children']
                except KeyError:
                    pass

# save all_EC_dict
with open('br08201_ECtop.json', 'w') as fp:
    json.dump(all_EC_dict, fp)
