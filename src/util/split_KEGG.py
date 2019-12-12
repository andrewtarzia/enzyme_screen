#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to convert KEGG JSON file into a file with EC numbers at the top
level.

Author: Andrew Tarzia

Date Created: 24 Nov 2018

"""

import sys
import json
from os.path import exists


def main():
    if (not len(sys.argv) == 2):
        raise ValueError('Usage: split_KEGG.py JSON')

    json_file = sys.argv[1]
    if exists(json_file):
        rxn_DB_file = json_file
    else:
        raise FileExistsError('Download KEGG BRITE JSON file')

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

    EC_list = all_EC_dict.keys()

    # save all_EC_dict
    with open(json_file.replace('.json', '_ECtop.json'), 'w') as fp:
        json.dump(all_EC_dict, fp)
    with open(json_file.replace('.json', '_EClist.txt'), 'w') as f:
        for line in EC_list:
            f.write(f'{line}\n')


if __name__ == '__main__':
    main()
