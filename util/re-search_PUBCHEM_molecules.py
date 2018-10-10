#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script used on 10/10/18 to re-search PUBCHEM for all molecules with multiple
entries for IUPAC, logP and complexity.

Author: Andrew Tarzia

Date Created: 10 Oct 2018

"""
from PUBCHEM_IO import get_logP_from_name, get_complexity_from_name, get_IUPAC_from_name
from molecule import yield_molecules, load_molecule
directory = '/home/atarzia/psp/molecule_DBs/atarzia/'

for j, i in enumerate(yield_molecules(directory=directory)):
    print('doing -----', j)
    print(i.name, i.XlogP, i.complexity, i.iupac_name)
    if '\n' in str(i.XlogP) or '\n' in str(i.complexity) or '\n' in str(i.iupac_name):
        print(i.name)
        print('old logp', i.XlogP)
        print('old comp', i.complexity)
        new_iupac = get_IUPAC_from_name(i.name)
        print('old_iupac', i.iupac_name)
        print('new iupac', new_iupac)
        i.iupac_name = new_iupac
        print('get values')
        new_logP = get_logP_from_name(i.name)
        new_comp = get_complexity_from_name(i.name)
        # try with iupac if those dont work
        if new_logP == 'not found' or new_logP is None:
            if i.iupac_name is not None:
                new_logP = get_logP_from_name(i.iupac_name)
        if new_comp == 'not found' or new_comp is None:
            if i.iupac_name is not None:
                new_comp = get_complexity_from_name(i.iupac_name)
        print('new logp', new_logP)
        print('new comp', new_comp)
        print(i.pkl)
        if type(new_logP) is float and type(new_comp) is float:
            print('save!')
            i.XlogP = new_logP
            i.complexity = new_comp
            i.save_object(i.pkl)

# SMILES:
from PUBCHEM_IO import get_logP_from_name, get_complexity_from_name, get_IUPAC_from_name
from molecule import yield_molecules, load_molecule
directory = '/home/atarzia/psp/molecule_DBs/atarzia/'

caught = []
for j, i in enumerate(yield_molecules(directory=directory)):
    print('doing -----', j)
    if '\n' in str(i.SMILES):
        print(i.name)
        print('SMILES', i.SMILES)
        caught.append(i)
len(caught)
for ca in caught:
    print(ca.name, ca.pkl)
    print(ca.SMILES)
    splits = ca.SMILES.split('\n')
    print(splits)
    break
