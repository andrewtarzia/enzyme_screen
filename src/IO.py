#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import pickle
import gzip
from os.path import exists, join
import pandas as pd

import molecule


def done_list_read(directory, file_name='collected_mols.txt'):
    """File that contains the names of molecules that have already been
    collected and their associated pkl file names.

    Returns the list.

    This function is not currently used, but could replace the molecule
    search functions.

    """
    names = []
    pkls = []
    with open(directory+file_name, 'r') as f:
        for line in f.readlines():
            names.append(line.rstrip().split('___')[0])
            pkls.append(line.rstrip().split('___')[1])
    return names, pkls


def done_list_write(
    new_name,
    pkl,
    directory,
    file_name='collected_mols.txt'
):
    """
    Appends (or writes) file with list of failed names.

    This function is not currently used, but could replace the molecule
    search functions.

    """
    if not exists(join(directory, file_name)):
        with open(join(directory, file_name), 'w') as f:
            f.write('\n')

    with open(join(directory, file_name), 'a') as f:
        f.write(new_name+"___"+pkl+'\n')


def fail_list_read(directory, file_name='failures.txt'):
    """
    File that contains the names of molecules that failed resolution to
    avoid double checking.

    Returns the list.

    """

    if not exists(join(directory, file_name)):
        with open(join(directory, file_name), 'w') as f:
            f.write('\n')

    names = []
    with open(join(directory, file_name), 'r') as f:
        for line in f.readlines():
            if line.rstrip != '':
                names.append(line.rstrip())
    return names


def fail_list_write(new_name, directory, file_name='failures.txt'):
    """
    Appends (or writes) file with list of failed names.

    """
    if not exists(join(directory, file_name)):
        with open(directory+file_name, 'w') as f:
            f.write('\n')

    with open(join(directory, file_name), 'a') as f:
        f.write(new_name+'\n')


def change_all_pkl_suffixes(directory):
    """
    Change the suffixes of pkl file names in all molecules in
    directory.

    For Debugging

    """
    for i in molecule.yield_molecules(directory=directory):
        i.pkl = i.pkl.replace('.pkl', '.gpkl')
        i.pkl = i.pkl.replace('.bpkl', '.gpkl')
        new_rs_pkls = []
        try:
            for j in i.rs_pkls:
                new_rs_pkls.append(
                    j.replace('.pkl', '.gpkl').replace(
                        '.bpkl', '.gpkl'
                    )
                )
        except AttributeError:
            pass
        i.rs_pkls = new_rs_pkls
        i.save_object(i.pkl)


def read_molecule_lookup_file(lookup_file):
    """
    Utility function the returns Pandas DataFrame of molecule lookups.

    """
    dataset = pd.read_table(lookup_file, delimiter='=', skiprows=[0],
                            names=['SMILES', 'iupac', 'name',
                                   'DB', 'DB_ID', 'KEGG_ID',
                                   'CHEBI_ID', 'InChiKey', 'pkl'],
                            engine='python')
    return dataset


def update_lookup_files_text():
    """
    Utility function that updates the KEGG translator and molecule
    search file.

    Uses text files and rewrites each update -- is SLOW!
    """
    print('updating lookup files')
    # lookup file
    lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
    with open(lookup_file, 'w') as f:
        f.write('SMILES___iupac___name___DB___DB_ID___KEGG_ID')
        f.write('___InChiKey___CHEBI_ID___pkl\n')

    # KEGG translator
    translator = '/home/atarzia/psp/molecule_DBs/KEGG/translator.txt'
    with open(translator, 'w') as f:
        f.write('')
    # iterate over all molecules in DB and if they have a KEGG ID then
    # write translation AND write information to lookup file
    directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
    for mol in molecule.yield_molecules(directory=directory):
        # KEGG translator
        if 'KEGG' in mol.DB_list:
            KID = mol.KEGG_ID
            pkl = mol.pkl
            with open(translator, 'a') as f:
                f.write(KID+'__'+pkl+'\n')
        # lookup file
        smiles = '-'
        iupac = '-'
        name = '-'
        DB = '-'
        DB_ID = '-'
        KEGG_ID = '-'
        CHEBI_ID = '-'
        IKEY = '-'
        pkl = '-'
        if mol.SMILES is not None:
            smiles = mol.SMILES
        if mol.iupac_name is not None:
            if type(mol.iupac_name) != list:
                iupac = mol.iupac_name
        if mol.name is not None:
            name = mol.name
        if mol.DB is not None:
            DB = mol.DB
        if mol.DB_ID is not None:
            DB_ID = mol.DB_ID
        if 'KEGG' in mol.DB_list:
            KEGG_ID = mol.KEGG_ID
        if mol.chebiID is not None:
            CHEBI_ID = mol.chebiID
        try:
            if mol.InChIKey is not None:
                IKEY = mol.InChIKey
        except AttributeError:
            pass
        if mol.pkl is not None:
            pkl = mol.pkl
        with open(lookup_file, 'a') as f:
            f.write(smiles+'___'+iupac+'___'+name+'___')
            f.write(DB+'___'+DB_ID+'___'+KEGG_ID+'___'+IKEY+'___')
            f.write(CHEBI_ID+'___'+pkl+'\n')


def write_lookup_row(mol):
    """
    Write row for lookup file and return as Pandas Dataframe.

    """
    # lookup file
    smiles = '-'
    iupac = '-'
    name = '-'
    DB = '-'
    DB_ID = '-'
    KEGG_ID = '-'
    CHEBI_ID = '-'
    IKEY = '-'
    pkl = '-'
    if mol.SMILES is not None:
        smiles = mol.SMILES
    if mol.iupac_name is not None:
        if type(mol.iupac_name) != list:
            iupac = mol.iupac_name
    if mol.name is not None:
        name = mol.name
    if mol.DB is not None:
        DB = mol.DB
    if mol.DB_ID is not None:
        DB_ID = mol.DB_ID
    if 'KEGG' in mol.DB_list:
        KEGG_ID = mol.KEGG_ID
    try:
        if mol.chebiID is not None:
            if isinstance(mol.chebiID, list):
                CHEBI_ID = ' '.join(mol.chebiID)
            else:
                CHEBI_ID = mol.chebiID
    except AttributeError:
        pass
    try:
        if mol.InChiKey is not None:
            IKEY = mol.InChiKey
    except AttributeError:
        pass
    if mol.pkl is not None:
        pkl = mol.pkl
    ROW_DF = pd.DataFrame(
        {
            'SMILES': smiles, 'iupac': iupac,
            'name': name, 'DB': DB, 'DB_ID': DB_ID,
            'KEGG_ID': KEGG_ID, 'CHEBI_ID': CHEBI_ID,
            'InChiKey': IKEY, 'pkl': pkl
        },
        index=[0]
    )
    return ROW_DF


def write_lookup_files(lookup_file, translator, directory):
    """
    Utility function that writes lookup files from scratch.

    """
    print('writing lookup files...')
    translation = {}
    lookup = pd.DataFrame(
        columns=[
            'SMILES', 'iupac', 'name', 'DB', 'DB_ID',
            'KEGG_ID', 'CHEBI_ID', 'InChiKey', 'pkl'
        ]
    )
    # iterate over all molecules in DB and if they have a KEGG ID then
    # write translation AND write information to lookup file
    for mol in molecule.yield_molecules(directory=directory):
        # KEGG translator
        if 'KEGG' in mol.DB_list:
            KID = mol.KEGG_ID
            pkl = mol.pkl
            write_translation_line(KID, pkl, translation)
        ROW_DF = write_lookup_row(mol)
        lookup = lookup.append(ROW_DF, ignore_index=True)
    # write translator
    with gzip.GzipFile(translator, 'wb') as output:
        pickle.dump(translation, output, pickle.HIGHEST_PROTOCOL)
    # write lookup dataframe
    lookup.to_csv(lookup_file, index=False, sep='=')
    print('done!')


def update_lookup_files(mol, unique):
    """
    Utility function that updates the KEGG translator and molecule
    search file (lookup.txt) with a new molecule.

    """
    print('updating lookup files...')
    # lookup file
    lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
    # KEGG translator
    translator = '/home/atarzia/psp/molecule_DBs/KEGG/translator.txt'
    # molecule DB directory
    directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
    if exists(lookup_file) is False or exists(translator) is False:
        # need to write from scratch
        write_lookup_files(lookup_file, translator, directory)
    # read files
    # read translator
    with gzip.GzipFile(translator, 'rb') as output:
        translation = pickle.load(output)
    # read lookup
    molecule_dataset = read_molecule_lookup_file(
        lookup_file=lookup_file
    )
    # molecule_dataset == lookup in other functions
    if unique is True:
        # need to append a new row to the existing lookup objects
        # KEGG translator
        if 'KEGG' in mol.DB_list:
            KID = mol.KEGG_ID
            pkl = mol.pkl
            translation[KID] = pkl
        # update lookup file
        ROW_DF = write_lookup_row(mol)
        molecule_dataset = molecule_dataset.append(
            ROW_DF,
            ignore_index=True
        )
    else:
        # may need to modify the lookup file (the KEGG translator
        # should not
        # change) -- add check if it does
        if 'KEGG' in mol.DB_list:
            KID = mol.KEGG_ID
            pkl = mol.pkl
            # update translation file
            try:
                print(translation[KID])
            except KeyError:
                write_translation_line(KID, pkl, translation)
            # within this if statement are some checks for
            # inconsistencies in order to fix them
            if KID != '':
                if ' ' in KID:
                    for i in KID.split(' '):
                        if translation[i] != pkl:
                            print(
                                KID, i, translation[i], pkl, mol.name
                            )
                            import sys
                            sys.exit(
                                '1a - there is a problem here '
                                '-- KEGG translation has changed'
                            )
                        for key in translation:
                            val = translation[key]
                            if val == pkl and key != i:
                                print(
                                    KID, i, translation[i],
                                    pkl, mol.name
                                )
                                import sys
                                sys.exit(
                                    '2a - there is a problem here'
                                    ' -- KEGG translation has changed'
                                )
                        break
                else:
                    if translation[KID] != pkl:
                        print(KID, translation[KID], pkl, mol.name)
                        import sys
                        sys.exit(
                            '1b - there is a problem here'
                            ' -- KEGG translation has changed'
                        )
                    for key in translation:
                        val = translation[key]
                        if val == pkl and key != KID:
                            print(KID, translation[KID], pkl, mol.name)
                            import sys
                            sys.exit(
                                '2b - there is a problem here'
                                ' -- KEGG translation has changed'
                            )
        # modify lookup file row associated with mol
        matching_pkl = molecule_dataset[
            molecule_dataset.pkl == mol.pkl
        ]
        if len(matching_pkl) == 0:
            import sys
            sys.exit(
                'problem here with pkl in molecule DB -- what is it?'
            )
        for idx, row in matching_pkl.iterrows():
            if row['pkl'] == mol.pkl:
                if row.SMILES == '-' and mol.SMILES is not None:
                    print('changing SMILES')
                    row.SMILES = mol.SMILES
                if row.DB == '-' and mol.DB is not None:
                    print('changing DB')
                    row.DB = mol.DB
                if row.iupac == '-' and mol.iupac_name is not None:
                    print('changing iupac')
                    row.iupac = mol.iupac_name
                if row.name == '-' and mol.name is not None:
                    print('changing name')
                    row.name = mol.name
                if row.DB_ID == '-' and mol.DB_ID is not None:
                    print('changing DB ID')
                    row.DB_ID = mol.DB_ID
                try:
                    if row.KEGG_ID == '-' and mol.KEGG_ID is not None:
                        print('changing KEGG')
                        row.KEGG_ID = mol.KEGG_ID
                except AttributeError:
                    pass
                try:
                    test1 = row.CHEBI_ID == '-'
                    if test1 and mol.CHEBI_ID is not None:
                        print('changing CHEBI')
                        row.CHEBI_ID = mol.CHEBI_ID
                except AttributeError:
                    pass
                try:
                    test1 = row.InChiKey == '-'
                    if test1 and mol.InChiKey is not None:
                        print('changing IKEY')
                        row.InChiKey = mol.InChiKey
                except AttributeError:
                    pass
                molecule_dataset.iloc[idx] = row
    # write translator
    with gzip.GzipFile(translator, 'wb') as output:
        pickle.dump(translation, output, pickle.HIGHEST_PROTOCOL)
    # write lookup dataframe
    molecule_dataset.to_csv(lookup_file, index=False, sep='=')


def read_molecule_list(file):

    mol_list = []
    with open(file, 'r') as f:
        for line in f.readlines():
            mol_list.append(line.rstrip())

    return mol_list
