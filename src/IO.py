#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining useful IO functions.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os.path import exists, join


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


def read_molecule_list(file):

    mol_list = []
    with open(file, 'r') as f:
        for line in f.readlines():
            mol_list.append(line.rstrip())

    return mol_list
