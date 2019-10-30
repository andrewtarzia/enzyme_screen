#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the reaction system class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import pickle
import gzip


class Reaction:
    """
    Class that defines a reaction system for all databases.

    """

    def __init__(self, EC, DB, DB_ID):
        # all non DB unique properties
        self.EC = EC
        self.DB = DB
        self.DB_ID = DB_ID
        # for unknwon EC tiers (given by '-'), use a known delimeter.
        EC_ul = EC.replace('.', '_').replace('-', 'XX')
        self.pkl = 'sRS-'+EC_ul+'-'+str(DB)+'-'+str(DB_ID)+'.gpkl'
        # need to have this as None by default
        self.UniprotID = None
        # allows for noting of skipped reaction
        self.skip_rxn = False
        # quick comment on why a rxn is skipped
        self.skip_reason = None
        # molecular components
        self.components = None
        self.mol_collected = False

    def save_object(self, filename):
        """Pickle reaction system object to file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        # Overwrites any existing file.
        with gzip.GzipFile(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(self, filename, verbose=True):
        """unPickle reaction system object from file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        if verbose:
            print('loading:', filename)
        with gzip.GzipFile(filename, 'rb') as input:
            self = pickle.load(input)
            return self


def yield_rxn_syst(output_dir, verbose=False):
    """Iterate over reaction systems for analysis.

    """
    react_syst_files = sorted(glob.glob(output_dir+'sRS-*.gpkl'))
    for rsf in react_syst_files:
        # try:
        rs = get_RS(
            filename=rsf,
            output_dir=output_dir,
            verbose=verbose
        )
        # except:
        #     print(
        #         'error loading:'
        #     )
        #     print(rsf)
        #     sys.exit()
        yield rs


def yield_rxn_syst_filelist(output_dir, filelist, verbose=False):
    """Iterate over reaction systems for analysis - uses a file list.

    """
    react_syst_files = []
    with open(filelist, 'r') as f:
        for line in f.readlines():
            react_syst_files.append(line.strip())
    for rsf in react_syst_files:
        # try:
        rs = get_RS(
            filename=rsf,
            output_dir=output_dir,
            verbose=verbose
        )
        # except:
        #     print('error loading:')
        #     print(rsf)
        #     sys.exit()
        yield rs


def change_all_pkl_suffixes_RS(directory):
    """For debugging. Changes all pkl attributes to be .gpkl

    """
    for rs in yield_rxn_syst(output_dir=directory):
        rs.pkl = rs.pkl.replace('.bpkl', '.gpkl')
        rs.save_object(directory+rs.pkl)
