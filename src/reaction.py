#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the reaction system class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os.path import exists
import glob
import pickle
import gzip


class Reaction:
    """
    Class that defines a reaction system for all databases.

    """

    def __init__(self, EC, DB, DB_ID, params):
        # all non DB unique properties
        self.EC = EC
        self.DB = DB
        self.DB_ID = DB_ID
        self.params = params
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

    def load_reaction(self, filename, verbose=True):
        """
        unPickle reaction system object from file.

        """
        filename = filename.replace('.pkl', '.gpkl')
        filename = filename.replace('.bpkl', '.gpkl')
        if verbose:
            print('loading:', filename)
        with gzip.GzipFile(filename, 'rb') as input:
            self = pickle.load(input)
            return self

        )


def yield_rxn_syst(output_dir, filelist=None, verbose=False):
    """
    Iterate over reaction systems for analysis.

    """
    if filelist is None:
        react_syst_files = sorted(glob.glob(output_dir+'sRS-*.gpkl'))
    else:
        react_syst_files = []
        with open(filelist, 'r') as f:
            for line in f.readlines():
                react_syst_files.append(line.strip())
    for rsf in react_syst_files:
        rs = get_RS(
            filename=rsf,
            output_dir=output_dir,
            verbose=verbose
        )
        yield rs



    """
