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

    def get_rxn_system(self):
        raise NotImplementedError()

    def fail_conversion(self, comp):
        print('could not be converted to MOL.')
        print('>>> comp:', comp)
        self.skip_rxn = True
        self.skip_reason = 'Component could not be converted to MOL'

    def fail_generic(self, comp):
        print('DB gave generic structure.')
        print('>>> comp:', comp)
        self.skip_rxn = True
        self.skip_reason = 'Rxn includes generic structure'

    def fail_polymer(self, comp):
        print('(n) in component implies polymeric species')
        print('>>> comp:', comp)
        self.skip_rxn = True
        self.skip_reason = 'Rxn includes polymeric species'

    def save_reaction(self, filename):
        """
        Pickle reaction system object to file.

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

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(EC={self.EC}, DB={self.DB}, '
            f'ID={self.DB_ID}, pkl={self.pkl})'
        )

    def __repr__(self):
        return str(self)


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


def get_RS(filename, output_dir, verbose=False):
    """
    Read in reaction system from filename.

    """
    _rsf = filename.replace(output_dir+'sRS-', '').replace('.gpkl', '')
    EC_, DB, DB_ID = _rsf.split('-')
    EC = EC_.replace("_", ".").replace('XX', '-')
    rs = Reaction(EC, DB, DB_ID)
    if exists(output_dir+rs.pkl) is False:
        raise FileNotFoundError(
            'you have not collected all reaction systems.'
            f'{output_dir+rs.pkl} does not exist!'
        )
    # load in rxn system
    if verbose:
        print(f'loading: {rs.pkl}')
    rs = rs.load_object(output_dir+rs.pkl, verbose=False)
    return rs
