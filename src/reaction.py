#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the reaction system class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os.path import exists, join
import glob
import json
import pickle
import gzip
import pandas as pd


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

    def fail_generic_smiles(self, comp):
        print('DB gave SMILEs with generic (*).')
        print('>>> comp:', comp)
        self.skip_rxn = True
        self.skip_reason = 'Rxn includes generic SMILEs'

    def fail_fail_list(self, comp):
        print('DB gave molecule that already failed.')
        print('>>> comp:', comp)
        self.skip_rxn = True
        self.skip_reason = 'Rxn includes failure'

    def fail_polymer(self, comp):
        print('(n) in component implies polymeric species')
        print('>>> comp:', comp)
        self.skip_rxn = True
        self.skip_reason = 'Rxn includes polymeric species'

    def fail_properties(self):
        print('no components above NHA threshold.')
        self.skip_rxn = True
        self.skip_reason = 'No components above NHA threshold'

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

    def get_component_properties(self):
        """
        Get properties of all components that define reaction system.

        """

        self.get_max_mid_diameter()
        self.get_logPs()
        self.get_logSs()
        self.get_SAs()

    def get_max_mid_diameter(self):
        """
        Get diffusion properties of all components.

        """
        max_min_mid_diam = 0

        for m in self.components:
            print(m)
            name = m.name
            diam_file = join(
                self.params['molec_dir'],
                name+'_size.csv'
            )
            results = pd.read_csv(diam_file)
            if len(results) == 0:
                continue
            min_mid_diam = min(results['diam2'])
            max_min_mid_diam = max([min_mid_diam, max_min_mid_diam])

        self.max_min_mid_diam = max_min_mid_diam

    def get_logPs(self):
        """
        Get logP (hydrophobicity) properties of all components.

        """

        self.min_logP = 1E10
        self.max_logP = -1E10
        for m in self.components:
            print(m)
            name = m.name
            prop_file = join(
                self.params['molec_dir'],
                name+'_prop.json'
            )

            if not exists(prop_file):
                raise FileNotFoundError(f'{prop_file} not found')

            with open(prop_file, 'r') as f:
                prop_dict = json.load(f)

            print(prop_dict)

            # Only include molecules with more than a certain number of
            # heavy atoms.
            if prop_dict['NHA'] < self.params['NHA_thresh']:
                continue

            self.min_logP = min([self.min_logP, prop_dict['logP']])
            self.max_logP = max([self.max_logP, prop_dict['logP']])

    def get_logSs(self):
        """
        Get logS (water solubility) properties of all components.

        """

        self.min_logS = 1E10
        self.max_logS = -1E10
        for m in self.components:
            print(m)
            name = m.name
            prop_file = join(
                self.params['molec_dir'],
                name+'_prop.json'
            )

            if not exists(prop_file):
                raise FileNotFoundError(f'{prop_file} not found')

            with open(prop_file, 'r') as f:
                prop_dict = json.load(f)

            print(prop_dict)

            # Only include molecules with more than a certain number of
            # heavy atoms.
            if prop_dict['NHA'] < self.params['NHA_thresh']:
                continue

            self.min_logS = min([self.min_logS, prop_dict['logS']])
            self.max_logS = max([self.max_logS, prop_dict['logS']])

    def get_SAs(self):
        """
        Get synthetic accessibility properties of all components.

        """
        self.r_max_SA = 0
        self.p_max_SA = 0
        self.delta_SA = 0
        for m in self.components:
            print(m)
            name = m.name
            prop_file = join(
                self.params['molec_dir'],
                name+'_prop.json'
            )

            if not exists(prop_file):
                raise FileNotFoundError(f'{prop_file} not found')

            with open(prop_file, 'r') as f:
                prop_dict = json.load(f)

            print(prop_dict)
            # Only include molecules with more than a certain number of
            # heavy atoms.
            if prop_dict['NHA'] < self.params['NHA_thresh']:
                continue

            if m.role == 'reactant':
                self.r_max_SA = max([
                    self.r_max_SA,
                    prop_dict['Synth_score']
                ])
            elif m.role == 'product':
                self.p_max_SA = max([
                    self.p_max_SA,
                    prop_dict['Synth_score']
                ])

        self.delta_SA = self.p_max_SA - self.r_max_SA

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(EC={self.EC}, DB={self.DB}, '
            f'ID={self.DB_ID}, pkl={self.pkl})'
        )

    def __repr__(self):
        return str(self)


def yield_rxn_syst(output_dir, pars, file=None, verbose=False):
    """
    Iterate over reaction systems for analysis.

    """
    if file is None:
        react_syst_files = sorted(glob.glob(join(
            output_dir,
            'sRS-*.gpkl'
        )))
    else:
        react_syst_files = []
        with open(file, 'r') as f:
            for line in f.readlines():
                if exists(line.strip()):
                    react_syst_files.append(line.strip())
    for rsf in react_syst_files:
        rs = get_RS(
            filename=rsf,
            output_dir=output_dir,
            pars=pars,
            verbose=verbose
        )
        yield len(react_syst_files), rs


def get_RS(filename, output_dir, pars, verbose=False):
    """
    Read in reaction system from filename.

    """

    prefix = join(output_dir, 'sRS-')
    _rsf = filename.replace(prefix, '').replace('.gpkl', '')
    EC_, DB, DB_ID = _rsf.split('-')
    EC = EC_.replace("_", ".").replace('XX', '-')
    rs = Reaction(EC, DB, DB_ID, params=pars)
    abs_filename = join(output_dir, rs.pkl)
    if exists(abs_filename) is False:
        raise FileNotFoundError(
            'you have not collected all reaction systems.'
            f'{abs_filename} does not exist!'
        )
    # load in rxn system
    if verbose:
        print(f'loading: {rs.pkl}')
    rs = rs.load_reaction(abs_filename, verbose=False)
    return rs
