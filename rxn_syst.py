#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the reaction system class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
# ensure cpickle usage
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


class reaction:
    """Class that defines a reaction system for all databases.

    """

    def __init__(self, EC, DB, DB_ID):
        # all non DB unique properties
        self.EC = EC
        self.DB = DB
        self.DB_ID = DB_ID
        EC_ul = EC.replace('.', '_')
        self.pkl = 'sRS-'+EC_ul+'-'+str(DB)+'-'+str(DB_ID)+'.pkl'
        self.skip_rxn = False  # allows for noting of skipped reaction
        self.components = None  # molecular components
        self.all_fit = None  # do all the components fit?
        # None implies that it is ambiguous/unknown.
        self.seed_MOF = None  # will the protein sequence seed MOF growth
        self.req_mod = None  # does it require modification to seed MOF growth

    def check_all_fit(self, threshold, molecule_output):
        """Check if all components of reaction system fit.

        """
        all_fit = True
        max_comp_size = 0
        for r in self.components:
            # is molecule in molecule_output?
            if r.name in list(molecule_output['name']):
                r_diam = float(molecule_output[molecule_output['name'] == r.name]['mid_diam'])
                max_comp_size = max([r_diam, max_comp_size])
                if r_diam > threshold or r_diam == 0:
                    all_fit = False
            else:
                # molecule not in output - assume it wasn't in database
                # dont report!
                all_fit = False
                break

        if all_fit is True:
            self.all_fit = True
            self.max_comp_size = max_comp_size
            self.print_rxn_system()
        else:
            self.all_fit = False
            self.max_comp_size = max_comp_size

    def save_object(self, filename):
        """Pickle reaction system object to file.

        """
        # Overwrites any existing file.
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(self, filename, verbose=True):
        """unPickle reaction system object from file.

        """
        if verbose:
            print('loading:', filename)
        with open(filename, 'rb') as input:
            self = pickle.load(input)
            return self

    def print_rxn_system(self):
        """Fancy print of reaction system.

        """
        print('-----------------------------------')
        print('EC:', self.EC)
        print('Database:', self.DB)
        print('Database ID:', self.DB_ID)
        print('-----------------------------------')
        if self.components is not None:
            for i in self.components:
                print(i.name, ' (ID:', i.cID+') as', i.role)
            print('-----------------------------------')
        if self.all_fit is True:
            print('All components will diffuse through!')
            print('-----------------------------------')


def get_reaction_systems(EC, DB, output_dir):
    """Get reaction system from SABIO reaction ID (rID).

    """
    if DB == 'SABIO':
        from SABIO_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir)
    elif DB == 'KEGG':
        from KEGG_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir)
    elif DB == 'BKMS':
        from BKMS_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir)
    elif DB == 'BRENDA':
        from BRENDA_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir)
