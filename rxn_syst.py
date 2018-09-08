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
import glob
import os


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
        self.UniprotID = None  # need to have this as None by default
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
            # use SMILES as check
            if r.SMILES in list(molecule_output['SMILE']):
                mol_frame = molecule_output[molecule_output['SMILE'] == r.SMILES]
                # are they multiple of the same molecule?
                if len(mol_frame) > 1:
                    # this was caused by tautomers of NADH
                    # can we discern by DB
                    DB_frame = mol_frame[mol_frame['DB'] == self.DB]
                    r_diam = float(DB_frame['mid_diam'])
                else:
                    r_diam = float(molecule_output[molecule_output['SMILE'] == r.SMILES]['mid_diam'])
                max_comp_size = max([r_diam, max_comp_size])
                r.mid_diam = r_diam
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
                print(i.name, ' (ID:', i.DB_ID+') as', i.role)
                print('SMILES:', i.SMILES)
            print('-----------------------------------')
        if self.all_fit is True:
            print('All components will diffuse through!')
            print('-----------------------------------')


def get_reaction_systems(EC, DB, output_dir, clean_system=False):
    """Get reaction system from SABIO reaction ID (rID).

    Keywords:
        EC (str) - Enzyme commision number (X.X.X.X)
        DB (str) - name of Database
        output_dir (str) - directory where all data should be saved
        clean_system (bool) - wipe the data in reaction systems for fresh start
            default = False

    """
    if DB == 'SABIO':
        from SABIO_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system)
    elif DB == 'KEGG':
        from KEGG_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system)
    elif DB == 'BKMS':
        from BKMS_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system)
    elif DB == 'BRENDA':
        from BRENDA_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system)


def yield_rxn_syst(output_dir):
    """Iterate over reaction systems for analysis.

    """
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
    for rsf in react_syst_files:
        _rsf = rsf.replace(output_dir+'sRS-', '').replace('.pkl', '')
        EC_, DB, DB_ID = _rsf.split('-')
        EC = EC_.replace("_", ".")
        rs = reaction(EC, DB, DB_ID)
        if os.path.isfile(output_dir+rs.pkl) is False:
            print('you have not collected all reaction systems.')
            print('Exitting.')
            import sys
            sys.exit()
        # load in rxn system
        rs = rs.load_object(output_dir+rs.pkl, verbose=False)
        yield rs
