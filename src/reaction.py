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
        print('no components had properties.')
        self.skip_rxn = True
        self.skip_reason = 'No components have properties'

    def fail_size(self):
        print(
            'a component is too big by MW or size could not be calc.'
        )
        self.skip_rxn = True
        self.skip_reason = 'A components size is too big, or failed.'

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
        self.get_purchasability_class()

    def get_max_mid_diameter(self):
        """
        Get diffusion properties of all components.

        """
        max_min_mid_diam = 0

        for m in self.components:
            name = m.name
            diam_file = join(
                self.params['molec_dir'],
                name+'_size.csv'
            )

            if exists(diam_file.replace('.csv', '.TOOBIG')):
                max_min_mid_diam = 0
                print(f'{m.name} too big based on MW')
                break
            if exists(diam_file.replace(
                'size.csv',
                'unopt.ETKDGFAILED'
            )):
                max_min_mid_diam = 0
                print(f'{m.name} failed ETKDG')
                break
            results = pd.read_csv(diam_file)
            min_mid_diam = min(results['diam2'])
            max_min_mid_diam = max([min_mid_diam, max_min_mid_diam])

        self.max_min_mid_diam = max_min_mid_diam

    def get_logPs(self):
        """
        Get logP (hydrophobicity) properties of all components.

        """

        self.min_logP = 1E10
        self.max_logP = -1E10
        max_NHA = 0
        for m in self.components:
            # Ignore water.
            if m.name in KEGG_IDs_to_ignore():
                continue
            prop_dict = m.read_prop_file()

            # Only get properties of the largest component.
            # Number of heavy atoms.
            NHA = prop_dict['NHA']
            if NHA >= max_NHA:
                self.min_logP = min([self.min_logP, prop_dict['logP']])
                self.max_logP = max([self.max_logP, prop_dict['logP']])
                max_NHA = NHA

    def get_logSs(self):
        """
        Get logS (water solubility) properties of all components.

        """

        self.min_logS = 1E10
        self.max_logS = -1E10
        max_NHA = 0
        for m in self.components:
            # Ignore water.
            if m.name in KEGG_IDs_to_ignore():
                continue
            prop_dict = m.read_prop_file()

            # Only get properties of the largest component.
            # Number of heavy atoms.
            NHA = prop_dict['NHA']
            if NHA >= max_NHA:
                self.min_logS = min([self.min_logS, prop_dict['logS']])
                self.max_logS = max([self.max_logS, prop_dict['logS']])
                max_NHA = NHA

    def get_purchasability_class(self):
        """
        Determine purchasability class of reaction.

        Class I
            All components purchasable.

        Class II
            A substrate not purchasable, A product is purchasable.

        Class III
            A substrate purchasable, A product not purchasable.

        Class IIII
            A substrate not purchasable, A product not purchasable.

        """
        self.p_class = 1
        r_purch = []
        p_purch = []
        for m in self.components:
            # Ignore water.
            if m.name in KEGG_IDs_to_ignore():
                continue
            prop_dict = m.read_prop_file()
            purchasability = prop_dict['purchasability']
            # print(m)
            # print(purchasability)
            # Only get properties of the largest component.
            # Number of heavy atoms.
            if m.role == 'reactant':
                r_purch.append(purchasability)
            elif m.role == 'product':
                p_purch.append(purchasability)

        # print(self.DB_ID)
        # print('rs', r_purch)
        # print(all(r_purch), not all(r_purch))
        # print('ps', p_purch)
        # print(all(p_purch), not all(p_purch))
        if all(r_purch) and all(p_purch):
            self.p_class = 1
        elif not all(r_purch) and all(p_purch):
            self.p_class = 2
        elif all(r_purch) and not all(p_purch):
            self.p_class = 3
        else:
            self.p_class = 4

        # print(self.p_class)
        # input('check the purchasability on zinc')

    def get_SAs(self):
        """
        Get synthetic accessibility properties of all components.

        For each reaction, only report the SA for the largest
        components of the reactants and products.

        """
        self.r_max_SA = 0
        self.p_max_SA = 0
        self.delta_SA = 0
        r_max_NHA = 0
        p_max_NHA = 0
        for m in self.components:
            # Ignore water.
            if m.name in KEGG_IDs_to_ignore():
                continue
            prop_dict = m.read_prop_file()

            # Number of heavy atoms.
            NHA = prop_dict['NHA']
            if m.role == 'reactant' and NHA >= r_max_NHA:
                self.r_max_SA = max([
                    self.r_max_SA,
                    prop_dict['Synth_score']
                ])
                r_max_NHA = NHA
            elif m.role == 'product' and NHA >= p_max_NHA:
                self.p_max_SA = max([
                    self.p_max_SA,
                    prop_dict['Synth_score']
                ])
                p_max_NHA = NHA

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


def KEGG_IDs_to_ignore():
    """
    A list of KEGG IDs to ignore in properties.

    Most of these are common solvents OR are easily derived from
    metal salts or similar compounds.

    Most of these come from those that end up being not purchasable,
    but are accessible to me (a computational chemist).

    Molecules to ignore:
        water, oxygen, carbon dioxide, H2O2, Mn2+, H+,
        Mn3+, cyanide ion, Hg2+, Hg, vinyl chloride, HF, HBr, methane,
        Fe2+, NO, HI, formate, nitroxyl, Na+, dimethylamine, II,
        I-, Co2+, Br-, DCM, ethene, hypochlorous acid, ethyne,
        nitrous oxide, ethanol, H2S, Cl-, hydrazine, nitrate,
        hydroxide ion, formamide, hydrogen cyanide, Fe+3,
        elemental sulfur, ethanal, chlorus acid, sulfur dioxide,
        methanol, superoxide anion, carbon monoxide, nitric oxide,
        F-, ammonia, ethyl chloride, Ni+2, hydrogen selenide,
        sulfurous acid, hyponitrous acid, cyanic acid, arsenite,
        methylamine, ethylamine, chlorate, glycol, nitrite, nitrogen,
        Mg2+, hydrogen, carboynl sulfide, HCl, azide, hydrogen sulfite,
        selenium, diimine, trimethylamine, formaldehyde,

    """
    ids = [
        'C00001', 'C00007', 'C00011', 'C00027', 'C19610', 'C00080',
        'C19610', 'C00177', 'C19611', 'C00703', 'C01319', 'C06793',
        'C16487', 'C13645', 'C01438', 'C14818', 'C00192', 'C05590',
        'C00058', 'C20937', 'C01330', 'C00543', 'C01382', 'C00708',
        'C00175', 'C01324', 'C02271', 'C06547', 'C19697', 'C01548',
        'C00887', 'C00469', 'C00283', 'C00698', 'C05361', 'C00244',
        'C01328', 'C00488', 'C01326', 'C14819', 'C00087', 'C00084',
        'C01486', 'C09306', 'C00132', 'C00704', 'C00237', 'C00533',
        'C00742', 'C00014', 'C18248', 'C19609', 'C01528', 'C00094',
        'C01818', 'C01417', 'C06697', 'C00218', 'C00797', 'C01485',
        'C15588', 'C00088', 'C00697', 'C00305', 'C00282', 'C07331',
        'C01327', 'C19935', 'C11481', 'C01529', 'C05360', 'C00565',
        'C00067',
    ]

    return ids
