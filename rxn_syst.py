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
import DB_functions
import Uniprot_IO
import pi_fn


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
        if self.seed_MOF is True:
            if self.req_mod is None:
                print('Has sequence that seeds MOF growth w/o modification!')
            elif self.req_mod is '1':
                print('Has sequence that seeds MOF growth w succinylation!')
            print('-----------------------------------')


def get_reaction_systems(EC, DB, output_dir, clean_system=False,
                         verbose=False):
    """Get reaction system from SABIO reaction ID (rID).

    Keywords:
        EC (str) - Enzyme commision number (X.X.X.X)
        DB (str) - name of Database
        output_dir (str) - directory where all data should be saved
        clean_system (bool) - wipe the data in reaction systems for fresh start
            default = False
        verbose (bool) - print update
            default = False

    """
    if DB == 'SABIO':
        from SABIO_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'KEGG':
        from KEGG_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'BKMS':
        from BKMS_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'BRENDA':
        from BRENDA_IO import get_rxn_systems
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)


def process_collection(EC, DB, search_output_dir,
                       search_redo, verbose):
    """Process the collection of new reaction systems.

    """
    # iterate over EC numbers of interest
    print('doing:', DB, 'EC:', EC)
    get_reaction_systems(EC, DB,
                         search_output_dir,
                         clean_system=search_redo,
                         verbose=verbose)


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


def percent_skipped(output_dir):
    """Print the percent of all reaction systems that will NOT be skipped.

    """
    # what percentage of reaction systems have skip_rxn = False
    count = 0
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is False:
            count += 1

    print('-----------------------------------')
    print(count, 'reaction systems of', len(react_syst_files),
          'are NOT skipped.')
    print('=>', round(count/len(react_syst_files), 4)*100, 'percent')
    print('-----------------------------------')


def percent_w_sequence(output_dir):
    """Print the percent of all reaction systems with a protein sequence.

    """
    # what percentage of reaction systems have skip_rxn = False
    count = 0
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
    for rs in yield_rxn_syst(output_dir):
        if rs.UniprotID != '':
            if rs.UniprotID is not None:
                count += 1

    print('-----------------------------------')
    print(count, 'reaction systems of', len(react_syst_files),
          'had a sequence.')
    print('=>', round(count/len(react_syst_files), 4)*100, 'percent')
    print('-----------------------------------')


def check_all_RS_diffusion(output_dir, mol_output_file, threshold,
                           vdwScale, boxMargin, spacing, N_conformers,
                           MW_thresh):
    """Check all reaction systems for their component diffusion.

    Keywords:
        output_dir (str) - directory to output molecule files
        mol_output_file (str) - molecule_output file
        threshold (float) - diffusion size threshold
        vdwScale (float) - Scaling factor for the radius of the atoms to
            determine the base radius used in the encoding
            - grid points inside this sphere carry the maximum occupancy
        boxMargin (float) - added margin to grid surrounding molecule
        spacing (float) - grid spacing
        N_conformers (int) - number of conformers to calculate diameter of
        MW_thresh (float) - Molecular Weight maximum

    """
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
    # iterate over reaction system files
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            continue
        # check if all_fit has already been done
        if rs.all_fit is not None:
            continue
        print('checking rxn', count, 'of', len(react_syst_files))
        # define reactants and products dict
        # name: (smile, DB, DB_ID, iupac_name, role)
        components_dict = {}

        # ignore any reactions with unknown components
        rs.skip_rxn = False
        for m in rs.components:
            if m.mol is None:
                rs.skip_rxn = True

        if rs.skip_rxn is True:
            print('skipping reaction - it is incomplete or generic')
            rs.save_object(output_dir+rs.pkl)
            continue
        for m in rs.components:
            # get IUPAC NAME
            if m.iupac_name is None:
                m.cirpy_to_iupac()
            # remove reactions with general atoms (given by '*' in SMILES)
            if "*" in m.SMILES:
                rs.skip_rxn = True
                print('skipping reaction - it is incomplete or generic')
                rs.save_object(output_dir+rs.pkl)
                break
            components_dict[m.name] = (m.SMILES, m.DB, m.DB_ID,
                                       m.iupac_name, m.role.lower())

        # calculate molecule size of all components
        # save to molecule output files
        molecule_output = DB_functions.initialize_mol_output_DF(
                            mol_output_file, overwrite=False)
        DB_functions.get_molecule_diameters(components_dict,
                                            molecule_output=molecule_output,
                                            mol_output_file=mol_output_file,
                                            out_dir=output_dir,
                                            vdwScale=vdwScale,
                                            boxMargin=boxMargin,
                                            spacing=spacing,
                                            N_conformers=N_conformers,
                                            MW_thresh=MW_thresh)

        molecule_output = DB_functions.initialize_mol_output_DF(
                            mol_output_file, overwrite=False)

        # get diameters (should alrady be calculated)
        # of all components of reaction
        rs.check_all_fit(threshold, molecule_output)

        # ignore any reactions with components with no sizes
        for m in rs.components:
            if m.mid_diam is None or m.mid_diam == 0:
                rs.skip_rxn = True
        rs.save_object(output_dir+rs.pkl)


def check_all_seedMOF(output_dir, pI_thresh):
    """Check all reaction systems an associated protein sequence and whether it
    seeds MOF growth.

    Keywords:
        output_dir (str) - directory to output molecule files
        pI_thresh (float) - thrshold pI for MOF growth

    """
    # iterate over reaction system files
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            continue
        # this is only possible for reaction systems with UNIPROT ID
        if rs.UniprotID is None or rs.UniprotID == '':
            continue
        # pI already checked?
        if rs.seed_MOF is not None:
            continue
        print('checking rxn', count, 'of', len(react_syst_files))
        # split UniprotID for the cases where multiple subunits exist
        IDs = rs.UniprotID.split(" ")
        print('Uniprot IDs:', IDs)
        if len(IDs) > 0:
            # iterate over all UniProtIDs
            # assume all sequences require pI < cutoff for MOF growth
            # this is done by collating all sequences
            total_sequence = ''
            for i in IDs:
                sequence = Uniprot_IO.get_sequence(i)
                total_sequence += sequence
            if len(total_sequence) > 0:
                rs = pi_fn.calculate_rxn_syst_pI(total_sequence, rs,
                                                 cutoff_pi=pI_thresh)
                print('seed MOF?', rs.seed_MOF)
            else:
                rs.pI = None
                rs.seed_MOF = 'unclear'
            rs.save_object(output_dir+rs.pkl)
        else:
            rs.save_object(output_dir+rs.pkl)
            print('seed MOF?', rs.seed_MOF)
