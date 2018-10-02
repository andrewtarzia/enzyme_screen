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
import pandas as pd


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


def load_molecule(filename, verbose=True):
    """unPickle molecule object from file.

    """
    if verbose:
        print('loading:', filename)
    with open(filename, 'rb') as input:
        mol = pickle.load(input)
        return mol


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
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'KEGG':
        from KEGG_IO import get_rxn_systems
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'BKMS':
        from BKMS_IO import get_rxn_systems
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'BRENDA':
        from BRENDA_IO import get_rxn_systems
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'ATLAS':
        from ATLAS_IO import get_rxn_systems
        # if 3rd tier EC only - do
        if EC.split('.')[3] != '-':
            return None
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


def collect_all_molecule_properties(output_dir, check=True):
    """Collect all molecule properties if they hadn't been collected during
    reaction system collection.

    """
    react_syst_files = glob.glob(output_dir+'sRS-*.pkl')
    # iterate over reaction system files
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            continue
        print('checking rxn', count, 'of', len(react_syst_files))
        for m in rs.components:
            m.get_properties(check)

        rs.save_object(output_dir+rs.pkl)


def process_molecule_collection(rs, output_dir, mol_db_dir, count,
                                react_syst_files):
    """Process the collection of molecule properties.

    """
    molecules = glob.glob(mol_db_dir+'ATRS_*.pkl')
    collect_RS_molecule_properties_parallel(rs=rs,
                                            output_dir=output_dir,
                                            mol_db_dir=mol_db_dir,
                                            count=count,
                                            molecules=molecules,
                                            react_syst_files=react_syst_files)


def collect_RS_molecule_properties_parallel(rs, output_dir, mol_db_dir,
                                            molecules,
                                            count=0, react_syst_files=[]):
    """Collect molecule properties from my database.

    """
    if rs.skip_rxn is True:
        return None
    if rs.components is None:
        rs.skip_rxn = True
        rs.save_object(output_dir+rs.pkl)
        return None
    print('checking rxn', count, 'of', len(react_syst_files))

    count_found = 0
    # collect properties from molecule DB
    for db_mol_pkl in molecules:
        if count_found == len(rs.components):
            break
        db_mol = load_molecule(db_mol_pkl, verbose=False)
        for m in rs.components:
            if db_mol.SMILES == m.SMILES:
                # copy DB object properties to RS
                # only overwrite None or NaN
                for key, val in db_mol.__dict__.items():
                    if key not in m.__dict__:
                        m.__dict__[key] = val
                    elif m.__dict__[key] is None and val is not None:
                        m.__dict__[key] = val
                rs.save_object(output_dir+rs.pkl)
                count_found += 1

    if count_found < len(rs.components):
        # if no match found for at least one molecule
        print('molecule not in database!')
        print('run molecule.py!')
        print('exitting...')
        sys.exit()


def collect_all_molecule_properties_parallel(rs, output_dir, check=True,
                                             count=0, react_syst_files=[]):
    """Collect all molecule properties if they hadn't been collected during
    reaction system collection.

    """
    if rs.skip_rxn is True:
        return None
    print('checking rxn', count, 'of', len(react_syst_files))
    for m in rs.components:
        m.get_properties(check)
    rs.save_object(output_dir+rs.pkl)


def process_RS_diffusion(rs, count, react_syst_files, output_dir,
                         threshold):
    """Process the check for diffusion of reaction system components.

    """
    check_RS_diffusion_parallel(rs=rs, count=count,
                                react_syst_files=react_syst_files,
                                output_dir=output_dir,
                                threshold=threshold)


def check_RS_diffusion_parallel(rs, count, react_syst_files,
                                output_dir, threshold):
    """Check all reaction systems for their component diffusion.

    All necessary properties should already be calculated!

    Keywords:
        rs (reaction object) - reaction system object
        count (int) - count of reactions tested
        react_syst_files (list) - list of RS files
        output_dir (str) - directory to output molecule files
        threshold (float) - diffusion size threshold

    """
    if rs.skip_rxn is True:
        return None
    # check if all_fit has already been done
    if rs.all_fit is not None:
        return None
    print('checking rxn', count, 'of', len(react_syst_files))
    # ignore any reactions with unknown components
    rs.skip_rxn = False
    for m in rs.components:
        if m.mol is None:
            rs.skip_rxn = True
    if rs.skip_rxn is True:
        print('skipping reaction - it is incomplete or generic')
        rs.save_object(output_dir+rs.pkl)
        return None

    all_fit = True
    max_comp_size = 0
    for m in rs.components:
        # remove reactions with general atoms (given by '*' in SMILES)
        if "*" in m.SMILES:
            rs.skip_rxn = True
            print('skipping reaction - it is incomplete or generic')
            rs.save_object(output_dir+rs.pkl)
            return None
        if m.mid_diam is None:
            print('molecule diameter has not been calculated yet')
            print('run molecule.py!')
            print('exitting...')
            sys.exit()
        # ignore any reactions with components with no sizes
        if m.mid_diam == 0:
            rs.skip_rxn = True
            print('molecule diameter could not be calculated for some reason.')
            print('skipping...')
            rs.save_object(output_dir+rs.pkl)
            return None
        max_comp_size = max([m.mid_diam, max_comp_size])
        if m.mid_diam > threshold or m.mid_diam == 0:
            all_fit = False

    if all_fit is True:
        rs.all_fit = True
        rs.max_comp_size = max_comp_size
        rs.print_rxn_system()
    else:
        rs.all_fit = False
        rs.max_comp_size = max_comp_size
    rs.save_object(output_dir+rs.pkl)


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
    molecule_output = DB_functions.initialize_mol_output_DF(
                        mol_output_file, overwrite=False)
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
                # check if IUPAC exists in molecule output
                if m.SMILES in list(molecule_output['SMILE']):
                    res_line = molecule_output[
                            molecule_output['SMILE'] == m.SMILES]
                    if str(res_line['iupac_name'].iloc[0]) != 'nan':
                        m.iupac_name = res_line['iupac_name'].iloc[0]
                # else use CIRPY
                else:
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


def check_all_solubility(output_dir):
    """Check all reaction systems for the solubility of their components.

    Defining solubility in terms of logP defined by Crippen et.al.

    Keywords:
        output_dir (str) - directory to output molecule files

    """

    # iterate over reaction system files
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            rs.min_logP = None
            rs.max_logP = None
            continue

        rs.min_logP = 100
        rs.max_logP = -100
        for m in rs.components:
            if m.logP is not None:
                rs.min_logP = min([rs.min_logP, m.logP])
                rs.max_logP = max([rs.max_logP, m.logP])

        if rs.min_logP == 100:
            rs.min_logP = None
        if rs.max_logP == -100:
            rs.max_logP = None

        rs.save_object(output_dir+rs.pkl)


def delta_sa_score(output_dir):
    """Get the change in maximum synthetic accessibility (sa) for all reactions.

    Keywords:
        output_dir (str) - directory to output molecule files

    """

    # iterate over reaction system files
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            rs.delta_sa = None
            rs.r_max_sa = None
            rs.p_max_sa = None
            continue

        rs.r_max_sa = -100
        rs.p_max_sa = -100
        for m in rs.components:
            if m.Synth_score is not None:
                if m.role == 'reactant':
                    rs.r_max_sa = max([rs.r_max_sa, m.Synth_score])
                elif m.role == 'product':
                    rs.p_max_sa = max([rs.p_max_sa, m.Synth_score])

        if rs.r_max_sa == -100:
            rs.r_max_sa = None
            rs.delta_sa = None
        elif rs.p_max_sa == -100:
            rs.p_max_sa = None
            rs.delta_sa = None
        else:
            rs.delta_sa = rs.p_max_sa - rs.r_max_sa

        rs.save_object(output_dir+rs.pkl)


def check_all_solubility_X(output_dir):
    """Check all reaction systems for the solubility of their components.

    Defining solubility in terms of XlogP from PUBCHEM

    Keywords:
        output_dir (str) - directory to output molecule files

    """

    # iterate over reaction system files
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            rs.min_XlogP = None
            rs.max_XlogP = None
            continue
        rs.min_XlogP = 100
        rs.max_XlogP = -100
        for m in rs.components:
            if m.XlogP is not None:
                rs.min_XlogP = min([rs.min_XlogP, float(m.XlogP)])
                rs.max_XlogP = max([rs.max_XlogP, float(m.XlogP)])

        if rs.min_XlogP == 100:
            rs.min_XlogP = None
        if rs.max_XlogP == -100:
            rs.max_XlogP = None

        rs.save_object(output_dir+rs.pkl)


def delta_complexity_score(output_dir):
    """Get the change in maximum synthetic accessibility (sa) for all reactions.

    Keywords:
        output_dir (str) - directory to output molecule files

    """

    # iterate over reaction system files
    count = 0
    for rs in yield_rxn_syst(output_dir):
        count += 1
        if rs.skip_rxn is True:
            rs.delta_comp = None
            rs.r_max_comp = None
            rs.p_max_comp = None
            continue

        rs.r_max_comp = -100000
        rs.p_max_comp = -100000
        for m in rs.components:
            if m.complexity is not None:
                if m.role == 'reactant':
                    rs.r_max_comp = max([rs.r_max_comp, float(m.complexity)])
                elif m.role == 'product':
                    rs.p_max_comp = max([rs.p_max_comp, float(m.complexity)])

        if rs.r_max_comp == -100000:
            rs.r_max_comp = None
            rs.delta_comp = None
        elif rs.p_max_comp == -100000:
            rs.p_max_comp = None
            rs.delta_comp = None
        else:
            rs.delta_comp = rs.p_max_comp - rs.r_max_comp

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


def get_ECs_from_file(EC_file):
    """Read in ECs to search from a file.

    """
    # get search EC numbers from file:
    # set EC numbers of interest
    # get from a data file - manually made from
    # https://enzyme.expasy.org/enzyme-byclass.html
    EC_DF = pd.read_table(EC_file, delimiter='__',
                          names=['EC_no', 'description'], engine='python')
    search_ECs = list(EC_DF['EC_no'])

    print(len(search_ECs), 'EC numbers to test')
    print('first EC:', search_ECs[0], '---- last EC:', search_ECs[-1])
    print('collect all reaction systems (ONLINE)...')
    return search_ECs


def wipe_reaction_properties(rs, output_dir):
    """Set attributes of rxn system to None.

    """
    rs.skip_rxn = False
    rs.all_fit = None  # do all the components fit?
    rs.max_comp_size = None
    rs.save_object(output_dir+rs.pkl)


if __name__ == "__main__":
    import sys
    import time
    from multiprocessing import Pool
    from molecule import molecule

    if (not len(sys.argv) == 5):
        print('Usage: rxn_syst.py run redo properties wipe\n')
        print('   run: T to run search for new rxn systems into current dir.')
        print('   redo: T to overwrite all rxn systems.')
        print('   properties: T to get properties of reaction systems in cwd.')
        print('   wipe: T to wipe properties of reaction systems in cwd.')
        sys.exit()
    else:
        run = sys.argv[1]
        redo = sys.argv[2]
        properties = sys.argv[3]
        wipe = sys.argv[4]

    if run == 'T':
        if redo == 'T':
            redo = True
        else:
            redo = False
        print('--------------------------------------------------------------')
        print('Screen new reactions')
        print('--------------------------------------------------------------')
        temp_time = time.time()
        search_DBs = ['BRENDA', 'SABIO', 'KEGG', 'BKMS', ]
        NP = 1  # number of processes
        search_EC_file = 'desired_EC.txt'
        print('settings:')
        print('    EC file:', search_EC_file)
        print('    Number of processes:', NP)
        print('    DBs to search:', search_DBs)
        inp = input('happy with these? (T/F)')
        if inp == 'F':
            sys.exit('change them in the source code')
        elif inp != 'T':
            sys.exit('I dont understand, T or F?')
        print('collect all reaction systems (ONLINE)...')
        search_ECs = get_ECs_from_file(EC_file=search_EC_file)
        search_output_dir = os.getcwd()+'/'
        print(search_output_dir)
        for DB in search_DBs:
            # iterate over EC numbers of interest
            # Create a multiprocessing Pool
            with Pool(NP) as pool:
                # process data_inputs iterable with pool
                # func(EC, DB, search_output_dir, search_redo, verbose)
                args = [(EC, DB, search_output_dir, redo, True)
                        for EC in search_ECs]
                pool.starmap(process_collection, args)
        percent_skipped(search_output_dir)
        print('---- time taken =', '{0:.2f}'.format(time.time()-temp_time),
              's')
    if wipe == 'T':
        print('--------------------------------------------------------------')
        print('Wipe reaction properties')
        print('--------------------------------------------------------------')
        search_output_dir = os.getcwd()+'/'
        react_syst_files = glob.glob(search_output_dir+'sRS-*.pkl')
        count = 0
        for rs in yield_rxn_syst(search_output_dir):
            print('wiping', count, 'of', len(react_syst_files))
            count += 1
            wipe_reaction_properties(rs, search_output_dir)
    if properties == 'T':
        print('--------------------------------------------------------------')
        print('Collect reaction properties')
        print('--------------------------------------------------------------')
        NP = 2  # number of processes
        pI_thresh = 6
        size_thresh = 4.2  # angstroms
        molecule_db_dir = '/home/atarzia/psp/molecule_DBs/atarzia/'
        print('settings:')
        print('    Number of processes:', NP)
        print('    pI threshold:', pI_thresh)
        print('    Diffusion threshold:', size_thresh, 'Angstrom')
        print('    Molecule database:', molecule_db_dir)
        inp = input('happy with these? (T/F)')
        if inp == 'F':
            sys.exit('change them in the source code')
        elif inp != 'T':
            sys.exit('I dont understand, T or F?')
        search_output_dir = os.getcwd()+'/'
        react_syst_files = glob.glob(search_output_dir+'sRS-*.pkl')
        print('collect component properties from molecule database...')
        # these should be calculated already using molecule.py
        # iterate over reaction systems
        # Create a multiprocessing Pool
        with Pool(NP) as pool:
            # process data_inputs iterable with pool
            # func(rs, output_dir, mol_db_dir, count, react_syst_files)
            args = [(rs,
                     search_output_dir, molecule_db_dir, i, react_syst_files)
                    for i, rs in enumerate(yield_rxn_syst(search_output_dir))]
            pool.starmap(process_molecule_collection, args)
        print('check all reaction systems for diffusion of components...')
        # iterate over reaction systems
        # Create a multiprocessing Pool
        with Pool(NP) as pool:
            # process data_inputs iterable with pool
            # func(rs, count, react_syst_files, search_output_dir, size_thrsh)
            args = [(rs, i, react_syst_files, search_output_dir,
                     size_thresh)
                    for i, rs in enumerate(yield_rxn_syst(search_output_dir))]
            pool.starmap(process_RS_diffusion, args)

        print('get subset of reactions with known protein sequences...')
        check_all_seedMOF(search_output_dir, pI_thresh)
        percent_w_sequence(search_output_dir)
        print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
        temp_time = time.time()
        print('determine solubility range of all reactions using logP...')
        check_all_solubility(output_dir=search_output_dir)
        print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
        temp_time = time.time()
        print('determine change in synthetic accessibility all reactions...')
        delta_sa_score(output_dir=search_output_dir)
        print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
        temp_time = time.time()
        print('determine solubility range of all reactions using XlogP...')
        check_all_solubility_X(output_dir=search_output_dir)
        print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
        temp_time = time.time()
        print('determine change in molecular complexity all reactions...')
        delta_complexity_score(output_dir=search_output_dir)
        print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')
