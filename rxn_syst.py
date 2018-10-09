#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the reaction system class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import pickle
import bz2
import glob
import os
import pandas as pd


class reaction:
    """Class that defines a reaction system for all databases.

    """

    def __init__(self, EC, DB, DB_ID):
        # all non DB unique properties
        self.EC = EC
        self.DB = DB
        self.DB_ID = DB_ID
        # for unknwon EC tiers (given by '-'), use a known delimeter.
        EC_ul = EC.replace('.', '_').replace('-', 'XX')
        self.pkl = 'sRS-'+EC_ul+'-'+str(DB)+'-'+str(DB_ID)+'.bpkl'
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
        with bz2.BZ2File(filename.replace('.pkl', '.bpkl'), 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(self, filename, verbose=True):
        """unPickle reaction system object from file.

        """
        if verbose:
            print('loading:', filename.replace('.pkl', '.bpkl'))
        with bz2.BZ2File(filename.replace('.pkl', '.bpkl'), 'rb') as input:
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
    with bz2.BZ2File(filename.replace('.pkl', '.bpkl'), 'rb') as input:
        mol = pickle.load(input)
        return mol


def get_reaction_systems(EC, DB, output_dir, molecule_dataset,
                         clean_system=False,
                         verbose=False):
    """Get reaction system from SABIO reaction ID (rID).

    Keywords:
        EC (str) - Enzyme commision number (X.X.X.X)
        DB (str) - name of Database
        output_dir (str) - directory where all data should be saved
        molecule_dataset (Pandas DataFrame) - look up for known molecules
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
        get_rxn_systems(EC, output_dir, molecule_dataset=molecule_dataset,
                        clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'KEGG':
        from KEGG_IO import get_rxn_systems
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, molecule_dataset=molecule_dataset,
                        clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'BKMS':
        from BKMS_IO import get_rxn_systems
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, molecule_dataset=molecule_dataset,
                        clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'BRENDA':
        from BRENDA_IO import get_rxn_systems
        # if 3rd tier EC only - skip
        if EC.split('.')[3] == '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, molecule_dataset=molecule_dataset,
                        clean_system=clean_system,
                        verbose=verbose)
    elif DB == 'ATLAS':
        from ATLAS_IO import get_rxn_systems
        # if 3rd tier EC only - do
        if EC.split('.')[3] != '-':
            return None
        # set DB specific properties
        get_rxn_systems(EC, output_dir, molecule_dataset=molecule_dataset,
                        clean_system=clean_system,
                        verbose=verbose)


def get_RS(filename, output_dir, verbose=False):
    """Read in reaction system from filename.

    """
    _rsf = filename.replace(output_dir+'sRS-', '').replace('.bpkl', '')
    EC_, DB, DB_ID = _rsf.split('-')
    EC = EC_.replace("_", ".").replace('XX', '-')
    rs = reaction(EC, DB, DB_ID)
    if os.path.isfile(output_dir+rs.pkl) is False:
        print('you have not collected all reaction systems.')
        print('Exitting.')
        import sys
        sys.exit()
    # load in rxn system
    if verbose:
        print('loading:', rs.pkl)
    rs = rs.load_object(output_dir+rs.pkl, verbose=False)
    return rs


def yield_rxn_syst(output_dir, verbose=False):
    """Iterate over reaction systems for analysis.

    """
    react_syst_files = glob.glob(output_dir+'sRS-*.bpkl')
    for rsf in react_syst_files:
        rs = get_RS(filename=rsf, output_dir=output_dir, verbose=verbose)
        yield rs


def percent_skipped(output_dir):
    """Print the percent of all reaction systems that will NOT be skipped.

    """
    # what percentage of reaction systems have skip_rxn = False
    count = 0
    count_atlas = 0
    count_brenda = 0
    count_bkms = 0
    count_kegg = 0
    count_sabio = 0
    react_syst_files = glob.glob(output_dir+'sRS-*.bpkl')
    rsf_atlas = glob.glob(output_dir+'sRS-*ATLAS*.bpkl')
    rsf_brenda = glob.glob(output_dir+'sRS-*BRENDA*.bpkl')
    rsf_bkms = glob.glob(output_dir+'sRS-*BKMS*.bpkl')
    rsf_kegg = glob.glob(output_dir+'sRS-*KEGG*.bpkl')
    rsf_sabio = glob.glob(output_dir+'sRS-*SABIO*.bpkl')
    for rs in yield_rxn_syst(output_dir):
        if rs.skip_rxn is False:
            count += 1
            if 'ATLAS' in rs.pkl:
                count_atlas += 1
            if 'BRENDA' in rs.pkl:
                count_brenda += 1
            if 'BKMS' in rs.pkl:
                count_bkms += 1
            if 'KEGG' in rs.pkl:
                count_kegg += 1
            if 'SABIO' in rs.pkl:
                count_sabio += 1

    print('-----------------------------------')
    print(count, 'reaction systems of', len(react_syst_files),
          'are NOT skipped.')
    print('=>', round(count/len(react_syst_files), 4)*100, 'percent')
    print('-----------------------------------')
    print(count_atlas, 'reaction systems of', len(rsf_atlas),
          'are NOT skipped in the ATLAS data set.')
    print('-----------------------------------')
    print(count_brenda, 'reaction systems of', len(rsf_brenda),
          'are NOT skipped in the BRENDA data set.')
    print('-----------------------------------')
    print(count_bkms, 'reaction systems of', len(rsf_bkms),
          'are NOT skipped in the BKMS data set.')
    print('-----------------------------------')
    print(count_kegg, 'reaction systems of', len(rsf_kegg),
          'are NOT skipped in the KEGG data set.')
    print('-----------------------------------')
    print(count_sabio, 'reaction systems of', len(rsf_sabio),
          'are NOT skipped in the SABIO data set.')
    print('-----------------------------------')


def collect_RS_molecule_properties(rs, output_dir, mol_db_dir,
                                   molecules,
                                   count=0, react_syst_files=[]):
    """Collect molecule properties from my database.

    """
    try:
        if rs.mol_collected is True:
            return None
    except AttributeError:
        rs.mol_collected = False
        rs.save_object(output_dir+rs.pkl)
    if rs.skip_rxn is True:
        return None
    if rs.components is None:
        rs.skip_rxn = True
        rs.save_object(output_dir+rs.pkl)
        return None
    # ignore any reactions with unknown components
    rs.skip_rxn = False
    for m in rs.components:
        print('component:', m.name)
        try:
            if m.mol is None:
                rs.skip_rxn = True
        except AttributeError:
            rs.skip_rxn = True
    if rs.skip_rxn is True:
        print('skipping reaction - it is incomplete or generic')
        rs.save_object(output_dir+rs.pkl)
        return None
    print('checking rxn', count, 'of', len(react_syst_files))
    print('-----')
    print(rs.pkl)
    count_found = 0
    # collect properties from molecule DB
    for db_mol_pkl in molecules:
        if count_found == len(rs.components):
            break
        db_mol = load_molecule(db_mol_pkl, verbose=False)
        for m in rs.components:
            if m.SMILES is None:
                continue
            if db_mol.SMILES == m.SMILES:
                print(db_mol_pkl)
                # copy DB object properties to RS
                # only overwrite None or NaN
                for key, val in db_mol.__dict__.items():
                    if key not in m.__dict__:
                        m.__dict__[key] = val
                    elif m.__dict__[key] is None and val is not None:
                        m.__dict__[key] = val
                m.pkl = db_mol_pkl
                rs.save_object(output_dir+rs.pkl)
                print('found -', m.name)
                count_found += 1

    if count_found < len(rs.components):
        # if no match found for at least one molecule
        print('molecule not in database!')
        print('run molecule.py!')
        print('exitting...')
        sys.exit()

    rs.mol_collected = True
    rs.save_object(output_dir+rs.pkl)


def RS_diffusion(rs, output_dir, threshold):
    """Check all reaction systems for their component diffusion.

    All necessary properties should already be calculated!

    Keywords:
        rs (reaction object) - reaction system object
        output_dir (str) - directory to output molecule files
        threshold (float) - diffusion size threshold

    """
    if rs.skip_rxn is True:
        return None
    # check if all_fit has already been done
    if rs.all_fit is not None:
        return None
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
    print('----')
    for m in rs.components:
        print(m.name, ':', m.mid_diam)
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


def RS_solubility(rs, output_dir):
    """Check solubility properties of a reaction system.

    Defining solubility in terms of logP defined by Crippen et.al.

    Keywords:
        rs (reaction) - reaction system
        output_dir (str) - directory to output molecule files

    """

    if rs.skip_rxn is True:
        rs.min_logP = None
        rs.max_logP = None
        return None
    rs.min_logP = 100
    rs.max_logP = -100
    for m in rs.components:
        # ignore reactions with unknown values
        if m.logP is None or m.logP == 'not found':
            rs.min_logP = 100
            rs.max_logP = -100
            break
        rs.min_logP = min([rs.min_logP, m.logP])
        rs.max_logP = max([rs.max_logP, m.logP])
    if rs.min_logP == 100:
        rs.min_logP = None
    if rs.max_logP == -100:
        rs.max_logP = None
    rs.save_object(output_dir+rs.pkl)


def RS_SAscore(rs, output_dir):
    """Get the change in maximum synthetic accessibility (sa) a reaction system.

    Keywords:
        rs (reaction) - reaction system
        output_dir (str) - directory to output molecule files

    """
    if rs.skip_rxn is True:
        rs.delta_sa = None
        rs.r_max_sa = None
        rs.p_max_sa = None
        return None
    rs.r_max_sa = -100
    rs.p_max_sa = -100
    for m in rs.components:
        # ignore reactions with unknown values
        if m.Synth_score is None or m.Synth_score == 'not found':
            rs.r_max_sa = -100
            rs.p_max_sa = -100
            break
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


def RS_solubility_X(rs, output_dir):
    """Check solubility properties of a reaction system.

    Defining solubility in terms of XlogP from PUBCHEM

    Keywords:
        rs (reaction) - reaction system
        output_dir (str) - directory to output molecule files

    """
    if rs.skip_rxn is True:
        rs.min_XlogP = None
        rs.max_XlogP = None
        return None
    rs.min_XlogP = 100
    rs.max_XlogP = -100
    for m in rs.components:
        # ignore reactions with unknown values
        if m.XlogP is None or m.XlogP == 'not found':
            rs.min_XlogP = 100
            rs.max_XlogP = -100
            break
        try:
            rs.min_XlogP = min([rs.min_XlogP, float(m.XlogP)])
            rs.max_XlogP = max([rs.max_XlogP, float(m.XlogP)])
        except ValueError:
            print('temporary fix for PUBCHEM errors')
            print('fix this!')
            rs.min_XlogP = 100
            rs.max_XlogP = -100
            m.XlogP = None
    if rs.min_XlogP == 100:
        rs.min_XlogP = None
    if rs.max_XlogP == -100:
        rs.max_XlogP = None
    rs.save_object(output_dir+rs.pkl)


def RS_complexity_score(rs, output_dir):
    """Get the change in maximum complexity for a reaction system.

    Keywords:
        rs (reaction) - reaction system
        output_dir (str) - directory to output molecule files

    """
    if rs.skip_rxn is True:
        rs.delta_comp = None
        rs.r_max_comp = None
        rs.p_max_comp = None
        return None
    rs.r_max_comp = -100000
    rs.p_max_comp = -100000
    for m in rs.components:
        # ignore reactions with unknown values
        if m.complexity is None or m.complexity == 'not found':
            rs.r_max_comp = -100000
            rs.p_max_comp = -100000
            break
        try:
            if m.role == 'reactant':
                rs.r_max_comp = max([rs.r_max_comp, float(m.complexity)])
            elif m.role == 'product':
                rs.p_max_comp = max([rs.p_max_comp, float(m.complexity)])
        except ValueError:
            print('temporary fix for PUBCHEM errors')
            print('fix this!')
            rs.r_max_comp = -100000
            rs.p_max_comp = -100000
            m.complexity = None
    if rs.r_max_comp == -100000:
        rs.r_max_comp = None
        rs.delta_comp = None
    elif rs.p_max_comp == -100000:
        rs.p_max_comp = None
        rs.delta_comp = None
    else:
        rs.delta_comp = rs.p_max_comp - rs.r_max_comp
    rs.save_object(output_dir+rs.pkl)


def parallel_analysis(rs, count, react_syst_files, output_dir, threshold,
                      mol_db_dir, molecules, done_pkls):
    """Run analysis in parallel.

    Keywords:
        count (int) - count of reactions tested
        react_syst_files (list) - list of RS files

    """
    print('checking rxn', count, 'of', len(react_syst_files))
    collect_RS_molecule_properties(rs=rs, output_dir=search_output_dir,
                                   mol_db_dir=mol_db_dir,
                                   molecules=molecules, count=count,
                                   react_syst_files=react_syst_files)
    if rs.pkl not in done_pkls:
        RS_diffusion(rs=rs, output_dir=output_dir,
                     threshold=threshold)
        RS_solubility(rs=rs, output_dir=output_dir)
        RS_SAscore(rs=rs, output_dir=output_dir)
        RS_solubility_X(rs=rs, output_dir=output_dir)
        RS_complexity_score(rs=rs, output_dir=output_dir)
        with open(output_dir+'prop_done.txt', 'a') as f:
            f.write(rs.pkl+'\n')


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
    print('wiping: skip_rxn, all_fit, max_comp_size.')
    rs.skip_rxn = False
    rs.all_fit = None  # do all the components fit?
    rs.max_comp_size = None
    rs.save_object(output_dir+rs.pkl)


def main_run(redo):
    """Run reaction system collection.

    """
    if redo == 'T':
        redo = True
    else:
        redo = False
    print('--------------------------------------------------------------')
    print('Screen new reactions')
    print('--------------------------------------------------------------')
    temp_time = time.time()
    DB_switch = input('biomin (1) or new (2)?')
    if DB_switch == '1':
        search_DBs = ['SABIO', 'BRENDA', 'KEGG', 'BKMS', ]
    elif DB_switch == '2':
        search_DBs = ['SABIO', 'ATLAS', 'BRENDA', 'KEGG', 'BKMS', ]
    else:
        print('answer correctly...')
        sys.exit()
    NP = 1  # number of processes
    search_EC_file = 'desired_EC.txt'
    lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
    molecule_dataset = pd.read_table(lookup_file, delimiter='___',
                                     skiprows=[0],
                                     names=['SMILES', 'iupac', 'name',
                                            'DB', 'DB_ID',
                                            'KEGG_ID', 'pkl'])
    print('settings:')
    print('    EC file:', search_EC_file)
    print('    Number of processes:', NP)
    print('    DBs to search:', search_DBs)
    print('    Molecule DB lookup file:', lookup_file)
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
        if NP > 1:
            # Create a multiprocessing Pool
            with Pool(NP) as pool:
                # process data_inputs iterable with pool
                # func(EC, DB, search_output_dir, mol dataset, search_redo,
                #      verbose)
                args = [(EC, DB, search_output_dir, molecule_dataset, redo,
                         True)
                        for EC in search_ECs]
                pool.starmap(get_reaction_systems, args)
        # in serial
        else:
            for EC in search_ECs:
                get_reaction_systems(EC=EC, DB=DB,
                                     output_dir=search_output_dir,
                                     molecule_dataset=molecule_dataset,
                                     clean_system=redo, verbose=True)
    percent_skipped(search_output_dir)
    print('---- time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')


def main_wipe():
    """Wipe reaction system properties to rerun analysis

    """
    print('--------------------------------------------------------------')
    print('Wipe reaction properties')
    print('--------------------------------------------------------------')
    inp = input('are you sure? (T/F)')
    if inp == 'F':
        sys.exit('')
    elif inp != 'T':
        sys.exit('I dont understand, T or F?')
    search_output_dir = os.getcwd()+'/'
    react_syst_files = glob.glob(search_output_dir+'sRS-*.bpkl')
    count = 0
    for rs in yield_rxn_syst(search_output_dir):
        print('wiping', count, 'of', len(react_syst_files))
        count += 1
        wipe_reaction_properties(rs, search_output_dir)


def main_analysis(prop_redo):
    """Analyse all reaction systems.

    """
    print('--------------------------------------------------------------')
    print('Collect reaction properties')
    print('--------------------------------------------------------------')
    NP = 1  # number of processes
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
    react_syst_files = glob.glob(search_output_dir+'sRS-*.bpkl')
    molecules = glob.glob(molecule_db_dir+'ATRS_*.bpkl')
    print('---------------------------------------------------------------')
    print('collect component properties and analyse reaction systems:')
    print('    - diffusion of components')
    print('    - solubility (logP) of components')
    print('    - change in synthetic accessibility of components')
    print('    - solubiity (XlogP) of components')
    print('    - change in complexity of components')
    print('---------------------------------------------------------------')
    temp_time = time.time()
    if prop_redo is True:
        with open(search_output_dir+'prop_done.txt', 'w') as f:
            f.write('pkls\n')
        done_pkls = []
    else:
        done_pkls = []
        with open(search_output_dir+'prop_done.txt', 'r') as f:
            for line in f.readlines():
                done_pkls.append(line.rstrip())
    # iterate over reaction systems
    if NP > 1:
        # Create a multiprocessing Pool
        with Pool(NP) as pool:
            # process data_inputs iterable with pool
            # func(rs, count, react_syst_files, search_output_dir, size_thrsh)
            args = [(rs, i, react_syst_files, search_output_dir,
                     size_thresh, molecule_db_dir, molecules, done_pkls)
                    for i, rs in enumerate(yield_rxn_syst(search_output_dir))]
            pool.starmap(parallel_analysis, args)
    # in serial
    else:
        for i, rs in enumerate(yield_rxn_syst(search_output_dir, verbose=True)):
            print('checking rxn', i, 'of', len(react_syst_files))
            collect_RS_molecule_properties(rs=rs, output_dir=search_output_dir,
                                           mol_db_dir=molecule_db_dir,
                                           molecules=molecules, count=i,
                                           react_syst_files=react_syst_files)
            if rs.pkl not in done_pkls:
                RS_diffusion(rs=rs, output_dir=search_output_dir,
                             threshold=size_thresh)
                RS_solubility(rs=rs, output_dir=search_output_dir)
                RS_SAscore(rs=rs, output_dir=search_output_dir)
                RS_solubility_X(rs=rs, output_dir=search_output_dir)
                RS_complexity_score(rs=rs, output_dir=search_output_dir)
                with open(search_output_dir+'prop_done.txt', 'a') as f:
                    f.write(rs.pkl+'\n')

    print('--- time taken =', '{0:.2f}'.format(time.time()-temp_time), 's')


if __name__ == "__main__":
    import sys
    import time
    from multiprocessing import Pool
    from molecule import molecule

    if (not len(sys.argv) == 7):
        print('Usage: rxn_syst.py run redo properties wipe\n')
        print('   run: T to run search for new rxn systems into current dir.')
        print('   redo: T to overwrite all rxn systems.')
        print('   properties: T to get properties of reaction systems in cwd.')
        print('   rerun properites?: T for rerun, F to read from prop_done.txt.')
        print('   wipe: T to wipe properties of reaction systems in cwd.')
        print('   skipped: T to see the number of skipped rxns in cwd.')
        sys.exit()
    else:
        run = sys.argv[1]
        redo = sys.argv[2]
        properties = sys.argv[3]
        prop_redo = sys.argv[4]
        wipe = sys.argv[5]
        skipped = sys.argv[6]

    if run == 'T':
        main_run(redo)
    if wipe == 'T':
        main_wipe()
    if properties == 'T':
        if prop_redo == 'T':
            main_analysis(prop_redo=True)
        elif prop_redo == 'F':
            main_analysis(prop_redo=False)
    if skipped == 'T':
        search_output_dir = os.getcwd()+'/'
        percent_skipped(search_output_dir)


    # out_dir = '/home/atarzia/psp/screening_results/new_reactions/'
    # filename = out_dir+'sRS-3_5_1_32-BRENDA-BR22.bpkl'
    # molecule_db_dir = '/home/atarzia/psp/molecule_DBs/atarzia/'
    # molecules = glob.glob(molecule_db_dir+'ATRS_*.bpkl')
    # rs = get_RS(filename=filename, output_dir=out_dir, verbose=False)
    # rs.__dict__
    # rs.mol_collected = False
    # collect_RS_molecule_properties(rs=rs, output_dir=out_dir,
    #                                mol_db_dir=molecule_db_dir,
    #                                molecules=molecules, count=0,
    #                                react_syst_files=[])
    # rs.__dict__
    # # rs.save_object(rs.pkl)
    # #
    # for m in rs.components:
    #     print(m.name)
    #     print(m.mid_diam)
    #     print(m.SMILES)
    #     print(m.XlogP)
    #     print(m.complexity)
    #     # m.complexity = None
    #     # m.XlogP = None
    #     print(m.pkl)

# for rs in yield_rxn_syst(output_dir='/home/atarzia/psp/screening_results/new_reactions/'):
#     if rs.pkl != 'sRS-4_99_1_1-BRENDA-BR12.bpkl':
#         continue
#     print(rs.pkl)
#     print(rs.pkl.replace('.pkl', '.bpkl'))
#     rs.pkl = rs.pkl.replace('.pkl', '.bpkl')
#     rs.save_object('/home/atarzia/psp/screening_results/biomin_search/'+rs.pkl)
