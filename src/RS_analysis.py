#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse all RS properties.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

import time
from multiprocessing import Pool
import glob
from os import getcwd
import pandas as pd
import sys
from ercollect.molecule import (
    read_molecule_lookup_file,
    search_molecule_by_ident,
    write_lookup_files
)


def collect_RS_molecule_properties(
    rs,
    output_dir,
    mol_db_dir,
    molecules,
    count=0,
    react_syst_files=None
):
    """Collect molecule properties from my database.

    """
    if rs.mol_collected is True:
        return None
    if rs.skip_rxn is True or rs.components is None:
        rs.skip_rxn = True
        rs.save_object(output_dir+rs.pkl)
        return None
    # ignore any reactions with unknown components
    rs.skip_rxn = False
    for m in rs.components:
        try:
            if m.mol is None:
                rs.skip_rxn = True
        except AttributeError:
            rs.skip_rxn = True
    if rs.skip_rxn is True:
        print('skipping reaction - it is incomplete or generic')
        rs.skip_reason = (
            'one component has no molecule -'
            ' rxn is incomplete or generic'
        )
        rs.save_object(output_dir+rs.pkl)
        return None
    print('>', rs.pkl)
    count_found = 0
    # collect properties from molecule DB
    for m in rs.components:
        lookup_file = (
            '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
        )
        dataset = read_molecule_lookup_file(lookup_file=lookup_file)
        print(m.name, m.SMILES)
        existing_pkl = search_molecule_by_ident(molec=m,
                                                dataset=dataset)
        if existing_pkl is not None:
            count_found += 1
            db_mol = load_molecule(existing_pkl, verbose=True)
            # copy DB object properties to RS
            # only overwrite None or NaN
            for key in db_mol.__dict__:
                val = db_mol.__dict__[key]
                if key not in m.__dict__:
                    m.__dict__[key] = val
                elif m.__dict__[key] is None and val is not None:
                    m.__dict__[key] = val
            m.pkl = existing_pkl
            rs.save_object(output_dir+rs.pkl)
    if count_found < len(rs.components):
        # if no match found for at least one molecule
        print('molecule not in database!')
        print('run molecule.py!')
        print('check temp_Failures.txt ---- exitting...')
        with open(output_dir+'temp_Failures.txt', 'a') as f:
            f.write(output_dir+rs.pkl+'\n')
        sys.exit()
    else:
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
    # ignore any reactions with unknown components
    rs.skip_rxn = False
    for m in rs.components:
        if m.mol is None:
            rs.skip_rxn = True
    if rs.skip_rxn is True:
        print('skipping reaction - it is incomplete or generic')
        rs.skip_reason = (
            'one component has no molecule -'
            ' rxn is incomplete or generic'
        )
        rs.save_object(output_dir+rs.pkl)
        return None

    max_comp_size = 0
    for m in rs.components:
        print('>>> ', m.name, ':', m.mid_diam, ':', m.pkl)
        # remove reactions with general atoms (given by '*' in SMILES)
        if "*" in m.SMILES:
            rs.skip_rxn = True
            rs.skip_reason = 'one component has wildcard SMILES'
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
            rs.skip_reason = (
                'one component could not have diameter calculated'
            )
            print(
                'molecule diameter could not'
                'be calculated for some reason.'
            )
            print('skipping...')
            rs.save_object(output_dir+rs.pkl)
            return None
        max_comp_size = max([m.mid_diam, max_comp_size])

    rs.max_comp_size = max_comp_size
    rs.save_object(output_dir+rs.pkl)


def except_from_solubility():
    """Returns a list of molecules (by SMILES) that we do not want to
    include
    in solubility calculations.

    """
    li = ['O', '[H+]', '[OH-]']
    return li


def RS_solubility(rs, output_dir):
    """Check solubility properties of a reaction system using logS.

    Defining solubility in terms of logS
    (https://github.com/PatWalters/solubility).

    Keywords:
        rs (reaction) - reaction system
        output_dir (str) - directory to output molecule files

    """

    if rs.skip_rxn is True:
        rs.min_logS = None
        rs.max_logS = None
        return None
    rs.min_logS = 100
    rs.max_logS = -100
    for m in rs.components:
        if m.SMILES in except_from_solubility():
            continue
        # ignore charged components
        if m.logS == 'charged':
            continue
        # ignore reactions with unknown values
        if m.logS is None or m.logS == 'not found':
            rs.min_logS = 100
            rs.max_logS = -100
            break
        rs.min_logS = min([rs.min_logS, m.logS])
        rs.max_logS = max([rs.max_logS, m.logS])
    if rs.min_logS == 100:
        rs.min_logS = None
    if rs.max_logS == -100:
        rs.max_logS = None
    rs.save_object(output_dir+rs.pkl)


def RS_hphobicity(rs, output_dir):
    """Check solubility properties of a reaction system using logP.

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
        if m.SMILES in except_from_solubility():
            continue
        # ignore charged components
        if m.logP == 'charged':
            sys.exit('there should be none of these')
            continue
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
    """Get the change in maximum synthetic accessibility (sa)
    a reaction system.

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
    # get the max no heavy atoms of reactants and products
    r_max_NHA = 0
    p_max_NHA = 0
    for m in rs.components:
        # ignore reactions with unknown values
        if m.Synth_score is None or m.Synth_score == 'not found':
            rs.r_max_sa = -100
            rs.p_max_sa = -100
            break
        if m.role == 'reactant':
            NHA = m.mol.GetNumHeavyAtoms()
            if NHA == r_max_NHA:
                rs.r_max_sa = max([rs.r_max_sa, m.Synth_score])
            elif NHA > r_max_NHA:
                rs.r_max_sa = m.Synth_score
                r_max_NHA = NHA
        elif m.role == 'product':
            NHA = m.mol.GetNumHeavyAtoms()
            if NHA == p_max_NHA:
                rs.p_max_sa = max([rs.p_max_sa, m.Synth_score])
            elif NHA > p_max_NHA:
                rs.p_max_sa = m.Synth_score
                p_max_NHA = NHA
    if rs.r_max_sa == -100:
        rs.r_max_sa = None
        rs.delta_sa = None
    elif rs.p_max_sa == -100:
        rs.p_max_sa = None
        rs.delta_sa = None
    else:
        rs.delta_sa = rs.p_max_sa - rs.r_max_sa
    rs.save_object(output_dir+rs.pkl)


def RS_hphobicity_X(rs, output_dir):
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
        if m.SMILES in except_from_solubility():
            continue
        # ignore charged components
        if m.XlogP == 'charged':
            continue
        # ignore reactions with unknown values
        if m.XlogP is None or m.XlogP == 'not found':
            rs.min_XlogP = 100
            rs.max_XlogP = -100
            break
        try:
            rs.min_XlogP = min([rs.min_XlogP, float(m.XlogP)])
            rs.max_XlogP = max([rs.max_XlogP, float(m.XlogP)])
        except ValueError:
            molecule_db_dir = '/home/atarzia/psp/molecule_DBs/atarzia/'
            molecules = glob.glob(molecule_db_dir+'ATRS_*.gpkl')
            rs.mol_collected = False
            collect_RS_molecule_properties(
                rs=rs,
                output_dir=output_dir,
                mol_db_dir=molecule_db_dir,
                molecules=molecules,
                count=0,
                react_syst_files=[]
            )
            print('temporary fix for PUBCHEM errors')
            print('fix this!')
            rs.min_XlogP = 100
            rs.max_XlogP = -100

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
                rs.r_max_comp = max([
                    rs.r_max_comp, float(m.complexity)
                ])
            elif m.role == 'product':
                rs.p_max_comp = max([
                    rs.p_max_comp, float(m.complexity)
                ])
        except ValueError:
            molecule_db_dir = '/home/atarzia/psp/molecule_DBs/atarzia/'
            molecules = glob.glob(molecule_db_dir+'ATRS_*.gpkl')
            rs.mol_collected = False
            collect_RS_molecule_properties(
                rs=rs, output_dir=output_dir,
                mol_db_dir=molecule_db_dir,
                molecules=molecules, count=0,
                react_syst_files=[]
            )
            print('temporary fix for PUBCHEM errors')
            print('fix this!')
            rs.r_max_comp = -100000
            rs.p_max_comp = -100000

    if rs.r_max_comp == -100000:
        rs.r_max_comp = None
        rs.delta_comp = None
    elif rs.p_max_comp == -100000:
        rs.p_max_comp = None
        rs.delta_comp = None
    else:
        rs.delta_comp = rs.p_max_comp - rs.r_max_comp
    rs.save_object(output_dir+rs.pkl)


def parallel_analysis(
    rs,
    count,
    react_syst_files,
    output_dir,
    threshold,
    mol_db_dir,
    molecules,
    done_pkls
):
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
        RS_hphobicity(rs=rs, output_dir=output_dir)
        RS_SAscore(rs=rs, output_dir=output_dir)
        RS_hphobicity_X(rs=rs, output_dir=output_dir)
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
    EC_DF = pd.read_table(
        EC_file,
        delimiter='__',
        names=['EC_no', 'description'],
        engine='python'
    )
    search_ECs = list(EC_DF['EC_no'])

    # remove all spaces within EC numbers
    search_ECs = [i.replace(' ', '') for i in search_ECs]

    # add check for '1' from '1.-.-.-'
    new_search_ECs = []
    for EC in search_ECs:
        if '-' in EC:
            new_search_ECs.append(EC.replace('.-', ''))
            new_search_ECs.append(EC)
        else:
            new_search_ECs.append(EC)

    print(len(search_ECs), 'EC numbers to test')
    print('first EC:', search_ECs[0], '---- last EC:', search_ECs[-1])
    print('collect all reaction systems (ONLINE)...')
    return new_search_ECs


def main_analysis(prop_redo, file_list):
    """Analyse all reaction systems.

    """
    print('-----------------------------------------------------')
    print('Collect reaction properties')
    print('-----------------------------------------------------')
    NP = 1  # number of processes
    pI_thresh = 6
    size_thresh = 4.2  # angstroms
    molecule_db_dir = '/home/atarzia/psp/molecule_DBs/atarzia/'
    print('settings:')
    print('    Number of processes:', NP)
    print('    pI threshold:', pI_thresh)
    print('    Diffusion threshold:', size_thresh, 'Angstrom')
    print('    Molecule database:', molecule_db_dir)
    search_output_dir = getcwd()+'/'
    react_syst_files = glob.glob(search_output_dir+'sRS-*.gpkl')
    molecules = glob.glob(molecule_db_dir+'ATRS_*.gpkl')
    lookup_file = '/home/atarzia/psp/molecule_DBs/atarzia/lookup.txt'
    translator = '/home/atarzia/psp/molecule_DBs/KEGG/translator.txt'
    molecule_DB_directory = '/home/atarzia/psp/molecule_DBs/atarzia/'
    # write molecule look up files based on molecule DB
    write_lookup_files(lookup_file, translator, molecule_DB_directory)
    print('------------------------------------------------------')
    print('collect component properties and analyse reaction systems:')
    print('    - diffusion of components')
    print('    - solubility (logP) of components')
    print('    - change in synthetic accessibility of components')
    print('    - solubiity (XlogP) of components')
    print('    - change in complexity of components')
    print('------------------------------------------------------')
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
            # func(rs, count, react_syst_files,
            # search_output_dir, size_thrsh)
            args = [
                (
                    rs,
                    i,
                    react_syst_files,
                    search_output_dir,
                    size_thresh,
                    molecule_db_dir,
                    molecules,
                    done_pkls
                )
                for i, rs in enumerate(
                    yield_rxn_syst(search_output_dir)
                )]
            pool.starmap(parallel_analysis, args)
    # in serial with trivial parallelisation
    else:
        # generator = yield_rxn_syst(search_output_dir, verbose=True)
        generator = yield_rxn_syst_filelist(
            search_output_dir,
            file_list,
            verbose=True
        )
        for i, rs in enumerate(generator):
            if 'KEGG' not in rs.pkl:
                continue
            print('---------------------------------------------')
            print('checking rxn', i, 'of', len(react_syst_files))
            # rs.mol_collected = False
            if rs.pkl not in done_pkls:
                # rs.mol_collected = False
                print('collect molecule properties...')
                collect_RS_molecule_properties(
                    rs=rs,
                    output_dir=search_output_dir,
                    mol_db_dir=molecule_db_dir,
                    molecules=molecules,
                    count=i,
                    react_syst_files=react_syst_files
                )
                if rs.mol_collected is False:
                    continue
                print('check diffusion...')
                RS_diffusion(rs=rs, output_dir=search_output_dir,
                             threshold=size_thresh)
                print('check solubility (logS)...')
                RS_solubility(rs=rs, output_dir=search_output_dir)
                print('check hphobicity (logP)...')
                RS_hphobicity(rs=rs, output_dir=search_output_dir)
                RS_SAscore(rs=rs, output_dir=search_output_dir)
                print('check synthetic accessibility...')
                print('check hphobicity (XlogP)...')
                RS_hphobicity_X(rs=rs, output_dir=search_output_dir)
                print('check complexity...')
                RS_complexity_score(
                    rs=rs,
                    output_dir=search_output_dir
                )
                with open(search_output_dir+'prop_done.txt', 'a') as f:
                    f.write(rs.pkl+'\n')

    print(
        '--- time taken =',
        '{0:.2f}'.format(time.time()-temp_time),
        's'
    )


def main():
    if (not len(sys.argv) == 4):
        print("""
Usage: RS_analysis.py properties rerun_properties prop_file
    properties: T to get properties of reaction systems in cwd.
    rerun properites?: T for rerun, F to read from prop_done.txt.
    prop_file: name of file containing list of RS.
""")
        sys.exit()
    else:
        properties = sys.argv[1]
        prop_redo = sys.argv[2]
        prop_file = sys.argv[3]

    if properties == 'T':
        if prop_redo == 'T':
            main_analysis(prop_redo=True, file_list=prop_file)
        elif prop_redo == 'F':
            main_analysis(prop_redo=False, file_list=prop_file)

    print('All done!')


if __name__ == "__main__":
    main()
