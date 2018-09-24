#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for biomineralisation literature reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""
import time
import pi_fn
import rdkit_functions
import plotting
from multiprocessing import Pool
import rxn_syst

# script/data set specific functions


def EC_sets():
    """Define target EC numbers and reactant for literature.

    """
    # set EC numbers of interest based on dataset
    # 'none' in species lists implies species was not given with literature
    # report.
    EC_set = {
        # EC : species list
        '1.11.1.5': ['equus caballus', ],
        '1.11.1.6': ['bos taurus', ],
        '1.11.1.7': ['armoracia rusticana', ],
        '1.9.3.1': ['equus caballus', ],
        '1.1.5.2': ['none', ],
        '3.5.1.5': ['canavalia ensiformis'],
        '1.1.3.4': ['aspergillus niger'],
        '1.13.12.4': ['none'],
        '3.2.1.26': ['none'],
        '3.1.1.3': ['thermomyces lanuginosus', 'alcaligenes sp.',
                    'pseudomonas fluorescens',
                    'rhizomucor miehei', 'candida antarctica',
                    'aspergillus niger'],
        '3.1.1.6': ['lactobacillus acidophilus'],
        '3.5.1.11': ['none'],
        }
    # set the molecules associated with each EC number in literature data set
    EC_mol_set = {
        # EC : species list : unique molecules
        '1.11.1.5': {'equus caballus': [
            ('hydrogen peroxide', 'OO'),
            ('pyrogallol', 'C1=CC(=C(C(=C1)O)O)O'),
            ('purporogallin', 'C1=CC(=O)C(=C2C(=C1)C=C(C(=C2O)O)O)O'),
        ]},
        '1.11.1.6': {'bos taurus': [
            ('hydrogen peroxide', 'OO'),
            ('water', 'O'),
            ('oxygen', 'O=O'),
            ('3-amino-1,2,4-triazole', 'C1=NNC(=N1)N'),
        ]},
        '1.11.1.7': {'armoracia rusticana': [
            ('pyrogallol', 'C1=CC(=C(C(=C1)O)O)O'),
            ('purporogallin', 'C1=CC(=O)C(=C2C(=C1)C=C(C(=C2O)O)O)O'),
            ('hydrogen peroxide', 'OO'),
            ('water', 'O'),
            ('oxygen', 'O=O'),
            ('ABTS', 'CCN1C2=C(C=C(C=C2)S(=O)(=O)[O-])SC1=NN=C3N(C4=C(S3)C=C(C=C4)S(=O)(=O)[O-])CC'),
        ]},
        '1.9.3.1': {'equus caballus': [
            ('hydrogen peroxide', 'OO'),
            ('water', 'O'),
            ('oxygen', 'O=O'),
            ('ABTS', 'CCN1C2=C(C=C(C=C2)S(=O)(=O)[O-])SC1=NN=C3N(C4=C(S3)C=C(C=C4)S(=O)(=O)[O-])CC'),
            ('Amplex Red', 'CC(=O)N1C2=C(C=C(C=C2)O)OC3=C1C=CC(=C3)O'),
            ('resorufin', 'C1=CC2=C(C=C1O)OC3=CC(=O)C=CC3=N2'),
            ('methyl ethyl ketone peroxide', 'CCC(C)(OO)OOC(C)(CC)OO'),
            ('tert-butyl hydroperoxide', 'CC(C)(C)OO'),
        ]},
        '1.1.5.2': {'none': [
            ('methosulfate', 'COS(=O)(=O)[O-]'),
            ('5-Methylphenazin-5-ium', 'C[N+]1=C2C=CC=CC2=NC3=CC=CC=C31'),
            ('2,6-dichloroindophenol', 'C1=CC(=O)C=CC1=NC2=CC(=C(C(=C2)Cl)O)Cl'),
        ]},
        '3.5.1.5': {'canavalia ensiformis': [
            ('Urea', 'C(=O)(N)N'),
            ('water', 'O'),
            ('carbon dioxide', 'C(=O)=O'),
            ('ammonia', 'N'),
        ]},
        '1.1.3.4': {'aspergillus niger': [
            ('D-glucose (chain)', 'C(C(C(C(C(C=O)O)O)O)O)O'),
            ('D-glucose (ring)', 'C(C1C(C(C(C(O1)O)O)O)O)O'),
            ('gluconic acid', 'C(C(C(C(C(C(=O)O)O)O)O)O)O'),
            ('hydrogen peroxide', 'OO'),
            ('Gluconolactone', 'C(C1C(C(C(C(=O)O1)O)O)O)O'),
        ]},
        '1.13.12.4': {'none': [
            ('hydrogen peroxide', 'OO'),
            ('pyruvate', 'CC(=O)C(=O)[O-]'),
            ('L-lactate', 'CC(C(=O)[O-])[O-]'),

        ]},
        '3.2.1.26': {'none': [
            ('sucrose', 'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O'),
            ('L-fructose', 'C(C(C(C(C(=O)CO)O)O)O)O'),
            ('D-fructose', 'C1C(C(C(C(O1)(CO)O)O)O)O'),
        ]},
        '3.1.1.3': {'thermomyces lanuginosus': [
                        ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                        ('p-nitrophenyl butyrate', 'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('butyric acid', 'CCCC(=O)O'),
                    ],
                    'alcaligenes sp.': [  # not sure about this one - the kinetic resolution may not belong to this enyzme?
                        ('2-octanol', 'CCCCCCC(C)O'),
                        ('vinyl acetate', 'CC(=O)OC=C'),
                        ('octyl acetate', 'CCCCCCCCOC(=O)C'),
                        ('p-nitrophenyl', 'C1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('p-nitrophenyl octanoate', 'CCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('octanoic acid', 'CCCCCCCC(=O)O'),
                    ],
                    'pseudomonas fluorescens': [
                        ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                        ('p-nitrophenyl butyrate', 'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('butyric acid', 'CCCC(=O)O'),
                    ],
                    'rhizomucor miehei': [
                        ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                        ('p-nitrophenyl butyrate', 'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('butyric acid', 'CCCC(=O)O'),
                    ],
                    'candida antarctica': [
                        ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                        ('p-nitrophenyl butyrate', 'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('butyric acid', 'CCCC(=O)O'),
                    ],
                    'aspergillus niger': [
                        ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                        ('p-nitrophenyl acetate', 'CC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
                        ('acetic acid', 'CC(=O)O'),
                    ]},
        '3.1.1.6': {'lactobacillus acidophilus': [
            ('p-nitrophenyl', 'C1=CC=C(C=C1)[N+](=O)[O-]'),
            ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
            ('p-nitrophenyl acetate', 'CC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
            ('acetic acid', 'CC(=O)O'),
            ('p-nitrophenyl phosphate', 'C1=CC(=CC=C1[N+](=O)[O-])OP(=O)(O)O'),
            ('phosphate acid', 'OP(=O)([O-])[O-]'),
            ('p-nitrophenyl butyrate', 'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
            ('butyric acid', 'CCCC(=O)O'),
            ('p-nitrophenyl hexanoate', 'CCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
            ('hexanoic acid', 'CCCCCC(=O)O'),
            ('p-nitrophenyl octanoate', 'CCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
            ('octanoic acid', 'CCCCCCCC(=O)O'),
            ('p-nitrophenyl decanoate', 'CCCCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
            ('decanoic acid', 'CCCCCCCCCC(=O)O'),
            ('p-nitrophenyl dodecanoate', 'CCCCCCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'),
            ('dodecanoic acid', 'CCCCCCCCCCCC(=O)O'),
        ]},
        '3.5.1.11': {'none': [
            ('penicillin-G', 'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C'),
        ]},
        }
    EC_descriptors = {
        # EC : species list
        '1.11.1.5': '',
        '1.11.1.6': '',
        '1.11.1.7': '',
        '1.9.3.1': '',
        '1.1.5.2': '',
        '3.5.1.5': '',
        '1.1.3.4': '',
        '1.13.12.4': '',
        '3.2.1.26': '',
        '3.1.1.3': '',
        '3.1.1.6': '',
        '3.5.1.11': '',
        }

    return EC_set, EC_mol_set, EC_descriptors


def get_molecule_DB(EC_mol_set, output_dir):
    """Get molecule dictionary + output 2D structures.

    """
    molecules = {}
    diameters = {}
    for i in EC_mol_set.keys():
        for j in EC_mol_set[i].keys():
            for mol in EC_mol_set[i][j]:
                molecules[mol[0]] = mol[1]
                diameters[mol[0]] = 0
                rdkit_functions.draw_smiles_to_svg(
                    mol[1], output_dir+mol[0].replace(' ', '_')+'_2d.svg')
    return molecules, diameters


if __name__ == "__main__":
    start = time.time()
    # set parameters
    EC_set, EC_mol_set, EC_descriptors = EC_sets()
    # pI
    pI_DB_dir = '/home/atarzia/psp/sequence_db/bio_min_dataset/'
    pI_output_dir = pI_DB_dir
    pI_csv = "output_data_pi.csv"
    redo_pI = False
    redo_pI_plots = False
    pI_thresh = 6
    # known molecule screening
    mol_DB_dir = '/home/atarzia/psp/screening_results/biomin_known/'
    mol_output_dir = mol_DB_dir
    vdwScale = 0.8
    boxMargin = 4.0
    spacing = 0.6
    show_vdw = False
    plot_ellip = False
    N_conformers = 50
    MW_thresh = 2000
    size_thresh = 4.2
    rerun_diameter_calc = False
    # reaction search
    search_DBs = ['BRENDA', 'SABIO', 'KEGG', 'BKMS', ]
    search_output_dir = '/home/atarzia/psp/screening_results/biomin_search/'
    search_ECs = ['1.11.1.5', '1.11.1.6', '1.11.1.7', '1.9.3.1',
                  '1.1.5.2', '3.5.1.5', '1.1.3.4', '1.13.12.4',
                  '3.2.1.26', '3.1.1.3', '3.1.1.6', '3.5.1.11']
    search_mol_output_file = search_output_dir+'screening_output.csv'
    search_MW_thresh = 250
    search_run = False
    search_redo = False
    NP = 2  # number of processes
    print('------------------------------------------------------------------')
    print('run parameters:')
    print('pI database dir:', pI_DB_dir)
    print('pI output dir:', pI_output_dir)
    print('Redo pI screening?:', redo_pI)
    print('Redo pI plotting?:', redo_pI_plots)
    print('molecule database dir:', mol_DB_dir)
    print('molecule output dir:', mol_output_dir)
    print('VDW scale:', vdwScale)
    print('Box Margin:', boxMargin, 'Angstrom')
    print('Grid spacing:', spacing, 'Angstrom')
    print('show VDW?:', show_vdw)
    print('Plot Ellipsoid?:', plot_ellip)
    print('No Conformers:', N_conformers)
    print('MW threshold:', MW_thresh, 'g/mol')
    print('pI threshold:', pI_thresh)
    print('Diffusion threshold:', size_thresh, 'Angstrom')
    print('Rerun diameter calculation?:', rerun_diameter_calc)
    print('Search output dir:', search_output_dir)
    print('Search molecule output file:', search_mol_output_file)
    print('Search MW threshold:', search_MW_thresh, 'g/mol')
    print('Run Search?:', search_run)
    print('Redo Search?:', search_redo)
    print('------------------------------------------------------------------')

    print('------------------------------------------------------------------')
    print('Screen pIs')
    print('------------------------------------------------------------------')
    temp_time = time.time()
    # prepare pI calculations
    database_names = pi_fn.prepare_pI_calc(database_directory=pI_DB_dir,
                                           redo_pi=redo_pI,
                                           output_dir=pI_output_dir,
                                           csv=pI_csv)
    # screen protein sequence from EC numbers
    print('--- calculate all pIs for target EC sequences...')
    pi_fn.screen_pIs(database_names, redo_pI_plots=redo_pI_plots,
                     redo_pI=redo_pI, pI_csv=pI_csv,
                     pI_output_dir=pI_output_dir, cutoff_pi=pI_thresh,
                     descriptors=EC_descriptors)

    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    print('------------------------------------------------------------------')
    print('Screen molecular size of compounds in known reactions')
    print('------------------------------------------------------------------')
    temp_time = time.time()
    # screen known reactant and product molecules
    print('--- get molecule DB + draw 2D structures...')
    molecules, diameters = get_molecule_DB(EC_mol_set=EC_mol_set,
                                           output_dir=mol_output_dir)

    # calculate all Molecule Weights
    print('--- calculate MWs...')
    rdkit_functions.calculate_all_MW(molecules)

    # calculate the size of the ellipsoid surroudning all molecules
    print('--- calculate molecular diameters...')
    rdkit_functions.calc_molecule_diameters(
                    molecules, out_dir=mol_output_dir, vdwScale=vdwScale,
                    boxMargin=boxMargin, spacing=spacing,
                    show_vdw=show_vdw, plot_ellip=plot_ellip,
                    N_conformers=N_conformers, MW_thresh=MW_thresh,
                    rerun=rerun_diameter_calc)

    # print results for each molecule
    print('--- print results and plot...')
    plotting.print_results(molecules,
                           threshold=size_thresh,
                           output_dir=mol_output_dir)

    # plotting
    plotting.categorical(molecules,
                         threshold=size_thresh,
                         output_dir=mol_output_dir)
    plotting.shapes(molecules,
                    threshold=size_thresh,
                    output_dir=mol_output_dir)

    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    print('------------------------------------------------------------------')
    print('Screen new reactions')
    print('------------------------------------------------------------------')
    temp_time = time.time()
    print('collect all reaction systems (ONLINE)...')
    if search_run is True:
        for DB in search_DBs:
            # iterate over EC numbers of interest
            # Create a multiprocessing Pool
            with Pool(NP) as pool:
                # process data_inputs iterable with pool
                # func(EC, DB, search_output_dir, search_redo, verbose)
                args = [(EC, DB, search_output_dir, search_redo, True)
                        for EC in search_ECs]
                pool.starmap(rxn_syst.process_collection, args)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    rxn_syst.percent_skipped(search_output_dir)
    print('check all reaction systems for diffusion of components (ONLINE)...')
    rxn_syst.check_all_RS_diffusion(output_dir=search_output_dir,
                                    mol_output_file=search_mol_output_file,
                                    threshold=size_thresh,
                                    vdwScale=vdwScale,
                                    boxMargin=boxMargin, spacing=spacing,
                                    N_conformers=N_conformers,
                                    MW_thresh=search_MW_thresh)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    rxn_syst.percent_skipped(search_output_dir)
    print('get subset of reactions with known protein sequences...')
    rxn_syst.check_all_seedMOF(search_output_dir, pI_thresh)
    rxn_syst.percent_w_sequence(search_output_dir)
    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    print('--- print results and plot...')
    # plot distribution of pI of all known sequences
    plotting.rs_pI_distribution(output_dir=search_output_dir,
                                cutoff_pI=pI_thresh)
    # plot max component size vs pI
    plotting.rs_size_vs_pI(output_dir=search_output_dir,
                           cutoff_pI=pI_thresh,
                           size_thresh=size_thresh)
    # plot number of new reactions as a function of size threshold
    plotting.number_rxns_vs_size(output_dir=search_output_dir,
                                 size_thresh=size_thresh,
                                 pI_thresh=pI_thresh)
    # categorize all molecules in mol output file
    plotting.categorical_moloutput(mol_output_file=search_mol_output_file,
                                   threshold=size_thresh,
                                   output_dir=search_output_dir)
    # print new reactions
    plotting.print_new_rxns(output_dir=search_output_dir)
    # plot a distribution of the change in molecule size due to reaction
    plotting.rs_delta_size(output_dir=search_output_dir)

    print('---- step time taken =', '{0:.2f}'.format(time.time()-temp_time),
          's')
    temp_time = time.time()
    end = time.time()
    print('---- total time taken =', '{0:.2f}'.format(end-start), 's')
