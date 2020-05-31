#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run screening for biomineralisation literature reactions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import sys
import time
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem as Chem

import rdkit_functions as rdkf
import plotting_fn as pfn
import plots_molecular as pm
import utilities


def EC_sets():
    """
    Define target EC numbers and reactant for literature.

    """
    # set EC numbers of interest based on dataset
    # 'none' in species lists implies species was not given with
    # literature
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
        '1.1.1.1': ['none'],
        '1.1.1.27': ['none'],
        '3.2.1.23': ['none'],
        '3.2.1.26': ['none'],
        '3.1.1.3': [
            'thermomyces lanuginosus',
            'alcaligenes sp.',
            'pseudomonas fluorescens',
            'rhizomucor miehei',
            'candida antarctica',
            'aspergillus niger'
        ],
        '3.1.1.6': ['lactobacillus acidophilus'],
        '3.5.1.11': ['none'],
    }
    # set the molecules associated with each EC number in literature
    # data set
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
            (
                'ABTS',
                'CCN1C2=C(C=C(C=C2)S(=O)(=O)[O-])SC1=NN=C3N'
                '(C4=C(S3)C=C(C=C4)S(=O)(=O)[O-])CC'
            ),
        ]},
        '1.9.3.1': {'equus caballus': [
            ('hydrogen peroxide', 'OO'),
            ('water', 'O'),
            ('oxygen', 'O=O'),
            (
                'ABTS',
                'CCN1C2=C(C=C(C=C2)S(=O)(=O)[O-])SC1=NN=C3N'
                '(C4=C(S3)C=C(C=C4)S(=O)(=O)[O-])CC'
            ),
            ('Amplex Red', 'CC(=O)N1C2=C(C=C(C=C2)O)OC3=C1C=CC(=C3)O'),
            ('resorufin', 'C1=CC2=C(C=C1O)OC3=CC(=O)C=CC3=N2'),
            ('methyl ethyl ketone peroxide', 'CCC(C)(OO)OOC(C)(CC)OO'),
            ('tert-butyl hydroperoxide', 'CC(C)(C)OO'),
        ]},
        '1.1.5.2': {'none': [
            (
                'D-glucose (ring)',
                'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
            ),
            # part of phenazine methosulfate
            ('methosulfate', 'COS(=O)(=O)[O-]'),
            # part of phenazine methosulfate
            (
                '5-methylphenazin-5-ium',
                'C[N+]1=C2C=CC=CC2=NC3=CC=CC=C31'
            ),
            (
                '2,6-dichloroindophenol',
                'C1=CC(=O)C=CC1=NC2=CC(=C(C(=C2)Cl)O)Cl'
            ),
        ]},
        '3.5.1.5': {'canavalia ensiformis': [
            ('urea', 'C(=O)(N)N'),
            ('water', 'O'),
            ('carbon dioxide', 'C(=O)=O'),
            ('ammonia', 'N'),
        ]},
        '1.1.3.4': {'aspergillus niger': [
            # ('D-glucose (chain)', 'C(C(C(C(C(C=O)O)O)O)O)O'),
            (
                'D-glucose (ring)',
                'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
            ),
            ('hydrogen peroxide', 'OO'),
            (
                'D-gluconolactone',
                'C([C@@H]1[C@H]([C@@H]([C@H](C(=O)O1)O)O)O)O'
            ),
            ('benzoquinone', 'C1=CC(=O)C=CC1=O'),
        ]},
        '1.13.12.4': {'none': [
            ('hydrogen peroxide', 'OO'),
            ('pyruvate', 'CC(=O)C(=O)[O-]'),
            ('L-lactate', 'CC(C(=O)[O-])[O-]'),
        ]},
        '1.1.1.1': {'none': [
            ('ethanol', 'CCO'),
            ('acetaldehyde', 'CC=O'),
        ]},
        '1.1.1.27': {'none': [
            ('pyruvate', 'CC(=O)C(=O)[O-]'),
            ('L-lactate', 'CC(C(=O)[O-])[O-]'),
        ]},
        '-.-.-.-': {'none': [
            # assume charge and non charged species are same size
            (
                'methylene blue+',
                'CN(C)C1=CC2=C(C=C1)N=C3C=CC(=[N+](C)C)C=C3S2'
            ),
            # (
            #     'fluorescein',
            #     'C1=CC=C2C(=C1)C(=O)OC23C4=C(C=C(C=C4)O)'
            #     'OC5=C3C=CC(=C5)O'
            # ),
            # ('hexanol', 'CCCCCCO'),
            # ('vinyl_acetate', 'CC(=O)OC=C'),
            # ('hexyl_acetate', 'CCCCCCOC(=O)C'),
            # ('1-phenylethanol', 'CC(C1=CC=CC=C1)O'),
            # ('1-phenylethyl_acetate', 'CC(C1=CC=CC=C1)OC(=O)C'),
            # ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
            # (
            #     'p-nitrophenyl butyrate',
            #     'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            # ),
            # ('3-bromo-1-propanol', 'C(CO)CBr'),
            # ('1,3-propanediol', 'C(CO)CO'),
            # ('1,2-dibromoethane', 'C(CBr)Br'),
            # ('2-bromoethanol', 'C(CBr)O'),
            # ('1,2-ethanediol', 'C(CO)O'),
            # ('1,3-dibromopropane', 'C(CBr)CBr'),
            # ('methanol', 'CO'),
            # ('formaldehyde', 'C=O'),
            # ('urea', 'C(=O)(N)N'),
        ]},
        '3.2.1.23': {'none': [
            (
                'beta-lactose',
                'C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]'
                '2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O'
            ),
            (
                'D-glucose (ring)',
                'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
            ),
            (
                'D-galactose (ring)',
                'C([C@@H]1[C@@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
            ),
            ('water', 'O'),
        ]},
        '3.2.1.26': {'none': [
            (
                'sucrose',
                'C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O'
            ),
            (
                'D-glucose (ring)',
                'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
            ),
            (
                'D-fructose (ring)',
                'C([C@@H]1[C@H]([C@@H](C(O1)(CO)O)O)O)O'
            ),
        ]},
        '3.1.1.3': {
            'thermomyces lanuginosus': [
                ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                (
                    'p-nitrophenyl butyrate',
                    'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                ('butyric acid', 'CCCC(=O)O'),
            ],
            # not sure about this one - the kinetic resolution may not
            #  belong to this enyzme?
            'alcaligenes sp.': [
                ('2-octanol', 'CCCCCCC(C)O'),
                ('vinyl acetate', 'CC(=O)OC=C'),
                ('octyl acetate', 'CCCCCCCCOC(=O)C'),
                ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                (
                    'p-nitrophenyl octanoate',
                    'CCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                ('octanoic acid', 'CCCCCCCC(=O)O'),
            ],
            'pseudomonas fluorescens': [
                ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                (
                    'p-nitrophenyl butyrate',
                    'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                ('butyric acid', 'CCCC(=O)O'),
            ],
            'rhizomucor miehei': [
                ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                (
                    'p-nitrophenyl butyrate',
                    'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                ('butyric acid', 'CCCC(=O)O'),
            ],
            'candida antarctica': [
                ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                (
                    'p-nitrophenyl butyrate',
                    'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                ('butyric acid', 'CCCC(=O)O'),
            ],
            'aspergillus niger': [
                ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
                (
                    'p-nitrophenyl acetate',
                    'CC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
                ),
                ('acetic acid', 'CC(=O)O'),
            ]
        },
        '3.1.1.6': {'lactobacillus acidophilus': [
            ('p-nitrophenol', 'C1=CC(=CC=C1[N+](=O)[O-])O'),
            (
                'p-nitrophenyl acetate',
                'CC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            ),
            ('acetic acid', 'CC(=O)O'),
            (
                'p-nitrophenyl phosphate',
                'C1=CC(=CC=C1[N+](=O)[O-])OP(=O)(O)O'
            ),
            ('phosphate acid', 'OP(=O)([O-])[O-]'),
            (
                'p-nitrophenyl butyrate',
                'CCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            ),
            ('butyric acid', 'CCCC(=O)O'),
            (
                'p-nitrophenyl hexanoate',
                'CCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            ),
            ('hexanoic acid', 'CCCCCC(=O)O'),
            (
                'p-nitrophenyl octanoate',
                'CCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            ),
            ('octanoic acid', 'CCCCCCCC(=O)O'),
            (
                'p-nitrophenyl decanoate',
                'CCCCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            ),
            ('decanoic acid', 'CCCCCCCCCC(=O)O'),
            (
                'p-nitrophenyl dodecanoate',
                'CCCCCCCCCCCC(=O)OC1=CC=C(C=C1)[N+](=O)[O-]'
            ),
            ('dodecanoic acid', 'CCCCCCCCCCCC(=O)O'),
        ]},
        '3.5.1.11': {'none': [
            (
                'penicillin-G',
                'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C'
            ),
        ]},
    }

    EC_descriptors = {
        # EC : descriptor
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


def biomin_known(molecules, output_dir, plot_suffix):
    """
    Scatter plot of all molecule sizes in dictionary.

    """
    fig, ax = plt.subplots(figsize=(8, 5))
    m_diams = []
    for name in molecules:
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        print('-----', name, mid_diam, '-----')
        m_diams.append(mid_diam)

    m_diams = np.asarray(m_diams)

    X_bins = np.arange(0.1, 21, 0.5)
    hist, bin_edges = np.histogram(a=m_diams, bins=X_bins)
    ax.bar(
        bin_edges[:-1],
        hist,
        align='edge',
        width=0.5,
        color='#2C3E50',
        edgecolor='k',
        alpha=0.8
    )
    ax.axvline(x=3.4, c='k')
    ax.axvspan(
        xmin=4.0,
        xmax=6.6,
        facecolor='k',
        alpha=0.25,
        # hatch="/"
    )
    # ax.axvspan(xmin=5.4, xmax=6.6, facecolor='k', alpha=0.2)
    pfn.define_standard_plot(
        ax,
        # xtitle='intermediate diameter [$\mathrm{\AA}$]',
        xtitle=r'$d$ [$\mathrm{\AA}$]',
        ytitle='count',
        xlim=(0, 15),
        ylim=(0, 15)
    )
    fig.tight_layout()
    fig.savefig(
        f"molecule_size_{plot_suffix}.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def n_phenyl_assay(output_dir):
    """
    Prepare figure showing the change in intermediate diameter for
    molecules
    commonly used in n-phenyl ester hydrolysis assays.

    """
    # the n-phenyl esters
    mol_list_1 = [
        'p-nitrophenyl acetate',
        'p-nitrophenyl butyrate',
        'p-nitrophenyl hexanoate',
        'p-nitrophenyl octanoate',
        'p-nitrophenyl decanoate',
        'p-nitrophenyl dodecanoate'
    ]

    # the products
    mol_list_2 = [
        'acetic acid',
        'butyric acid',
        'hexanoic acid',
        'octanoic acid',
        'decanoic acid',
        'dodecanoic acid'
    ]

    # no Cs
    no_Cs = [2, 4, 6, 8, 10, 12]

    fig, ax = plt.subplots()
    for i, C in enumerate(no_Cs):
        # ester
        name = mol_list_1[i]
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        print(C, mol_list_1[i], mid_diam)
        ax.scatter(
            C, mid_diam,
            c='r',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=100
        )
        # acid
        name = mol_list_2[i]
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        print(C, mol_list_2[i], mid_diam)
        ax.scatter(
            C, mid_diam,
            c='b',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=120
        )
    ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2, hatch="/")
    # n-phenol
    name = 'p-nitrophenol'
    out_file = (
        f"{output_dir}/"
        f"{name.replace(' ', '_').replace('/', '__')}"
        '_diam_result.csv'
    )
    if os.path.exists(out_file) is False:
        import sys
        sys.exit('calc molecule diameters!')
    results = pd.read_csv(out_file)
    mid_diam = min(results['diam2'])
    ax.axhline(y=mid_diam, c='purple', alpha=1)
    pfn.define_standard_plot(
        ax,
        xtitle='no. carbons',
        ytitle=r'$d$ [$\mathrm{\AA}$]',
        xlim=(1, 14),
        ylim=(2.5, 8)
    )
    # decoy legend
    ax.scatter(
        -100, -100, c='r',
        edgecolors='k', marker='o', alpha=1.0,
        s=100, label='ester'
    )
    ax.scatter(
        -100, -100, c='b',
        edgecolors='k', marker='o', alpha=1.0,
        s=100, label='acid'
    )
    ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig(
        "ester_comp.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def cyt_C_perox_assay(output_dir):
    """
    Prepare figure showing the change in intermediate diameter for 3
    peroxide molcules degraded by Cyt-C in ZIF-8 (One-Pot Synthesis of
    Protein-Embedded Metal–Organic Frameworks with Enhanced Biological
    Activities, DOI:10.1021/nl5026419)

    """
    # the n-phenyl esters
    mol_list_1 = [
        'hydrogen peroxide',
        'methyl ethyl ketone peroxide',
        'tert-butyl hydroperoxide'
    ]
    smiles_list_1 = [
        'OO',
        'CCC(C)(OO)OOC(C)(CC)OO',
        'CC(C)(C)OO'
    ]
    fig, ax = plt.subplots()
    for i, name in enumerate(mol_list_1):
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles_list_1[i]))
        MW = Descriptors.MolWt(mol)
        print(name, mol_list_1[i], MW, mid_diam)
        ax.scatter(
            MW,
            mid_diam,
            c='k',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=100
        )

    ax.axhspan(ymin=4.0, ymax=6.6, facecolor='k', alpha=0.2, hatch="/")
    pfn.define_standard_plot(
        ax,
        xtitle='molecular weight [g/mol]',
        ytitle=r'$d$ [$\mathrm{\AA}$]',
        xlim=(10, 250),
        ylim=(2.5, 8)
    )
    fig.tight_layout()
    fig.savefig(
        "cytC_comp.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def HOF_examples(output_dir):
    """
    Prepare figure showing the value of d for all molecules used in the
    BioHOFs from: 10.1021/jacs.9b06589

    """
    # the n-phenyl esters
    mol_list_1 = [
        'fluorescein',
        'hydrogen_peroxide',
        'methanol',
        'formaldehyde',
        'urea'
    ]
    smiles_list_1 = [
        'C1=CC=C2C(=C1)C(=O)OC23C4=C(C=C(C=C4)O)OC5=C3C=CC(=C5)O',
        'OO',
        'CO',
        'C=O',
        'C(=O)(N)N'
    ]
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, name in enumerate(mol_list_1):
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles_list_1[i]))
        MW = Descriptors.MolWt(mol)
        print(name, mol_list_1[i], MW, mid_diam)
        ax.scatter(
            MW,
            mid_diam,
            c='#5499C7',
            edgecolors='k',
            marker='o',
            alpha=1.0,
            s=140
        )

    # ax.axhline(y=11.8, c='k', alpha=0.2)
    pfn.define_standard_plot(
        ax,
        xtitle='molecular weight [g/mol]',
        ytitle=r'$d$ [$\mathrm{\AA}$]',
        xlim=(10, 500),
        ylim=(2.5, 15)
    )
    fig.tight_layout()
    fig.savefig(
        "HOF_examples.pdf",
        dpi=720,
        bbox_inches='tight'
    )


def Tash_esters(output_dir):
    """


    """
    # the n-phenyl esters
    mol_list_1 = [
        'hexanol',
        'vinyl_acetate',
        'hexyl_acetate',
        '1-phenylethanol',
        '1-phenylethyl_acetate',
        'p-nitrophenol',
        'p-nitrophenyl butyrate',
        '3-bromo-1-propanol',
        '1,3-propanediol',
        '1,2-dibromoethane',
        '2-bromoethanol',
        '1,2-ethanediol',
        '1,3-dibromopropane'
    ]
    print('===== Tash examples =====')
    for i, name in enumerate(mol_list_1):
        out_file = (
            f"{output_dir}/"
            f"{name.replace(' ', '_').replace('/', '__')}"
            '_diam_result.csv'
        )
        if os.path.exists(out_file) is False:
            continue
        results = pd.read_csv(out_file)
        mid_diam = min(results['diam2'])
        print(f'{name}: d = {round(mid_diam, 2)} Angstrom')


def get_molecule_DB(EC_mol_set, output_dir):
    """
    Get molecule dictionary + output 2D structures.

    """

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    molecules = {}
    diameters = {}
    for i in EC_mol_set.keys():
        for j in EC_mol_set[i].keys():
            for mol in EC_mol_set[i][j]:
                molecules[mol[0]] = mol[1]
                diameters[mol[0]] = 0
                print(mol[0]+'&'+mol[1])
                rdkf.draw_smiles_to_svg(
                    mol[1],
                    os.path.join(
                        output_dir,
                        mol[0].replace(' ', '_')+'_2d.svg'
                    )
                )
    return molecules, diameters


def main():
    if (not len(sys.argv) == 3):
        print("""
    Usage: biomin_screening.py

        rerun_diameter_calc

        param_file

        """)
        sys.exit()
    else:
        rerun_diameter_calc = True if sys.argv[1] == 't' else False
        pars = utilities.read_params(sys.argv[2])

    start = time.time()
    # set parameters
    EC_set, EC_mol_set, EC_descriptors = EC_sets()

    # Ignore MW restrictions.
    pars['MW_thresh'] = 2000

    print('------------------------------------------------')
    print('Screen molecular size of compounds in known reactions')
    print('------------------------------------------------')

    # screen known reactant and product molecules
    print('--- get molecule DB + draw 2D structures...')
    molecules, diameters = get_molecule_DB(
        EC_mol_set=EC_mol_set,
        output_dir='2d_'
    )

    # calculate the size of the ellipsoid surroudning all molecules
    if rerun_diameter_calc:
        print('--- calculate molecular diameters...')
        rdkf.calc_molecule_diameters(
            molecules,
            pars=pars,
            out_dir='biomin_sizes',
        )

    # print results for each molecule
    print('--- print results and plot...')
    pm.print_results(
        molecules,
        threshold=pars['size_thresh'],
        output_dir='biomin_sizes'
    )

    # plotting
    biomin_known(
        molecules,
        output_dir='biomin_sizes',
        plot_suffix='biomin_known'
    )
    pm.categorical(
        molecules,
        threshold=pars['size_thresh'],
        output_dir='biomin_sizes',
        plot_suffix='biomin_known'
    )
    pm.shapes(
        molecules,
        threshold=pars['size_thresh'],
        output_dir='biomin_sizes',
        plot_suffix='biomin_known'
    )

    n_phenyl_assay(output_dir='biomin_sizes')
    cyt_C_perox_assay(output_dir='biomin_sizes')
    HOF_examples(output_dir='biomin_sizes')
    Tash_esters(output_dir='biomin_sizes')

    end = time.time()
    print(f'---- total time taken = {round(end-start, 2)} s')


if __name__ == "__main__":
    main()
