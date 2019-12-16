#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to populate the properties of all molecules in database.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""

from os.path import exists
import sys
import glob
import json
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

import IO
import rdkit_functions as rdkf
import plots_molecular as pm
import utilities


def populate_all_molecules(params, redo, mol_file=None):
    """
    Populate all molecules in pickle files in directory.

    """

    vdwScale = params['vdwScale']
    boxMargin = params['boxMargin']
    spacing = params['spacing']
    show_vdw = params['show_vdw']
    plot_ellip = params['plot_ellip']
    N_conformers = int(params['N_conformers'])
    MW_thresh = params['MW_thresh']
    seed = int(params['seed'])

    fail_list = IO.fail_list_read(
        directory=params['molec_dir'],
        file_name='failures.txt'
    )

    if mol_file is None:
        molecule_list = glob.glob('*_unopt.mol')
    else:
        molecule_list = IO.read_molecule_list(mol_file)

    print(f'{len(molecule_list)} molecules in DB.')

    count = 0
    for mol in sorted(molecule_list):
        count += 1

        name = mol.replace('_unopt.mol', '')
        if name in fail_list:
            continue
        print('--------------------------------------------------')
        print(f'populating {mol}, {count} of {len(molecule_list)}')

        opt_file = name+'_opt.mol'
        etkdg_fail = name+'_unopt.ETKDGFAILED'
        diam_file = name+'_size.csv'
        prop_file = name+'_prop.json'
        smiles = rdkf.read_structure_to_smiles(mol)

        # Check for generics.
        if '*' in smiles:
            IO.fail_list_write(
                new_name=name,
                directory=params['molec_dir'],
                file_name='failures.txt'
            )
            fail_list.append(name)
            continue

        # Get molecular properties from 2D structure.
        if not exists(prop_file) or redo:
            print('>> calculating molecule descriptors')
            prop_dict = {}
            rdkitmol = Chem.MolFromSmiles(smiles)
            rdkitmol = Chem.AddHs(rdkitmol)
            rdkitmol.Compute2DCoords()
            prop_dict['logP'] = Descriptors.MolLogP(
                rdkitmol,
                includeHs=True
            )
            prop_dict['logS'] = rdkf.get_logSw(rdkitmol)
            prop_dict['Synth_score'] = rdkf.get_SynthA_score(rdkitmol)
            prop_dict['NHA'] = rdkitmol.GetNumHeavyAtoms()
            prop_dict['NRB'] = CalcNumRotatableBonds(rdkitmol)
            with open(prop_file, 'w') as f:
                json.dump(prop_dict, f)

        # Get a 3D representation of all molecules using ETKDG.
        if (not exists(opt_file) or redo) and not exists(etkdg_fail):
            print('>> optimising molecule')
            rdkit_mol = rdkf.ETKDG(mol, seed=seed)
            if rdkit_mol is not None:
                rdkf.write_structure(opt_file, rdkit_mol)

        # Only property to determine at the moment is the molecular
        # size. This produces a csv for all conformers, which will be
        # used in the analysis.
        if (not exists(diam_file) or redo) and not exists(etkdg_fail):
            print('>> getting molecular size')
            _ = rdkf.calc_molecule_diameter(
                name,
                smiles,
                out_file=diam_file,
                vdwScale=vdwScale,
                boxMargin=boxMargin,
                spacing=spacing,
                MW_thresh=MW_thresh,
                show_vdw=show_vdw,
                plot_ellip=plot_ellip,
                N_conformers=N_conformers,
                rSeed=seed
            )
            del _


def main():
    if (not len(sys.argv) == 5):
        print("""
Usage: molecule_population.py param_file redo mol_file
    param_file:
    redo:
        t to overwrite all molecules.
    plot:
        t to plot distributions of molecule properties.
    mol_file :
        file name of list of molecules to allow for trivial
        parallelisation, `f` if not specified, where all `mol`
        files are populated.
""")
        sys.exit()
    else:
        params = utilities.read_params(sys.argv[1])
        redo = True if sys.argv[2] == 't' else False
        plot = True if sys.argv[3] == 't' else False
        mol_file = None if sys.argv[4] == 'f' else sys.argv[4]

    print('settings:')
    print('    Molecule file:', mol_file)
    print(
        'populate the properties attributes for all '
        'molecules in DB...'
    )

    populate_all_molecules(params=params, mol_file=mol_file, redo=redo)

    if plot:
        pm.mol_parity(
            propx='logP',
            propy='logS',
            file='logPvslogS',
            xtitle='logP',
            ytitle='logS'
        )
        pm.mol_parity(
            propx='NHA',
            propy='Synth_score',
            file=f"NHAvsSA_{params['file_suffix']}",
            xtitle='no. heavy atoms',
            ytitle='SAScore'
        )
        pm.mol_all_dist(plot_suffix=params['file_suffix'])


if __name__ == "__main__":
    main()
