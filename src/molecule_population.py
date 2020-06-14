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
from rdkit.Chem.GraphDescriptors import BertzCT

import IO
import rdkit_functions as rdkf
import plots_molecular as pm
import utilities
import chemcost_IO


def populate_all_molecules(
    params,
    redo_size,
    redo_prop,
    mol_file=None
):
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
    print(mol_file)

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
        if not exists(prop_file) or redo_prop:
            print(
                f'>> calculating molecule descriptors of {mol}, '
                f'{count} of {len(molecule_list)}'
            )
            prop_dict = {}
            rdkitmol = Chem.MolFromSmiles(smiles)
            rdkitmol.Compute2DCoords()
            Chem.SanitizeMol(rdkitmol)
            prop_dict['logP'] = Descriptors.MolLogP(
                rdkitmol,
                includeHs=True
            )
            prop_dict['logS'] = rdkf.get_logSw(rdkitmol)
            prop_dict['Synth_score'] = rdkf.get_SynthA_score(rdkitmol)
            prop_dict['NHA'] = rdkitmol.GetNumHeavyAtoms()
            prop_dict['MW'] = Descriptors.MolWt(rdkitmol)
            nrb = CalcNumRotatableBonds(rdkitmol)
            nb = rdkitmol.GetNumBonds(onlyHeavy=True)
            prop_dict['NRB'] = nrb
            if nb == 0:
                prop_dict['NRBr'] = 0.0
            else:
                prop_dict['NRBr'] = nrb / nb
            prop_dict['bertzCT'] = BertzCT(rdkitmol)
            prop_dict['purchasability'] = chemcost_IO.is_purchasable(
                name=name,
                smiles=smiles
            )

            with open(prop_file, 'w') as f:
                json.dump(prop_dict, f)

        # Get a 3D representation of all molecules using ETKDG.
        chk1 = (not exists(opt_file) or redo_size)
        if chk1 and not exists(etkdg_fail):
            print(
                f'>> optimising {mol}, '
                f'{count} of {len(molecule_list)}'
            )
            rdkit_mol = rdkf.ETKDG(mol, seed=seed)
            if rdkit_mol is not None:
                rdkf.write_structure(opt_file, rdkit_mol)

        # Only property to determine at the moment is the molecular
        # size. This produces a csv for all conformers, which will be
        # used in the analysis.
        chk2 = (not exists(diam_file) or redo_size)
        if chk2 and not exists(etkdg_fail):
            print(
                f'>> getting molecular size of {mol}, '
                f'{count} of {len(molecule_list)}'
            )
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
    if (not len(sys.argv) == 6):
        print("""
Usage: molecule_population.py param_file redo mol_file
    param_file:
    redo size:
        t to overwrite SIZE of all molecules.
    redo rest:
        t to overwrite properties of all molecules.
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
        redo_size = True if sys.argv[2] == 't' else False
        redo_prop = True if sys.argv[3] == 't' else False
        plot = True if sys.argv[4] == 't' else False
        mol_file = None if sys.argv[5] == 'f' else sys.argv[5]

    print('settings:')
    print('    Molecule file:', mol_file)
    print(
        'populate the properties attributes for all '
        'molecules in DB...'
    )

    populate_all_molecules(
        params=params,
        mol_file=mol_file,
        redo_size=redo_size,
        redo_prop=redo_prop,
    )

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
        pm.mol_parity(
            propx='NHA',
            propy='logP',
            file=f"NHAvslogP_{params['file_suffix']}",
            xtitle='no. heavy atoms',
            ytitle='logP'
        )
        pm.mol_parity(
            propx='NHA',
            propy='logS',
            file=f"NHAvslogS_{params['file_suffix']}",
            xtitle='no. heavy atoms',
            ytitle='logS'
        )
        pm.mol_categ(
            propx='purchasability',
            propy='Synth_score',
            file=f"purchvsSA_{params['file_suffix']}",
            xtitle='is purchasable',
            ytitle='SAScore'
        )
        pm.mol_categ(
            propx='purchasability',
            propy='bertzCT',
            file=f"purchvsbCT_{params['file_suffix']}",
            xtitle='is purchasable',
            ytitle='BertzCT'
        )
        pm.mol_categ(
            propx='purchasability',
            propy='size',
            file=f"purchvssize_{params['file_suffix']}",
            xtitle='is purchasable',
            ytitle=r'$d$ [$\mathrm{\AA}$]'
        )
        pm.mol_all_dist(plot_suffix=params['file_suffix'])


if __name__ == "__main__":
    main()
