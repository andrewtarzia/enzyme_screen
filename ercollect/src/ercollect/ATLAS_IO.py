#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of ATLAS DB.

Author: Andrew Tarzia

Date Created: 02 Oct 2018

"""
import os
import pandas as pd
from ercollect import DB_functions
from ercollect import rxn_syst
from ercollect.molecule import molecule, iterate_rs_components, load_molecule
from ercollect.KEGG_IO import check_translator, KEGGID_to_CHEBIID


def get_rxn_systems(EC, output_dir, molecule_dataset,
                    clean_system=False, verbose=False):
    """Get reaction systems from KEGG entries in one EC and output to Pickle.

    """
    top_tier = EC.split('.')[0]
    DB_prop = DB_functions.get_DB_prop('ATLAS')
    # read in JSON file of whole DB
    rxn_DB_file = DB_prop[0]+DB_prop[1]['ATLAS_CSV_'+top_tier]
    no_rxns = DB_prop[2]['ATLAS_CSV_'+top_tier]
    # read in chunks
    ATLAS = pd.read_csv(rxn_DB_file, chunksize=1000)
    # iterate over chunks over reactions
    count = 0
    for chunk in ATLAS:
        for idx, row in chunk.iterrows():
            ATLAS_ID = row['ATLAS']
            rxn_ECs = [i.lstrip().rstrip() for i in row['REACTIONRULE'].split("|")]
            if EC not in rxn_ECs:
                continue
            # initialise reaction system object
            rs = rxn_syst.reaction(EC, 'ATLAS', ATLAS_ID)
            if os.path.isfile(output_dir+rs.pkl) is True and clean_system is False:
                count += 1
                continue
            if verbose:
                print('DB: ATLAS - EC:', EC, '-',
                      'DB ID:', ATLAS_ID, '-', count, 'of', no_rxns)
            # there are no ATLAS specific properties (for now)
            # could append pathways and orthologies from KEGG
            # using BridgIT results

            # get reaction system using DB specific function
            rxn_string = row['REACTION']
            rs = get_rxn_system(rs, rs.DB_ID, rxn_string)
            if rs.skip_rxn is False:
                # append compound information
                iterate_rs_components(rs, molecule_dataset=molecule_dataset)
            # pickle reaction system object to file
            # prefix (sRS for SABIO) + EC + EntryID .pkl
            rs.save_object(output_dir+rs.pkl)
            count += 1


def get_rxn_system(rs, ID, rxn_string):
    """Get reaction system from ATLAS reaction CSV.

    Uses KEGG API to convert Compound ID to molecule - online.

    Keywords:
        rs (class) - reaction system object
        ID (str) - DB reaction ID
        rxn_string (str) - reaction string from ATLAS CSV

    """
    # collate request output
    rs.components = []
    if '<=>' in rxn_string:
        # implies it is reversible
        reactants, products = rxn_string.split("<=>")
        print("reaction is reversible")
    else:
        reactants, products = rxn_string.split("=")

    reactants = [i.lstrip().rstrip() for i in reactants.split("+")]
    products = [i.lstrip().rstrip() for i in products.split("+")]

    # collect KEGG Compound/Glycan ID
    # check if reactant or product are compound or glycan
    # remove stoichiometry for now
    comp_list = []
    for r in reactants:
        if 'G' in r:
            # is glycan
            KID = 'G'+r.split('G')[1].rstrip()
        elif 'C' in r:
            # is compound
            KID = 'C'+r.split('C')[1].rstrip()
        comp_list.append((KID, 'reactant'))
    for r in products:
        if 'G' in r:
            # is glycan
            KID = 'G'+r.split('G')[1].rstrip()
        elif 'C' in r:
            # is compound
            KID = 'C'+r.split('C')[1].rstrip()
        comp_list.append((KID, 'product'))

    for comp in comp_list:
        # check for CHEBI ID in KEGG translations file
        translated = check_translator(comp[0])
        if translated is not None:
            pkl = translated
            print('collecting KEGG molecule using translator:', comp[0])
            new_mol = load_molecule(pkl, verbose=True)
            new_mol.KEGG_ID = comp[0]
            new_mol.translated = True
            rs.components.append(new_mol)
        else:
            chebiID = KEGGID_to_CHEBIID(KEGG_ID=comp[0])
            if chebiID is None:
                print('CHEBI ID not available - skipping whole reaction.')
                rs.skip_rxn = True
                rs.skip_reason = 'CHEBI ID of a component not available'
                return rs
            if chebiID is not None:
                new_mol = molecule(comp[0], comp[1], 'KEGG', chebiID)
                # add new_mol to reaction system class
                new_mol.KEGG_ID = comp[0]
                new_mol.translated = False
                new_mol.chebiID = chebiID
                rs.components.append(new_mol)
            else:
                new_mol = molecule(comp[0], comp[1], 'KEGG', comp[0])
                # add new_mol to reaction system class
                new_mol.KEGG_ID = comp[0]
                new_mol.translated = False
                new_mol.chebiID = None
                rs.components.append(new_mol)
    return rs


if __name__ == "__main__":
    print('testing and debugging here')
    import pandas as pd
    ATLAS_CSV = '/home/atarzia/psp/molecule_DBs/atlas/ATLAS-FULL.csv'
    ATLAS = pd.read_csv(ATLAS_CSV, chunksize=1000)
    for chunk in ATLAS:
        for idx, row in chunk.iterrows():
            print(idx)
            print(row)
            EC_list = row['REACTIONRULE'].split("|")
            print(EC_list)
            DB = 'ATLAS'
            DB_ID = row['ATLAS']
            print(DB, DB_ID)
            # add KEGG to DB list if present
            rxn = row['REACTION']
            print(rxn)
            break
        break
