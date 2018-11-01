#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of CHEBI DB.

Using downloaded flat files. Offline.

Author: Andrew Tarzia

Date Created: 30 Aug 2018
"""

import pandas as pd
from libchebipy import ChebiEntity
from libchebipy import search as chebi_search
import DB_functions
from rdkit.Chem import AllChem as Chem
from re import sub
from json.decoder import JSONDecodeError
from molecule import charge_except


def search_for_compound_by_name(file, cmpd):
    """Search compounds.tsv file for matching name.

    """
    # read file line by line
    # header = 'ID', 'STATUS', 'CHEBI_ACCESSION',
    #          'SOURCE', 'PARENT_ID', 'NAME',
    #          'DEFINITION', 'MODIFIED_ON', 'CREATED_BY', 'STAR'
    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if name == cmpd or name.lower() == cmpd.lower():
                print('cmpd:', line_split)
                # check line for deprotonated carboxylic acid
                # if it is, get new molecule properties
                new_prop = check_line_for_carboxylate(line_split)
                print(new_prop)
                if new_prop is not None:
                    print('checking for new name')
                    for line2 in lines:
                        line_split2 = line2.split('\t')
                        ID, _, _, _, parent_id, name, _, _, _, star = line_split2
                        if ID == 'ID':
                            continue
                        if name == new_prop or name.lower() == new_prop.lower():
                            print(name, new_prop, ID, parent_id)
                            input('happy with change?')
                            return ID, parent_id, name, star
                    return None
                return ID, parent_id, name, star
    return None


def search_for_compound_by_id(file, chebiID):
    """Search compounds.tsv file for matching name.

    """
    # read file line by line
    # header = 'ID', 'STATUS', 'CHEBI_ACCESSION',
    #          'SOURCE', 'PARENT_ID', 'NAME',
    #          'DEFINITION', 'MODIFIED_ON', 'CREATED_BY', 'STAR'
    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if ID == str(chebiID):
                return ID, parent_id, name, star
    return None


def search_for_name_by_id(file, chebiID):
    """Search names.tsv file for matching CHEBI ID to get name.

    """
    # read file line by line
    # header = 'ID', 'COMPOUND_ID', 'TYPE',
    #          'SOURCE', 'NAME', 'ADAPTED', 'LANGUAGE'

    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, C_ID, _, _, name, _, _ = line_split
            if ID == 'ID':
                continue
            if C_ID == str(chebiID):
                return C_ID, name

    return None


def search_for_name_by_name(file, cmpd):
    """Search names.tsv file for matching name to get CHEBI ID.

    """
    # read file line by line
    # header = 'ID', 'COMPOUND_ID', 'TYPE',
    #          'SOURCE', 'NAME', 'ADAPTED', 'LANGUAGE'

    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, C_ID, _, _, name, _, _ = line_split
            if ID == 'ID':
                continue
            if name == cmpd or name.lower() == cmpd.lower():
                print('name:', line_split)
                return C_ID, name

    return None


def check_line_for_carboxylate(line):
    """Check a line from the CHEBI compounds file for carboxylic acid. Returns
    name of the protonated form.

    """
    # USING OFFLINE FILE -- NOT COMPLETE
    # ID, _, _, _, parent_id, name, notes, _, _, star = line
    # print(name)
    # print(notes)
    # if 'Conjugate base of ' in notes:
    #     if name[-3:] == 'ate':
    #         acid_name = notes.replace("Conjugate base of ", '')
    #         acid_name = acid_name.split('acid')[0] + 'acid'
    #         print('Conjugate base of:', acid_name, '<<<<<<<<<<')
    #         return acid_name
    # Using libchebipy -- Online:
    ID, _, _, _, parent_id, name, notes, _, _, star = line
    entity = ChebiEntity(ID)
    outgoings = entity.get_outgoings()
    for out in outgoings:
        type = out.get_type()
        id = out.get_target_chebi_id()
        if type == 'is_conjugate_base_of':
            print(out)
            # use regular expression to remove any non alphabetic characters
            # that may follow the name
            only_alph = sub('[^A-Za-z]', '', name)
            print(only_alph)
            if only_alph[-3:] == 'ate' and only_alph[-5:] != 'phate':
                acid_ID = id.replace("CHEBI:", "")
                print('acid ID:', acid_ID)
                acid_entity = ChebiEntity(acid_ID)
                acid_name = acid_entity.get_name()
                print('Conjugate base of:', acid_name, '<<<<<<<<<<')
                return acid_name
    return None


def check_entity_for_carboxylate(entity):
    """Check an entity from libchebipy for carboxylic acid. Returns
    name of the protonated form.

    The 'ate' term will also handle phosphates - so add check for that.

    """
    # Using libchebipy -- Online:
    name = entity.get_name()
    outgoings = entity.get_outgoings()
    for out in outgoings:
        type = out.get_type()
        id = out.get_target_chebi_id()
        if type == 'is_conjugate_base_of':
            # use regular expression to remove any non alphabetic characters
            # that may follow the name
            only_alph = sub('[^A-Za-z]', '', name)
            if only_alph[-3:] == 'ate' and only_alph[-5:] != 'phate':
                acid_ID = id.replace("CHEBI:", "")
                acid_entity = ChebiEntity(acid_ID)
                acid_name = acid_entity.get_name()
                # check if new SMILES is charged - don't change if it is
                acid_smiles = acid_entity.get_smiles()
                print('>>> new SMILES:', acid_smiles)
                if '-' in acid_smiles or '+' in acid_smiles:
                    if charge_except(acid_smiles) is False:
                        continue
                print(acid_name, acid_entity)
                print('---- Conjugate base of:', acid_name, '<<<<<<<<<<')
                return acid_name, acid_entity
    return None, None


def convert_nameID_to_parent(file, nameID):
    """Make sure ID extracted for the name is parent - get CHEBI ID of parent
    if not.

    """
    # read file line by line
    # header = 'ID', 'STATUS', 'CHEBI_ACCESSION',
    #          'SOURCE', 'PARENT_ID', 'NAME',
    #          'DEFINITION', 'MODIFIED_ON', 'CREATED_BY', 'STAR'
    # dont use pandas - code worked though:
    # cmpds_data = pd.read_table(compounds_file, delimiter='\t',
    #                            chunksize=1000)
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            # check line for carboxylic acid
            # if it is, get new molecule properties
            import sys
            sys.exit()
            new_prop = check_line_for_carboxylate(line_split)
            if new_prop is not None:
                # change properties
                print('need to change properties')
                import sys
                sys.exit()
            if nameID == ID:
                if parent_id != 'null':
                    desired_parent_id = parent_id
                    break
                else:
                    return ID, parent_id, name, star

    # now search for parent if necessary
    with open(file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            ID, _, _, _, parent_id, name, _, _, _, star = line_split
            if ID == 'ID':
                continue
            if ID == desired_parent_id:
                if parent_id == 'null':
                    return ID, parent_id, name, star

    return None


def get_structure(file, C_ID):
    """Get structure from CHEBI structures.csv file.

    """
    # use pandas (in chunks to avoid memory issues)
    # because it is a nice CSV file.
    structs = pd.read_csv(file, chunksize=1000)

    for chunk in structs:
        CID_frame = chunk[chunk['COMPOUND_ID'] == int(C_ID)]
        if len(CID_frame) > 0:
            # may be a 2D or 3D structure
            # output in hierarchy
            types = list(CID_frame['TYPE'])
            if 'SMILES' in types:
                row = CID_frame[CID_frame['TYPE'] == 'SMILES']
                structure = row['STRUCTURE'].iloc[0]
                s_type = row['TYPE'].iloc[0]
            elif 'InChI' in types:
                row = CID_frame[CID_frame['TYPE'] == 'InChI']
                structure = row['STRUCTURE'].iloc[0]
                s_type = row['TYPE'].iloc[0]
            elif 'mol' in types:
                row = CID_frame[CID_frame['TYPE'] == 'mol']
                structure = row['STRUCTURE'].iloc[0]
                s_type = row['TYPE'].iloc[0]
            else:
                structure = None
                s_type = None
            return structure, s_type

    return None, None


def get_chebiID_offline(mol_name):
    """Convert molecule name to chebiID using CHEBI DB files.

    Offline.
    Used by BKMS and BRENDA.

    Keywords:
        mol_name (str) - molecule name

    Returns:
        ID (str) - chebiID
    """
    DB_prop = DB_functions.get_DB_prop('CHEBI')
    compounds_file = DB_prop[0]+DB_prop[1]['cmpds_file']
    names_file = DB_prop[0]+DB_prop[1]['names_file']

    # search for name in compound file
    res = search_for_compound_by_name(compounds_file, mol_name)
    if res is None:
        # search for name in names file
        res = search_for_name_by_name(names_file, mol_name)
        if res is None:
            print('no match in DB')
            return None
        else:
            ID, name = res
            parent_id = None
    else:
        ID, parent_id, name, star = res

    # make sure is parent compound
    if parent_id != 'null':
        res = convert_nameID_to_parent(compounds_file, nameID=ID)
        if res is None:
            print("this should not happen - error with cross reference")
            print('check this!')
            import sys
            sys.exit()
        ID, parent_id, name, star = res
    print('chebiID:', ID)
    return ID


def clean_up_ID(ID):
    """Apply some clean up steps to Chebi IDs obtained from search.

    """
    print('>> CHEBI ID', ID)
    # check for carboxylate
    new_name, new_entity = check_entity_for_carboxylate(
                    entity=ChebiEntity(ID))
    if new_name is not None and new_entity is not None:
        ID = new_entity.get_id().replace("CHEBI:", '')
    # check for parent ID
    parent_ID = ChebiEntity(ID).get_parent_id()
    while parent_ID is not None:
        ID = parent_ID
        parent_ID = ChebiEntity(ID).get_parent_id()
    return ID


def find_synonym(results, target):
    """See if 'name' exists in results from a non-exact Chebi search.

    """
    for res in results:
        entity = ChebiEntity(res.get_id())
        names = entity.get_names()
        for name in names:
            if name.get_name().lower() == target.lower():
                print('found synonym - give ID')
                return res.get_id().replace("CHEBI:", '')
    return None


def get_chebiID(mol_name, iupac_name=False):
    """Get ChebiID using libchebipy API.

    Online.
    Used by BKMS and BRENDA.

    Keywords:
        mol_name (str) - molecule name

    Returns:
        ID (str) - chebiID
    """
    # set an order of search conditions
    # (identifier, exact search)
    search_conditions = [(str(mol_name), True),
                         (str(mol_name).lower(), True),
                         (str(mol_name), False),
                         (str(mol_name).lower(), False),
                         (str(iupac_name), True),
                         (str(iupac_name).lower(), True),
                         (str(iupac_name), False),
                         (str(iupac_name).lower(), False)]
    for search in search_conditions:
        if search[0] != str(False) and search[0] != str(False).lower():
            try:
                search_result = chebi_search(search[0], search[1])
            except (JSONDecodeError, KeyError):
                print('failed search due to internet(?) - except')
                try:
                    search_result = chebi_search(search[0], search[1])
                except (JSONDecodeError, KeyError):
                    print('failed search due to internet(?) - except')
                    continue
            except ValueError:
                print('error with CHEBI DB file - skip.')
                continue
            if len(search_result) == 1:
                ID = search_result[0].get_id().replace("CHEBI:", '')
                ID = clean_up_ID(ID=ID)
                return ID
            elif len(search_result) > 1:
                print('multiple matches to exact search', search)
                # search through synonyms for our target name
                ID = find_synonym(results=search_result, target=mol_name)
                return ID
            else:
                print('no match in DB for search:', search)
    return None


def get_cmpd_information_offline(molec):
    """Get information from CHEBI Database of a compound from CHEBI ID.

    Done Offline unless necessary.
    molec must have attribute 'chebiID' as integer.

    """

    DB_prop = DB_functions.get_DB_prop('CHEBI')
    compounds_file = DB_prop[0]+DB_prop[1]['cmpds_file']
    names_file = DB_prop[0]+DB_prop[1]['names_file']
    structures_file = DB_prop[0]+DB_prop[1]['strct_file']

    # set name by searching compound file
    res = search_for_compound_by_id(compounds_file, molec.chebiID)
    if res is None:
        print('chebiID not found:', molec.chebiID)
        print('no match in DB - this should not happen for CHEBI ID search')
        print('check this!')
        print('Exitting....')
        import sys
        sys.exit()
    else:
        ID, parent_id, name, star = res
        molec.name = name
        molec.change_name = False

    # make sure is parent compound
    if parent_id != 'null':
        res = convert_nameID_to_parent(compounds_file, nameID=ID)
        if res is None:
            print("this should not happen - error with cross reference")
            print('check this!')
            print('Exitting....')
            import sys
            sys.exit()
        ID, parent_id, name, star = res
        molec.name = name
        molec.change_name = False
        molec.chebiID = int(ID)

    # get structure using CHEBI ID
    # structures.csv - read in, get COMPOUND ID match then extract the
    # get SMILES
    structure, s_type = get_structure(structures_file, molec.chebiID)
    print(structure, s_type)
    if structure is not None:
        # is structure a MolBlock or Smiles
        if s_type == 'mol':
            # convert structure to SMILEs
            rdkitmol = Chem.MolFromMolBlock(structure)
            if rdkitmol is None:
                print('structure could not be deciphered')
                smile = None
                molec.SMILES = smile
                molec.mol = None
                print('probably a polymeric structure - skipping.')
            else:
                rdkitmol.Compute2DCoords()
                smile = Chem.MolToSmiles(rdkitmol)
                molec.SMILES = smile
                # remove molecules with generalised atoms
                if '*' in smile:
                    molec.mol = None
                else:
                    molec.mol = rdkitmol
        elif s_type == 'SMILES':
            smile = structure
            rdkitmol = Chem.MolFromSmiles(smile)
            if rdkitmol is None:
                print('structure could not be deciphered')
                molec.SMILES = smile
                molec.mol = None
            else:
                rdkitmol.Compute2DCoords()
                molec.SMILES = smile
                # remove molecules with generalised atoms
                if '*' in smile:
                    molec.mol = None
                else:
                    molec.mol = rdkitmol
        elif s_type == 'InChI':
            rdkitmol = Chem.MolFromInchi(structure)
            rdkitmol.Compute2DCoords()
            smile = Chem.MolToSmiles(rdkitmol)
            molec.SMILES = smile
            # remove molecules with generalised atoms
            if '*' in smile:
                molec.mol = None
            else:
                molec.mol = rdkitmol
        elif s_type == 'InChIKey':
            rdkitmol = Chem.MolFromInchi(structure)
            rdkitmol.Compute2DCoords()
            smile = None
            molec.SMILES = smile
            molec.mol = None
            print('molecule given as InChIKey - ambiguous')
            print('probably a generic structure - skipping.')
    else:
        # try using the CHEBI API
        # libChEBIpy (https://github.com/libChEBI/libChEBIpy)
        print('testing libchebipy...')
        entity = ChebiEntity(molec.chebiID)
        smile = entity.get_smiles()
        print('libchebipy result:', smile)
        if smile is not None:
            rdkitmol = Chem.MolFromSmiles(smile)
            if rdkitmol is None:
                print('structure could not be deciphered')
                molec.SMILES = smile
                molec.mol = None
            else:
                rdkitmol.Compute2DCoords()
                molec.SMILES = smile
                # remove molecules with generalised atoms
                if '*' in smile:
                    molec.mol = None
                else:
                    molec.mol = rdkitmol
        elif smile is None:
            molec.SMILES = smile
            molec.mol = None
            print('molecule does not have recorded structure in CHEBI DB')
            print('probably a generic structure - skipping.')
        # save InChiKey
        iKEY = entity.get_inchi_key()
        if iKEY is not None:
            molec.InChiKey = iKEY


def get_cmpd_information(molec):
    """Get information from CHEBI Database of a compound from CHEBI ID.

    Online using libChEBIpy (https://github.com/libChEBI/libChEBIpy)

    """
    if molec.chebiID is None and molec.iupac_name is not None:
        # try one more time for chebi ID
        chebiID = get_chebiID(mol_name=molec.name, iupac_name=molec.iupac_name)
        if chebiID is None:
            print('cannot get structure from chebi')
            return None
        molec.chebiID = chebiID
    # get entity with chebiID
    entity = ChebiEntity(molec.chebiID)
    # check for parent ID
    ID = molec.chebiID
    parent_ID = entity.get_parent_id()
    while parent_ID is not None:
        ID = parent_ID.replace("CHEBI:", '')
        parent_ID = ChebiEntity(ID).get_parent_id()
    # change chebiID to parent ID
    if ID != molec.chebiID:
        molec.chebiID = ID
        entity = ChebiEntity(molec.chebiID)
    # set name if name is only a code at this point
    try:
        if molec.change_name is True:
            molec.name = entity.get_name()
            molec.change_name = False
    except AttributeError:
        molec.change_name = False
    # check for deprotonated carboxylate
    new_name, new_entity = check_entity_for_carboxylate(
            entity=entity)
    if new_name is not None and new_entity is not None:
        entity = new_entity
        molec.chebiID = entity.get_id().replace("CHEBI:", "")

    # get structure
    # SMILES
    smile = entity.get_smiles()
    print('libchebipy result:', smile)
    if smile is not None:
        rdkitmol = Chem.MolFromSmiles(smile)
        if rdkitmol is None:
            print('structure could not be deciphered')
            molec.SMILES = smile
            molec.mol = None
        else:
            rdkitmol.Compute2DCoords()
            molec.SMILES = smile
            # remove molecules with generalised atoms
            if '*' in smile:
                molec.mol = None
            else:
                molec.mol = rdkitmol
    elif smile is None:
        molec.SMILES = smile
        molec.mol = None
        print('molecule does not have recorded structure in CHEBI DB')
        print('probably a generic structure - skipping.')
    # save InChiKey
    iKEY = entity.get_inchi_key()
    if iKEY is not None:
        molec.InChiKey = iKEY
    # save inchi
    inchi = entity.get_inchi()
    if inchi is not None:
        molec.InChi = inchi
