#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O with PUBCHEM DATABASE API (PUG REST).

Online.

Author: Andrew Tarzia

Date Created: 10 Sep 2018


"""

import requests
import sys
import pubchempy as pcp


def run_request(query, smiles=False, option=False):
    """Query URL and handle errors. New line errors can be handled if recieving
    SMILES.

    """
    request = requests.post(query)
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            if option is not False:
                print('picking option:', option)
                return request.text.rstrip().split('\n')[option]
            elif smiles is False:
                print(request.text.rstrip())
                print('new line found!')
                raise ValueError('new line found')
            elif smiles is True:
                return request.text.rstrip(), True
        else:
            return request.text.rstrip().split('\n')[0]
    else:
        return None


def hier_name_search(molecule, property, option=False):
    """Search for molecule property in PUBCHEM using a hierarchy of name spaces

    Order:
        1 - pubchem ID
        2 - KEGG ID
        3 - chebiID
        4 - chebiID to InChIKey
        5 - IUPAC name
        6 - name

    Properties:
        CanononicalSMILES
        IUPACName
        XLogP
        complexity
        InChiKey

    Tutorial: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial$_Toc458584413
    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'

    try:
        if molecule.PubChemID is not None:
            QUERY_URL_fin = QUERY_URL + molecule.PubChemID
            QUERY_URL_fin += '/property/'+property+'/TXT'
            result = run_request(query=QUERY_URL_fin)
            if result is not None:
                print('passed pubchemID')
                return result
    except (AttributeError, ValueError):
        print('failed pubchemID')
        pass
    try:
        if molecule.KEGG_ID is not None:
            QUERY_URL_fin = QUERY_URL + molecule.KEGG_ID
            QUERY_URL_fin += '/property/'+property+'/TXT'
            result = run_request(query=QUERY_URL_fin)
            if result is not None:
                print('passed KEGG ID')
                return result
    except (AttributeError, ValueError):
        print('failed KEGG ID')
        pass
    try:
        if molecule.chebiID is not None:
            QUERY_URL_fin = QUERY_URL + 'chebi:'+molecule.chebiID
            QUERY_URL_fin += '/property/'+property+'/TXT'
            result = run_request(query=QUERY_URL_fin)
            if result is not None:
                print('passed chebiID')
                return result
    except (AttributeError, ValueError):
        print('failed chebiID')
        pass
    try:
        if molecule.chebiID is not None:
            if molecule.InChiKey is None:
                # try using the CHEBI API
                # libChEBIpy (https://github.com/libChEBI/libChEBIpy)
                print('using libchebipy to get InChiKey...')
                from libchebipy import ChebiEntity
                entity = ChebiEntity(molecule.chebiID)
                iKEY = entity.get_inchi_key()
                print(iKEY)
            else:
                iKEY = molecule.InChiKey
            QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
            QUERY_URL_fin = QUERY_URL + iKEY
            QUERY_URL_fin += '/property/'+property+'/TXT'
            result = run_request(query=QUERY_URL_fin)
            if result is not None:
                print('passed chebiID/inchiKey')
                return result
    except (AttributeError, ValueError):
        print('failed chebiID/inchiKey')
        pass
    try:
        if molecule.InChiKey is not None:
            iKEY = molecule.InChiKey
            QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/'
            QUERY_URL_fin = QUERY_URL + iKEY
            QUERY_URL_fin += '/property/'+property+'/TXT'
            result = run_request(query=QUERY_URL_fin)
            if result is not None:
                print('passed inchiKey')
                return result
    except (AttributeError, ValueError):
        print('failed inchiKey')
        pass
    try:
        if molecule.iupac_name is not None:
            QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
            QUERY_URL_fin = QUERY_URL + molecule.iupac_name
            QUERY_URL_fin += '/property/'+property+'/TXT'
            if property == 'CanonicalSMILES':
                result = run_request(query=QUERY_URL_fin, smiles=True)
                print('res', result)
                if type(result) == tuple:
                    # handle new line errors in SMILES
                    text, boolean = result
                    if boolean is True:
                        # pick the uncharged SMILES
                        for option, smi in enumerate(text.split('\n')):
                            print('smiles1:', smi)
                            if '-' in smi or '+' in smi:
                                # charged
                                continue
                            return smi, option
                elif type(result) == str and result is not None:
                    print('passed name')
                    return result
            else:
                result = run_request(query=QUERY_URL_fin, option=option)
            if result is not None:
                print('passed IUPAC name')
                return result
    except (AttributeError, ValueError):
        print('failed IUPAC name')
        pass
    try:
        if molecule.name is not None:
            QUERY_URL_fin = QUERY_URL + molecule.name
            QUERY_URL_fin += '/property/'+property+'/TXT'
            if property == 'CanonicalSMILES':
                result = run_request(query=QUERY_URL_fin, smiles=True)
                print('res', result)
                if type(result) == tuple:
                    # handle new line errors in SMILES
                    text, boolean = result
                    if boolean is True:
                        # pick the uncharged SMILES
                        for option, smi in enumerate(text.split('\n')):
                            print('smiles1:', smi)
                            if '-' in smi or '+' in smi:
                                # charged
                                continue
                            return smi, option
                elif type(result) == str and result is not None:
                    print('passed name')
                    return result
            else:
                result = run_request(query=QUERY_URL_fin, option=option)
            if result is not None:
                print('passed name')
                return result
    except (AttributeError, ValueError):
        print('failed name')
        import sys
        sys.exit()

    return None


def get_SMILES_from_name(name):
    """Search for Canononical SMILES from PUBCHEM Compound using name.

    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/CanonicalSMILES/TXT'
    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            sys.exit('has new line - exitting...')
        SMILES = request.text.rstrip().split('\n')[0]
    else:
        SMILES = None

    return SMILES


def get_IUPAC_from_name(name):
    """Search for IUPAC name from PUBCHEM Compound using name.
    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/IUPACName/TXT'

    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            sys.exit('has new line - exitting...')
        iupac_name = request.text.rstrip().split('\n')[0]
    else:
        iupac_name = None

    return iupac_name


def get_logP_from_name(name):
    """Search for logP from PUBCHEM Compound using name.

    Computationally generated octanol-water partition coefficient or
    distribution coefficient. XLogP is used as a measure of hydrophilicity or
    hydrophobicity of a molecule.

    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/XLogP/TXT'
    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            sys.exit('has new line - exitting...')
        try:
            items = request.text.rstrip().split('\n')
            print(items)
            if len(items) > 1:
                # all equal?
                se = list(set(items))
                print(se)
                if len(se) == 1:
                    XLogP = float(items[0])
                else:
                    XLogP = items
            elif len(items) > 0:
                XLogP = float(items[0])
        except ValueError:
            XLogP = 'not found'
    else:
        XLogP = 'not found'

    return XLogP


def get_complexity_from_name(name):
    """Search for complexity from PUBCHEM Compound using name.

    The molecular complexity rating of a compound, computed using the
    Bertz/Hendrickson/Ihlenfeldt formula.
    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/complexity/TXT'
    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            sys.exit('has new line - exitting...')
        try:
            items = request.text.rstrip().split('\n')
            print(items)
            if len(items) > 1:
                # all equal?
                se = list(set(items))
                print(se)
                if len(se) == 1:
                    complexity = float(items[0])
                else:
                    complexity = items
            elif len(items) > 0:
                complexity = float(items[0])
        except ValueError:
            complexity = 'not found'
    else:
        complexity = 'not found'

    return complexity


def get_logP_from_SMILES(SMILES):
    """Search for logP from PUBCHEM Compound using SMILES.

    Computationally generated octanol-water partition coefficient or
    distribution coefficient. XLogP is used as a measure of hydrophilicity or
    hydrophobicity of a molecule.
    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'
    QUERY_URL += SMILES
    QUERY_URL += '/property/XLogP/TXT'
    request = requests.post(QUERY_URL)
    print(SMILES)
    print(QUERY_URL)
    print(request.text.rstrip())
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            sys.exit('has new line - exitting...')
        try:
            items = request.text.rstrip().split('\n')
            print(items)
            if len(items) > 1:
                # all equal?
                se = list(set(items))
                print(se)
                if len(se) == 1:
                    XLogP = float(items[0])
                else:
                    XLogP = items
            elif len(items) > 0:
                XLogP = float(items[0])
        except ValueError:
            XLogP = 'not found'
    else:
        XLogP = 'not found'

    return XLogP


def get_complexity_from_SMILES(SMILES):
    """Search for complexity from PUBCHEM Compound using SMILES.

    The molecular complexity rating of a compound, computed using the
    Bertz/Hendrickson/Ihlenfeldt formula.
    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'
    QUERY_URL += SMILES
    QUERY_URL += '/property/complexity/TXT'
    request = requests.post(QUERY_URL)
    print(SMILES)
    print(QUERY_URL)
    print(request.text.rstrip())
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            sys.exit('has new line - exitting...')
        try:
            items = request.text.rstrip().split('\n')
            print(items)
            if len(items) > 1:
                # all equal?
                se = list(set(items))
                print(se)
                if len(se) == 1:
                    complexity = float(items[0])
                else:
                    complexity = items
            elif len(items) > 0:
                complexity = float(items[0])
        except ValueError:
            complexity = 'not found'
    else:
        complexity = 'not found'

    return complexity


def run_request_pcp(ident, namespace,
                    smiles=False, option=False):
    """Query PubChem PUG API and handle errors.

    """
    result = pcp.get_compounds(identifier=ident, namespace=namespace)
    if len(result) > 1:
        print('multiple hits', smiles)
        if option is not False:
            print('picking option:', option)
            return result[option]
        elif smiles is False:
            print(result)
            print('new line found!')
            raise ValueError('new line found')
        elif smiles is True:
            return result, True
    elif len(result) == 1:
        print('one hit')
        return result[0]
    else:
        return None


def extract_property(property, result):
    """Extract the desired property from the pubchempy result.

    """
    if property == 'CanonicalSMILES':
        return result.canonical_smiles
    if property == 'IUPACName':
        return result.iupac_name
    if property == 'XLogP':
        return result.xlogp
    if property == 'complexity':
        return result.complexity
    if property == 'InChiKey':
        return result.inchikey
    if property == 'PubChemID':
        return result.cid
    if property == 'synonyms':
        return result.synonyms


def hier_name_search_pcp(molecule, property, option=False):
    """Search for molecule property in PUBCHEM using a hierarchy of name spaces
    using pubchempy.

    Property can now be a list.

    Order:
        1 - pubchem ID
        2 - KEGG ID
        3 - chebiID
        4 - chebiID to InChIKey
        5 - IUPAC name
        6 - name

    Properties:
        CanononicalSMILES
        IUPACName
        XLogP
        complexity
        InChiKey
        PubChemID
        synonyms

    """
    if type(property) is not list:
        property = [property]
    try:
        if molecule.PubChemID is not None:
            result = run_request_pcp(ident=molecule.PubChemID,
                                     namespace='cid')
            if result is not None:
                print('passed pubchemID')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed pubchemID')
        pass
    try:
        if molecule.KEGG_ID is not None:
            result = run_request_pcp(ident=molecule.KEGG_ID,
                                     namespace='name')
            if result is not None:
                print('passed KEGG ID')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed KEGG ID')
        pass
    try:
        if molecule.chebiID is not None:
            result = run_request_pcp(ident='chebi:'+molecule.chebiID,
                                     namespace='name')
            if result is not None:
                print('passed chebiID')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed chebiID')
        pass
    try:
        if molecule.chebiID is not None:
            if molecule.InChiKey is None:
                # try using the CHEBI API
                # libChEBIpy (https://github.com/libChEBI/libChEBIpy)
                print('using libchebipy to get InChiKey...')
                from libchebipy import ChebiEntity
                entity = ChebiEntity(molecule.chebiID)
                iKEY = entity.get_inchi_key()
                print(iKEY)
            else:
                iKEY = molecule.InChiKey
            result = run_request_pcp(ident=iKEY,
                                     namespace='inchikey')
            if result is not None:
                print('passed chebiID/inchiKey')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed chebiID/inchiKey')
        pass
    try:
        if molecule.InChiKey is not None:
            result = run_request_pcp(ident=molecule.InChiKey,
                                     namespace='inchikey')
            if result is not None:
                print('passed inchiKey')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed inchiKey')
        pass
    try:
        if molecule.iupac_name is not None:
            if property == 'CanonicalSMILES':
                result = run_request_pcp(ident=molecule.iupac_name,
                                         namespace='name',
                                         smiles=True)
                print('res', result, type(result))
                if type(result) == tuple:
                    # handle new line errors in SMILES
                    text, boolean = result
                    if boolean is True:
                        # pick the uncharged SMILES
                        for option, smi in enumerate(text.split('\n')):
                            print('smiles1:', smi)
                            if '-' in smi or '+' in smi:
                                # charged
                                continue
                            return smi, option
                elif type(result) == str and result is not None:
                    print('passed IUPAC name')
                    print('I am interested in what this result is:')
                    print(result)
                    import sys
                    sys.exit()
                    if len(property) > 1:
                        return [extract_property(i, result) for i in property]
                    else:
                        return [extract_property(i, result) for i in property][0]
            else:
                result = run_request_pcp(ident=molecule.iupac_name,
                                         namespace='name',
                                         option=option)
            if result is not None:
                print('passed IUPAC name')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed IUPAC name')
        pass
    try:
        print('trying name!', property)
        if molecule.name is not None:
            if property == 'CanonicalSMILES':
                result = run_request_pcp(ident=molecule.name,
                                         namespace='name',
                                         smiles=True)
                print('res', result, type(result))
                if type(result) == tuple:
                    # handle new line errors in SMILES
                    text, boolean = result
                    if boolean is True:
                        # pick the uncharged SMILES
                        for option, Compound in enumerate(text):
                            synon = [i.lower() for i in Compound.synonyms]
                            if molecule.name.lower() in synon:
                                # ignore charged species
                                smi = Compound.canonical_smiles
                                if '-' in smi or '+' in smi:
                                    continue
                                print('smiles1:', smi)
                                return smi, option
                elif type(result) == str and result is not None:
                    print('passed name')
                    print('I am interested in what this result is:')
                    print(result)
                    import sys
                    sys.exit()
                    return result
            else:
                result = run_request_pcp(ident=molecule.name,
                                         namespace='name',
                                         option=option)
            if result is not None:
                print('passed name')
                if len(property) > 1:
                    return [extract_property(i, result) for i in property]
                else:
                    return [extract_property(i, result) for i in property][0]
    except (AttributeError, ValueError):
        print('failed name')
    return None


def pubchem_synonym(mol):
    """Search PUBCHEM for molecule names and check resultant synonyms for CHEBI
    ID.

    """
    print('search pubchem by name for compound with synonym chebiID')
    for Compound in pcp.get_compounds(mol.name, 'name'):
        synon = [i.lower() for i in Compound.synonyms]
        if mol.name.lower() in synon:
            # ignore charged species
            smi = Compound.canonical_smiles
            if '-' in smi or '+' in smi:
                continue
            for syn in Compound.synonyms:
                if 'CHEBI:' in syn:
                    chebiID = syn.replace("CHEBI:", '')
                    return chebiID
    return None


def pubchem_check_smiles(mol):
    """Search PUBCHEM for molecule names and return SMILES for molecules with
    no ChebiID.

    """
    print('collecting SMILES from PUBCHEM in BKMS with Chebi == None')
    smiles_search = hier_name_search_pcp(mol,
                                         'CanonicalSMILES')
    print('search result', smiles_search)
    if smiles_search is not None:
        if len(smiles_search) == 2:
            smiles = smiles_search[0]
            option = smiles_search[1]
        elif len(smiles_search) > 0:
            smiles = smiles_search
            option = 0
        else:
            smiles = None
    else:
        smiles = None
    if smiles is not None:
        # implies we got the SMILES from a PUBCHEM search
        # collect other properties from PUBCHEM using the option if
        # a new line was found
        mol.SMILES = smiles
        mol.InChiKey = hier_name_search_pcp(mol,
                                            'InChiKey',
                                            option=option)
        print('IKEY:', mol.InChiKey)
        mol.iupac_name = hier_name_search_pcp(mol,
                                              'IUPACName',
                                              option=option)
        print('iupac_name', mol.iupac_name)
        return mol, smiles
    else:
        return mol, None


if __name__ == "__main__":
    name = 'guanidine'
    a = get_complexity_from_name(name)
    print(a)
