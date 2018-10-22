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


def run_request(query):
    """Query URL and handle errors.

    """
    request = requests.post(query)
    # status 200 is success
    if request.status_code == 200:
        if '\n' in request.text.rstrip():
            print(request.text.rstrip())
            raise ValueError('new line found')
        else:
            return request.text.rstrip().split('\n')[0]
    else:
        return None


def hier_name_search(molecule, property):
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
            if molecule.InChiKey is not None:
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
        if molecule.iupac_name is not None:
            QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
            QUERY_URL_fin = QUERY_URL + molecule.iupac_name
            QUERY_URL_fin += '/property/'+property+'/TXT'
            result = run_request(query=QUERY_URL_fin)
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
            result = run_request(query=QUERY_URL_fin)
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


if __name__ == "__main__":
    name = 'guanidine'
    a = get_complexity_from_name(name)
    print(a)
