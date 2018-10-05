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


def get_SMILES_from_name(name):
    """Search for Canononical SMILES from PUBCHEM Compound using name.

    Tutorial: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial$_Toc458584413

    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/CanonicalSMILES/TXT'
    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        SMILES = request.text.rstrip().split('\n')[0]
    else:
        SMILES = None

    return SMILES


def get_IUPAC_from_name(name):
    """Search for IUPAC name from PUBCHEM Compound using name.

    Tutorial: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial$_Toc458584413

    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/IUPACName/TXT'

    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        iupac_name = request.text.rstrip().split('\n')[0]
    else:
        iupac_name = None

    return iupac_name


def get_logP_from_name(name):
    """Search for logP from PUBCHEM Compound using name.

    Computationally generated octanol-water partition coefficient or
    distribution coefficient. XLogP is used as a measure of hydrophilicity or
    hydrophobicity of a molecule.

    Tutorial: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial$_Toc458584413

    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/XLogP/TXT'
    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        try:
            XLogP = float(request.text.rstrip().split('\n')[0])
        except ValueError:
            XLogP = 'not found'
    else:
        XLogP = 'not found'

    return XLogP


def get_complexity_from_name(name):
    """Search for complexity from PUBCHEM Compound using name.

    The molecular complexity rating of a compound, computed using the
    Bertz/Hendrickson/Ihlenfeldt formula.

    Tutorial: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial$_Toc458584413

    """
    QUERY_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    QUERY_URL += name
    QUERY_URL += '/property/complexity/TXT'
    request = requests.post(QUERY_URL)
    # status 200 is success
    if request.status_code == 200:
        try:
            complexity = float(request.text.rstrip().split('\n')[0])
        except ValueError:
            complexity = 'not found'
    else:
        complexity = 'not found'

    return complexity


if __name__ == "__main__":
    name = 'phenylacetate'
    a = get_complexity_from_name(name)
    print(a)
