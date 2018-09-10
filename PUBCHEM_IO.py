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
    if request.status_code != 404:
        SMILES = request.text.rstrip()
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
    if request.status_code != 404:
        iupac_name = request.text.rstrip()
    else:
        iupac_name = None

    return iupac_name
