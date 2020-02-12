#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining use of chemcost module.

Author: Andrew Tarzia

Date Created: 04 Feb 2020

"""

from chemcost import PriceScraper
import requests


def is_purchasable(name, smiles):
    """
    Use chemcost module to determine if a molecule is purchasable.

    """
    p_scraper = PriceScraper()
    zinc_id = p_scraper.get_zinc_id(smiles)
    if len(zinc_id) > 1 and not isinstance(zinc_id, str):
        print(f'{len(zinc_id)} zinc IDs found.')
        # Check KEGG CAS numbers.
        KEGG_cas_numbers = get_cas_number(name)
        found = False
        for zid, ztitle in zinc_id:
            print(f'checking {zid} for {KEGG_cas_numbers}')
            cas_nums = p_scraper.get_cas_numbers(zinc_id=zid)
            if cas_nums is None:
                continue
            for KEGG_cas in KEGG_cas_numbers:
                if KEGG_cas in cas_nums:
                    zinc_id = zid
                    print(f'found matching {zinc_id}!')
                    found = True
                    break
            if found:
                break

    input('check')
    # If zinc_id never got updated, then cas search failed.
    if isinstance(zinc_id, list):
        print(f'no zinc ID found!')
        zinc_id = None

    if zinc_id is None:
        return False
    else:
        print(f'ZINC ID: {zinc_id}')
        return p_scraper.is_purchasable(
            zinc_id
        )


def get_cas_number(KEGG_ID):
    """
    Use KEGG API to collect CAS number using KEGG API.

    Examples: https://www.kegg.jp/kegg/rest/keggapi.html#get

    """
    # get compound information from KEGG API
    URL = 'http://rest.kegg.jp/get/'+KEGG_ID+'/'
    request = requests.post(URL)
    # KEGG API will give a 404 if MOL file cannot be collected
    if request.status_code == 404:
        print(KEGG_ID, 'has no mol file associated with it')
        return None
    elif request.status_code == 200:
        output = request.text
        print('convert', KEGG_ID, 'to CAS Number...')
        result_line = [i for i in output.split('\n') if 'CAS:' in i]
        if len(result_line) == 0:
            return None
        res_line = result_line[0].split('CAS:')[1].lstrip().rstrip()
        result = [
            i for i in res_line.split(' ')
        ]
        return result
    elif request.status_code == 400:
        # implies a bad request/syntax
        print(KEGG_ID, 'has bad syntax')
        return None
    else:
        print("haven't come across this yet.")
        raise ValueError(
            f'KEGG ID: {KEGG_ID} gives this error: '
            f'{request.status_code}'
        )
