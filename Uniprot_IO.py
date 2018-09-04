#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O with Uniprot DB.

Author: Andrew Tarzia

Date Created: 04 Sep 2018

"""

import requests


def get_sequence(UniprotID):
    """Get protein sequence from UniProt Fasta (using REST API).

    """
    # collect UniProtID data
    fasta_URL = 'https://www.uniprot.org/uniprot/'+UniprotID+'.fasta'
    request = requests.post(fasta_URL)
    request.raise_for_status()
    fasta_string = request.text.split('\n')
    sequence = ''.join(fasta_string[1:])
    return sequence
