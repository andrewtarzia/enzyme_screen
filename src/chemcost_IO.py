#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining use of chemcost module.

Author: Andrew Tarzia

Date Created: 04 Feb 2020

"""

from chemcost import PriceScraper


def is_purchasable(name, smiles):
    """
    Use chemcost module to determine if a molecule is purchasable.

    """

    p_scraper = PriceScraper(wait_time=5)
    is_purchasable = p_scraper.is_purchasable(smiles)
    print(f'{name} ({smiles}) purchasable?: {is_purchasable}')
    return is_purchasable
