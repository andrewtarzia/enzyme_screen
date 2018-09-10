#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module defining the molecule class.

Author: Andrew Tarzia

Date Created: 05 Sep 2018

"""
import cirpy


class molecule:
    """Class that defines the molecules extracted from a database.

    """

    def __init__(self, name, role, DB, DB_ID):
        self.name = name
        if role == 'Substrate':
            self.role = 'reactant'
        else:
            self.role = role.lower()
        self.DB = DB
        self.DB_ID = DB_ID
        self.InChi = None
        self.iupac_name = None

    def get_compound(self):
        """Get reaction system from SABIO reaction ID (rID).

        """
        if self.DB == 'SABIO':
            from SABIO_IO import get_cmpd_information
            # set DB specific properties
            self.cID = self.DB_ID
            get_cmpd_information(self)
        elif self.DB == 'KEGG':
            from CHEBI_IO import get_cmpd_information
            # set DB specific properties
            self.chebiID = self.DB_ID
            self.change_name = True  # name is set to KEGG C-ID at this point
            get_cmpd_information(self)
        elif self.DB == 'BKMS':
            from CHEBI_IO import get_cmpd_information
            # set DB specific properties
            self.chebiID = self.DB_ID
            get_cmpd_information(self)
        elif self.DB == 'BRENDA':
            from BRENDA_IO import get_cmpd_information
            # set DB specific properties
            get_cmpd_information(self)

    def cirpy_to_iupac(self):
        """Attempt to resolve IUPAC molecule name using CIRPY.

        Returns None if not possible.

        """
        self.iupac_name = cirpy.resolve(self.name, 'iupac_name')
