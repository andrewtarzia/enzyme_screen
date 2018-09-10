#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of BRENDA DB.

Author: Andrew Tarzia

Date Created: 24 Apr 2018

"""

# ensure cpickle usage
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import molecule
import rxn_syst
import os


class Reaction_system:
    """Class that defines a singular reaction system extracted from BRENDA.

    """

    def __init__(self, EC_no, prop, PR, count):
        self.target_PR = PR
        self.react_mol = []
        self.prod_mol = []
        self.activating_mol = []
        self.cofactor_mol = []
        self.meta = ''
        self.reaction_cat = ''
        self.reaction_type = ''
        self.assoc_PR = []
        self.assoc_refs = []
        # number of sub units is None if unknown
        self.no_sub_units = None
        # flags for reaction system
        # True if confirmed by BRENDA
        # None if data not available in BRENDA
        self.reversible = None
        self.PTM = None  # post translation mods?
        self.cofactors = None  # requires cofactors?
        self.pickle_name = 'RS-'+EC_no+'-'+prop+'_'+PR+'_'+count+'.pkl'

    def save_object(self, filename):
        """Pickle reaction system object to file.

        """
        # Overwrites any existing file.
        with open(filename, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def load_object(filename):
        """unPickle reaction system object from file.

        """
        with open(filename, 'rb') as input:
            self = pickle.load(input)
            return self

    def print_rxn_system(self):
        """Fancy print of reaction system.

        """
        print('--------------------------')
        print('Reaction system in:', self.pickle_name)
        print('Reaction Catalysed:')
        print(self.reaction_cat)
        print('Reaction Type:')
        print(self.reaction_type)
        print('--------------------------')
        string = ''
        for i in self.react_mol:
            if i != self.react_mol[-1]:
                string += i + ' + '
            else:
                string += i
        print('Reactants:', string, '-------->')
        string = ''
        for i in self.prod_mol:
            if i != self.prod_mol[-1]:
                string += i + ' + '
            else:
                string += i
        print('Products:', string)
        string = ''
        for i in self.activating_mol:
            if i != self.activating_mol[-1]:
                string += i + ' OR '
            else:
                string += i
        print('Activating Molecules:', string)
        print('--------------------------')
        print('No. Sub units =', self.no_sub_units)
        print('Reversible?:', self.reversible)
        print('Co factors?:', self.cofactors)
        string = ''
        for i in self.cofactor_mol:
            if i != self.cofactor_mol[-1]:
                string += i + ' OR '
            else:
                string += i
        print('Cofactor Molecules:', string)
        print('Post Translational Mods?:', self.PTM)
        print('--------------------------')
        print('Meta (in full):')
        print(self.meta)
        print('References:', self.assoc_refs)
        print('--------------------------')

    def extract_subunit_info(self, br_data, PR):
        """Extract all sub unit information for the PR of this reaction system.
            It is possible that there is no information, or multiple entries.

        """
        nprop = 'SU'
        output_nprop = []
        list_of_nprop = br_data[nprop]
        EC_nprop_PR = get_prop_PR_codes(list_of_nprop)
        for i, entry in enumerate(list_of_nprop):
            if PR in EC_nprop_PR[i]:
                # collect reported value
                # by splitting the string using the
                # "#" at the end of the PR codes and the
                # "(" at the start of the notes
                value = entry.split("#")[2].split("(")[0].lstrip().rstrip()
                output_nprop.append(value)

        self.no_sub_units = output_nprop

    def extract_PTM(self, br_data, PR):
        """Extract all PTM information for the PR of this reaction system.
            It is possible that there is no information, or multiple entries.

        """
        nprop = 'PM'
        output_nprop = []
        list_of_nprop = br_data[nprop]
        EC_nprop_PR = get_prop_PR_codes(list_of_nprop)
        for i, entry in enumerate(list_of_nprop):
            if PR in EC_nprop_PR[i]:
                # collect reported value
                # by splitting the string using the
                # "#" at the end of the PR codes and the
                # "(" at the start of the notes
                value = entry.split("#")[2].split("(")[0].lstrip().rstrip()
                output_nprop.append(value)

        self.PTM = output_nprop

    def extract_cofactor_info(self, br_data, PR):
        """Extract all cofactor information for the PR of this reaction system.
            It is possible that there is no information, or multiple entries.

        """
        nprop = 'CF'
        output_nprop = []
        list_of_nprop = br_data[nprop]
        EC_nprop_PR = get_prop_PR_codes(list_of_nprop)
        for i, entry in enumerate(list_of_nprop):
            if PR in EC_nprop_PR[i]:
                # collect reported value
                # by splitting the string using the
                # "#" at the end of the PR codes and the
                # "(" at the start of the notes
                value = entry.split("#")[2].split("(")[0].lstrip().rstrip()
                output_nprop.append(value)

        self.cofactor_mol = output_nprop
        self.cofactors = True

    def extract_activating_mol(self, br_data, PR):
        """Extract all activating molecules for the PR of this reaction system.
            It is possible that there is no information, or multiple entries.

        """
        nprop = 'AC'
        output_nprop = []
        list_of_nprop = br_data[nprop]
        EC_nprop_PR = get_prop_PR_codes(list_of_nprop)
        for i, entry in enumerate(list_of_nprop):
            if PR in EC_nprop_PR[i]:
                # collect reported value
                # by splitting the string using the
                # "#" at the end of the PR codes and the
                # "(" at the start of the notes
                value = entry.split("#")[2].split("(")[0].lstrip().rstrip()
                output_nprop.append(value)

        self.activating_mol = output_nprop

    def extract_general_rxn_info(self, br_data):
        """Extract all activating molecules for the PR of this reaction system.
            It is possible that there is no information, or multiple entries.

        """
        nprop = 'RT'
        output_nprop = []
        list_of_nprop = br_data[nprop]
        for i, entry in enumerate(list_of_nprop):
            value = entry.rstrip()
            output_nprop.append(value)

        self.reaction_type = output_nprop
        nprop = 'RE'
        output_nprop = []
        list_of_nprop = br_data[nprop]
        for i, entry in enumerate(list_of_nprop):
            value = entry.rstrip()
            output_nprop.append(value)

        self.reaction_cat = output_nprop


def initialise_br_dict():
    """
    Initialise the brenda data dictionary and the assoc. symbol definitions.

    Returns:
        br_symbols (dict) - dictionary of brenda symbols
        br_data (dict) - dictionary of brenda data
    """
    # get BRENDA symbols:
    # initials: ("HEADER", "shortname")
    br_symbols = {
        "ID": ("", "ec_no"),
        "PR": ("PROTEIN", 'protein'),
        "RN": ('RECOMMENDED_NAME', 'reco_name'),
        "SN": ('SYSTEMATIC_NAME', 'syst_name'),
        "SY": ('SYNONYMS', 'syn'),
        "RE": ('REACTION', 'rxn'),
        "RT": ('REACTION_TYPE', 'rxn_type'),
        "ST": ('SOURCE_TISSUE', 'src_tis'),
        "NSP": ('NATURAL_SUBSTRATE_PRODUCT', 'nat_subs_prod'),
        "SP": ('SUBSTRATE_PRODUCT', 'subs_prod'),
        "TN": ('TURNOVER_NUMBER', 'turnover'),
        "KM": ('KM_VALUE', 'km'),
        "PHO": ('PH_OPTIMUM', 'ph_opt'),
        "SA": ('SPECIFIC_ACTIVITY', 'spec_act'),
        "TO": ('TEMPERATURE_OPTIMUM', 'T_opt'),
        "CF": ('COFACTOR', 'co_fact'),
        "AC": ('ACTIVATING_COMPOUND', 'act_comp'),
        "IN": ('INHIBITORS', 'inhib'),
        "ME": ('METALS_IONS', 'metals'),
        "MW": ('MOLECULAR_WEIGHT', 'mw'),
        "PM": ('POSTTRANSLATIONAL_MODIFICATION', 'mods'),
        "SU": ('SUBUNITS', 'sub_units'),
        "AP": ('APPLICATION', 'appl'),
        "CL": ('CLONED', 'cloned'),
        "PU": ('PURIFICATION', 'puri'),
        "OSS": ('ORGANIC_SOLVENT_STABILITY', 'solvent_stab'),
        "PHS": ('PH_STABILITY', 'ph_stab'),
        "SS": ('STORAGE_STABILITY', 'store_stab'),
        "TS": ('TEMPERATURE_STABILITY', 'T_stab'),
        "RF": ('REFERENCE', 'refs'),
        "CR": ('CRYSTALLIZATION', 'cryst'),
        "EN": ('ENGINEERING', 'eng'),
        "EXP": ('EXPRESSION', 'exp'),
        "GS": ('GENERAL_STABILITY', 'gen_stab'),
        "GI": ('GENERAL_STABILITY', 'gen_info'),
        "IC50": ('IC50_VALUE', 'ic50'),
        "KKM": ('KKM_VALUE', 'kkm'),
        "KI": ('KI_VALUE', 'ki'),
        "LO": ('LOCALIZATION', 'local'),
        "OS": ('OXIDATION_STABILITY', 'oxid_stab'),
        "PHR": ('PH_RANGE', 'ph_range'),
        "PI": ('PI_VALUE', 'pi'),
        "REN": ('RENATURED', 'renat'),
                 }

    # collect data for brenda txt file
    br_data = {}
    for i in br_symbols.keys():
        br_data[i] = []

    return br_symbols, br_data


def get_brenda_dict(br_file):
    """
    Convert brenda EC text file into dictionary.

    Arguments:
        br_file (str) - name of txt file

    Returns:
        br_symbols (dict) - dictionary of brenda symbols
        br_data (dict) - dictionary of brenda data

    """
    br_symbols, br_data = initialise_br_dict()
    # read in brenda datafile and remove all empty lines
    with open(br_file, 'r') as f:
        lines = f.readlines()
    # remove empty lines
    for l in range(len(lines)-1, -1, -1):
        if lines[l] == '\n':
            del lines[l]
    # collect a single set of lines associated with one of the known
    # 'initial' values
    ind_lines_saved = []
    for i, l in enumerate(lines):
        if i in ind_lines_saved:
            continue
        initial = l.split('\t')[0]
        if initial in br_data.keys():
            # check for wrapped lines
            aux_lines = [l.split('\t')[1]]
            ind_lines_saved.append(i)
            end_wrap = False
            # check following lines
            # 1000 is obscene - but it should stop before then.
            for ind in range(1, 1000):
                if end_wrap is True:
                    break
                if '\t' not in lines[i+ind]:
                    end_wrap = True
                    break
                else:
                    # is the initial part of this line equivalent to a known
                    # 'initial' value?
                    n_initial = lines[i+ind].split('\t')[0]
                    if n_initial in br_data.keys():
                        # yes - then end wrapping
                        end_wrap = True
                        break
                    else:
                        # no - append to previous line
                        aux_lines.append(lines[i+ind].split('\t')[1])
                        ind_lines_saved.append(i+ind)
            if end_wrap:
                final_line = aux_lines[0]
                if len(aux_lines) > 0:
                    for aux in aux_lines[1:]:
                        final_line += aux
                # check for duplicate lines in BRENDA file
                # they do exist!
                if final_line not in br_data[initial]:
                    br_data[initial].append(final_line)

    return br_symbols, br_data


def check_new_lines_and_split(string):
    """Check for new line symbols in list of strings and split string to list.

    """
    new_list = []
    for i in string:
        if '\n' in i:
            # split the item
            n_i = i.split('\n')
            for I in n_i:
                new_list.append(I)
        else:
            new_list.append(i)

    return new_list


def collect_PR_from_line(line):
    """Collect all unique protein (PR) numbers from between "#" in BRENDA.

    """
    split_l = line.split("#")
    PRs = split_l[1]
    # first_sec = line.split(" ")[0]
    # check formatting is consistent
    # if first_sec[0] != '#' or first_sec[-1] != '#':
    #    print('formatting is off -- exitting...')
    #    print('problem line:')
    #    print(line)
    #    from sys import exit
    #    exit()
    # first_sec = first_sec.replace('#', '')
    PR_codes = PRs.split(",")
    return PR_codes


def is_species_reported(spec, EC_pi_data, br_data, verbose=True):
    """Check if target species in sequences tested and BRENDA.

    """
    # is species in sequence data?
    spec_in_seq = False
    all_species = list(set(sorted(EC_pi_data['species'])))
    for i in all_species:
        if spec in i.lower():
            spec_in_seq = True
            break
    # is species in brenda data
    spec_in_br = False
    for i in sorted(br_data['PR']):
        if spec in i.lower():
            spec_in_br = True
            break
    if verbose:
        print('species in sequence data:', spec_in_seq)
        print('species in BRENDA data:', spec_in_br)
    return spec_in_seq, spec_in_br


def get_prop_PR_codes(list_of_int):
    """Get list of proteins with property of interest in BRENDA.

    """
    # get list of protein codes with property of interest
    EC_prop_PR_codes = {}
    for i, line in enumerate(list_of_int):
        PR_codes = collect_PR_from_line(line)
        if len(PR_codes) > 0:
            EC_prop_PR_codes[i] = []
        for pr in PR_codes:
            if '\n' in pr:
                # split the item
                n_i = pr.split('\n')
                for I in n_i:
                    EC_prop_PR_codes[i].append(I)
            else:
                EC_prop_PR_codes[i].append(pr)
    return EC_prop_PR_codes


def get_cmpd_information(molec):
    """Get information from SABIO Database of a compound with ID cID.

    """
    QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchCompoundDetails'

    # input: SabioCompoundID
    # valid output fields: "fields[]":["Name","ChebiID",
    #                                  "PubChemID","InChI",
    #                                  "SabioCompoundID","KeggCompoundID"]
    params = {"SabioCompoundID": molec.cID,
              "fields[]": ["Name", "ChebiID", "PubChemID", "InChI"]}
    if molec.InChi is None:
        request = requests.post(QUERY_URL, params=params)
        request.raise_for_status()
        if request.text == 'No results found for query':
            print(request.text)
            molec.mol = None
        else:
            # results
            txt = request.text.split('\n')[1].split('\t')
            _, molec.chebiID, molec.PubChemId, molec.InChi = txt

    if molec.InChi != 'null':
        molec.mol = get_rdkit_mol_from_InChi(molec.InChi)
        smiles = Chem.MolToSmiles(Chem.RemoveHs(molec.mol))
        molec.SMILES = smiles
    else:
        molec.mol = None
        molec.SMILES = None


def get_rxn_systems(EC, output_dir, clean_system=False):
    """Get reaction systems from SABIO entries in one EC and output to Pickle.

    """
    # get all SABIO entries
    entries = get_entries_per_EC(EC)
    # iterate over entries
    count = 0
    for eID in entries:
        print('DB: SABIO - EC:', EC, '-',
              'DB ID:', eID, '-', count, 'of', len(entries))
        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'SABIO', eID)
        if os.path.isfile(output_dir+rs.pkl) is True and clean_system is False:
            print('-----------------------------------')
            count += 1
            continue
        # get reaction ID
        # SABIO specific properties
        rs.organism, rs.rID, rs.UniprotID = get_rxnID_from_eID(eID)
        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.rID)
        if rs.skip_rxn is False:
            # append compound information
            for m in rs.components:
                m.get_compound()

        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        rs.save_object(output_dir+rs.pkl)
        print('-----------------------------------')
        count += 1


def get_rxn_system(rs, ID):
    """Get reaction system from SABIO reaction ID (rID).

    Uses SABIO API - online.

    Keywords:
        rs (class) - reaction system object
        ID (str) - DB reaction ID

    """
    QUERY_URL = 'http://sabiork.h-its.org/testSabio/sabioRestWebServices/searchReactionParticipants'
    # input: SabioReactionID
    # valid output fields: "fields[]":
    #    ["Name","Role","SabioCompoundID","ChebiID",
    #     "PubChemID","KeggCompoundID","InChI"]

    params = {"SabioReactionID": ID,
              "fields[]": ["Name", "Role", "SabioCompoundID"]}
    # , "ChebiID",
    # "PubChemID", "KeggCompoundID", 'UniprotID']}
    request = requests.post(QUERY_URL, params=params)
    request.raise_for_status()
    if request.text == 'No results found for query':
        print(request.text)
        rs.skip_rxn = True
        return rs
    # collate request output
    rs.components = []
    for i in request.text.split('\n')[1:]:
        if len(i) > 1:
            mol, role, cID = i.split('\t')
            new_mol = molecule.molecule(mol, role, 'SABIO', cID)
            # add new_mol to reaction system class
            rs.components.append(new_mol)
    return rs
