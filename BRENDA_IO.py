#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for I/O of BRENDA DB.

Author: Andrew Tarzia

Date Created: 24 Apr 2018

"""

import molecule
import rxn_syst
import os
import DB_functions
import CHEBI_IO


def extract_subunit_info(br_data, PR):
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

    no_sub_units = output_nprop

    return no_sub_units


def extract_PTM(br_data, PR):
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

    PTM = output_nprop

    return PTM


def extract_cofactor_info(br_data, PR):
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

    cofactor_mol = output_nprop
    cofactors = True

    return cofactors, cofactor_mol


def extract_activating_mol(br_data, PR):
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

    activating_mol = output_nprop

    return activating_mol


def extract_general_rxn_info(br_data):
    """Extract all activating molecules for the PR of this reaction system.
        It is possible that there is no information, or multiple entries.

    """
    nprop = 'RT'
    output_nprop = []
    list_of_nprop = br_data[nprop]
    for i, entry in enumerate(list_of_nprop):
        value = entry.rstrip()
        output_nprop.append(value)

    reaction_type = output_nprop
    nprop = 'RE'
    output_nprop = []
    list_of_nprop = br_data[nprop]
    for i, entry in enumerate(list_of_nprop):
        value = entry.rstrip()
        output_nprop.append(value)

    reaction_cat = output_nprop

    return reaction_type, reaction_cat


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


def get_chebiID_for_BRENDA(mol_name, ont):
    """Convert molecule name to chebiID using CHEBI Ontology files.

    Offline.

    Keywords:
        mol_name (str) - molecule name
        ont (pront ontology object) - CHEBI ontology

    Returns:
        ID (str) - chebiID
    """

    DB_prop = DB_functions.get_DB_prop('CHEBI')
    compounds_file = DB_prop[0]+DB_prop[1]['cmpds_file']
    names_file = DB_prop[0]+DB_prop[1]['names_file']
    structures_file = DB_prop[0]+DB_prop[1]['strct_file']

    # search for name in compound file
    res = CHEBI_IO.search_for_compound_by_name(compounds_file, mol_name)
    if res is None:
        # search for formula in names file
        res = CHEBI_IO.search_for_name_by_name(names_file, mol_name)
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
        res = CHEBI_IO.convert_nameID_to_parent(compounds_file, nameID=ID)
        if res is None:
            print("this should not happen - error with cross reference")
            print('check this!')
            import sys
            sys.exit()
        ID, parent_id, name, star = res

    return ID


def define_BRENDA_file(EC):
    """Define the location of BRENDA data files.

    """
    DB_prop = DB_functions.get_DB_prop('BRENDA')
    file = DB_prop[0]+'brenda_download_'+EC.replace('.', '_')+'.txt'
    return file


def get_rxn_systems(EC, output_dir, clean_system=False):
    """Get reaction systems from BRENDA file of one EC and output to Pickle.

    """
    # read in CHEBI ontology file for usage
    ont = 1  # CHEBI_IO.read_ontology()
    # get EC BRENDA data
    br_datafile = define_BRENDA_file(EC)
    br_sym, br_data = get_brenda_dict(br_datafile)

    entries = br_data['SP']

    # iterate over entries
    count = 0
    eID_c = 1
    for e in entries:
        eID = 'BR'+str(eID_c)
        eID_c += 1
        print('DB: BRENDA - EC:', EC, '-',
              'DB ID:', eID, '-', count, 'of', len(entries))
        # initialise reaction system object
        rs = rxn_syst.reaction(EC, 'BRENDA', eID)
        if os.path.isfile(output_dir+rs.pkl) is True and clean_system is False:
            print('-----------------------------------')
            count += 1
            continue
        # BRENDA specific properties
        # items collected in get_rxn_system:
            # associated PR codes -- rs.assoc_PR
            #   (numbers for linking within BRENDAfiles )
            # associated references -- rs.assoc_refs
            #   (numbers for linking within BRENDAfiles )
            # meta -- rs.meta
            # reversibility -- rs.reversible (bool)
        # removed currently because it unnecessarily complicates the code.
        # rs.activating_compounds = 1  # defines a list of activating compounds
        # rs.cofactors = 1  # defines a list of cofactors

        # get reaction system using DB specific function
        rs = get_rxn_system(rs, rs.DB_ID, e, ont)
        if rs.skip_rxn is False:
            # append compound information
            for m in rs.components:
                m.get_compound()

        # pickle reaction system object to file
        # prefix (sRS for SABIO) + EC + EntryID .pkl
        print('skip?', rs.skip_rxn)
        rs.save_object(output_dir+rs.pkl)
        print('-----------------------------------')
        count += 1

    # remove ontology data from memory
    del ont


def get_rxn_system(rs, ID, entry, ont):
    """Get reaction system from BRENDA data.

    Offline.

    Keywords:
        rs (class) - reaction system object
        ID (str) - DB reaction ID
        entry (str) - entry in BRENDA data file
        ont (pront ontology object) - CHEBI ontology

    """
    # string manipulations
    entry = entry.replace('\n', ' ')
    PR_sect = entry.split("# ")[0]+"#"
    entry_2 = entry.replace(PR_sect, '')
    rxn_sect = entry_2.split(" ( ")[0].split(" (#")[0].split(" <")[0]
    entry_3 = entry_2.replace(rxn_sect, '')
    meta_sect = entry_3

    # collect some BRENADA specific information
    t_assoc_PR = PR_sect.split("#")[1].split(',')
    # check for new lines and split string into list
    rs.assoc_PR = check_new_lines_and_split(t_assoc_PR)
    rs.meta = meta_sect
    # references for an entry in BRENDA are within "<" and ">"
    t_assoc_refs = meta_sect.split(" <")[-1]
    t_assoc_refs = t_assoc_refs.split(">")[0].split(',')
    # check for new lines and split string into list
    rs.assoc_refs = check_new_lines_and_split(t_assoc_refs)
    # reversible?
    # for the "SP" entries - a 'r' enclosed in "{" and "}"
    # implies reversible
    if '{r}' in rs.meta:
        rs.reversible = True

    # skip reactions with unknown components
    if '?' in rxn_sect:
        rs.skip_rxn = True
        return rs

    # get reactants and products
    react, prod = rxn_sect.replace("{", "").replace("}", "").split(" = ")
    print(react, prod)
    # separate react and prod into molecules by "+"
    r_mol = react.split(" + ")
    p_mol = prod.split(" + ")
    # remove preceding and succeeding white space from all molecule names
    r_mol = [i.lstrip().rstrip() for i in r_mol]
    p_mol = [i.lstrip().rstrip() for i in p_mol]

    # remove stoichiometry
    new_r = []
    for r in r_mol:
        # have stoich?
        if r.split(' ')[0].isnumeric() is True:
            # yes
            new_r.append(' '.join(r.split(' ')[1:]))
        else:
            # no
            new_r.append(r)
    new_p = []
    for r in p_mol:
        # have stoich?
        if r.split(' ')[0].isnumeric() is True:
            # yes
            new_p.append(' '.join(r.split(' ')[1:]))
        else:
            # no
            new_p.append(r)

    print(new_r, new_p)
    # define component list
    comp_list = []
    for i in new_r:
        comp_list.append((i, 'reactant'))
    for i in new_p:
        comp_list.append((i, 'product'))

    rs.components = []
    for comp in comp_list:
        chebiID = get_chebiID_for_BRENDA(comp[0], ont)
        print(comp[0], chebiID)
        if chebiID is None:
            rs.skip_rxn = True
            continue
        new_mol = molecule.molecule(comp[0], comp[1], 'BRENDA', chebiID)
        # add new_mol to reaction system class
        rs.components.append(new_mol)

    return rs
