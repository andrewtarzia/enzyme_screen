#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for I/O of BRENDA DB.

Author: Andrew Tarzia

Date Created: 24 Apr 2018

License:


"""


def initialise_br_dict():
    """
    Initialise the brenda data dictionary and the associated symbol definitions.

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
            EC_prop_PR_codes[i].append(pr)
    return EC_prop_PR_codes
