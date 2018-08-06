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
        br_lines = f.readlines()
    
    # remove empty lines and use '___' as delimiter.
    br_lines = [i.rstrip() for i in br_lines]
    br_lines_new = []
    for i in br_lines:
        if i != '':
            br_lines_new.append(i.replace('\t', '___'))
    
    # append lines to appropriate dictionary lists
    for i, line in enumerate(br_lines_new):
        if "___" in line:
            l_split = line.split("___")
            initial = l_split[0]
            if initial in br_data.keys():
                if l_split[1] not in br_data[initial]:
                    br_data[initial].append(l_split[1])
            else:
                # implies this line is spill over from the last line
                # so we add to that data
                prev_initial = br_lines_new[i-1].split("___")[0]
                if prev_initial in br_data.keys():
                    prev_data = br_data[prev_initial][-1]
                    prev_data += l_split[1]
                    br_data[prev_initial][-1] = prev_data
    
    
    return br_symbols, br_data