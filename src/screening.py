#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script for plotting data that is used for screening.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import os
import sys
import pandas as pd

import plots_reaction as pr
import utilities


def case_studies(string, pars):

    CS = {
        'production': {
            'file_suffix': 'production',
        },
        'biomin': {
            'file_suffix': 'biomin',
            'EC_file': os.path.join(
                os.path.dirname(__file__),
                '../data/desired_EC_biomin.txt'
            ),
        },
        'AO': {
            'file_suffix': 'AO',
            'EC_file': os.path.join(
                os.path.dirname(__file__),
                '../data/desired_EC_AO.txt'
            ),
        },
        'linb': {
            'file_suffix': 'linb',
            'EC_file': os.path.join(
                os.path.dirname(__file__),
                '../data/desired_EC_linb.txt'
            ),
        },
        'esterases': {
            'file_suffix': 'esterases',
            'EC_file': os.path.join(
                os.path.dirname(__file__),
                '../data/desired_EC_esterases.txt'
            ),
        },
        'lipases': {
            'file_suffix': 'lipases',
            'EC_file': os.path.join(
                os.path.dirname(__file__),
                '../data/desired_EC_lipases.txt'
            ),
        },
    }

    if string is None:
        print(f'{CS.keys()}')
        return
    elif string not in CS:
        raise KeyError(
            f'{string} not a defined case study. Options: {CS.keys()}'
        )

    for i in CS[string]:
        pars[i] = CS[string][i]

    return pars


def main():
    if (not len(sys.argv) == 3):
        print('Usage: screening.py param_file case\n')
        print("""
        param_file (str) :
        case (str) :
            define the case study to use.
        """)
        case_studies(None, None)
        sys.exit()
    else:
        pars = utilities.read_params(sys.argv[1])
        case = sys.argv[2]

    pars = case_studies(string=case, pars=pars)

    if not os.path.exists(pars['file_suffix']):
        os.mkdir(pars['file_suffix'])

    # Get all EC numbers.
    search_EC_file = pars['EC_file']
    search_ECs = utilities.get_ECs_from_file(
        EC_file=search_EC_file
    )

    # Iterate through all reactions in directory.
    search_output_dir = os.getcwd()
    prop_output_file = os.path.join(
        search_output_dir,
        'rs_properties.csv'
    )

    if os.path.exists(prop_output_file):
        output_data = pd.read_csv(prop_output_file)
    else:
        raise FileNotFoundError(
            f'{prop_output_file} with all data is missing'
        )

    target_data = pd.DataFrame(columns=output_data.columns)
    for i, row in output_data.iterrows():
        if row['ec'] in search_ECs:
            target_data = target_data.append(row)
    pr.no_rxns_vs_size(
        data=target_data,
        params=pars,
        plot_suffix=pars['file_suffix']
    )
    pr.save_candidates(
        data=target_data,
        params=pars,
        filename=(
            f"{pars['file_suffix']}/"
            f"candidates_{pars['file_suffix']}.csv"
        )
    )

    plots_to_do = [
        # Column, type (dist, stacked, pie), xtitle, xlim, width
        ('minlogs', 'dist', 'min. logS', None, 0.5),
        ('minlogs', 'stacked', 'min. logS', None, 0.5),
        ('maxlogp', 'dist', 'max. logP', None, 0.5),
        ('maxlogp', 'stacked', 'max. logP', None, 0.5),
        ('nr', 'dist', 'no reactants', None, 0.5),
        ('np', 'dist', 'no products', None, 0.5),
        ('PC_class', 'pie', 'purchasability class', None, 1),
        (
            'max_mid_diam',
            'dist',
            r'$d$ of largest component [$\mathrm{\AA}$]',
            None,
            0.5
        ),
        (
            'max_mid_diam',
            'stacked',
            r'$d$ of largest component [$\mathrm{\AA}$]',
            None,
            0.5
        ),
        ('deltasa', 'dist', r'$\Delta$ SAscore', (-10, 10), 0.5),
        ('deltasa', 'stacked', r'$\Delta$ SAscore', (-10, 10), 0.5),
    ]

    for pl in plots_to_do:
        col, type, xtitle, xlim, width = pl
        if type == 'stacked':
            fig, ax = pr.stacked_dist(
                data=target_data,
                col=col,
                xtitle=xtitle,
                xlim=xlim,
                width=width
            )

            fig.tight_layout()
            fig.savefig(
                fname=(
                    f"{pars['file_suffix']}/"
                    f"stacked_{col}_{pars['file_suffix']}.pdf"
                ),
                dpi=720,
                bbox_inches='tight'
            )
            fig, ax = pr.violinplot(
                data=target_data,
                col=col,
                ytitle=xtitle,
                ylim=xlim
            )
            fig.tight_layout()
            fig.savefig(
                fname=(
                    f"{pars['file_suffix']}/"
                    f"violin_{col}_{pars['file_suffix']}.pdf"
                ),
                dpi=720,
                bbox_inches='tight'
            )
        elif type == 'dist':
            fig, ax = pr.dist(
                X=target_data[col],
                xtitle=xtitle,
                xlim=xlim,
                width=width
            )
            fig.tight_layout()
            fig.savefig(
                fname=(
                    f"{pars['file_suffix']}/"
                    f"dist_{col}_{pars['file_suffix']}.pdf"
                ),
                dpi=720,
                bbox_inches='tight'
            )
        elif type == 'pie':
            fig, ax = pr.pie(
                X=target_data[col],
                xtitle=xtitle,
                xlim=xlim,
                width=width
            )
            fig.tight_layout()
            fig.savefig(
                fname=(
                    f"{pars['file_suffix']}/"
                    f"pie_{col}_{pars['file_suffix']}.pdf"
                ),
                dpi=720,
                bbox_inches='tight'
            )


if __name__ == "__main__":
    main()
