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
import reaction
import utilities


def main():
    if (not len(sys.argv) == 3):
        print('Usage: screening.py plot_suffix target_EC\n')
        print("""
        plot_suffix: string to put at the end of plot file names.
        target_EC: file containing a target EC number (f if all EC
            should be used)
        param_file:

        """)
        sys.exit()
    else:
        plot_suffix = sys.argv[1]
        target_EC_file = sys.argv[2] if sys.argv[2] != 'f' else None
        pars = utilities.read_params(sys.argv[3])

    if target_EC_file is None:
        # Get all EC numbers.
        search_EC_file = pars['EC_file']
        search_ECs = utilities.get_ECs_from_file(
            EC_file=search_EC_file
        )
    else:
        # Read ECs from file.
        search_ECs = utilities.get_ECs_from_file(
            EC_file=target_EC_file
        )

    print(search_ECs)

    # Iterate through all reactions in directory.
    search_output_dir = os.get_cwd()
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

    print(output_data)
    target_data = output_data[output_data['ec'] in search_ECs]
    print(target_data)
    print(
        len(target_data),
        len(output_data),
        len(target_data)/len(output_data)
    )

    sys.exit()

    pr.no_rxns_vs_size(
        data=target_data,
        params=pars,
        ECs=search_ECs,
        plot_suffix=plot_suffix
    )
    pr.save_candidates(
        data=target_data,
        params=pars,
        filename=f'candidates_{plot_suffix}.csv'
    )

    sys.exit()

    plots_to_do = [
        # Column, stacked, xtitle, xlim, width
        ('minlogs', False, 'min. logS', None, 0.5),
        ('minlogs', True, 'min. logS', None, 0.5),
        ('maxlogp', False, 'max. logP', None, 0.5),
        ('maxlogp', True, 'max. logP', None, 0.5),
        ('nr', False, 'no reactants', None, 0.5),
        ('np', False, 'no products', None, 0.5),
        (
            'max_mid_diam',
            False,
            r'$d$ of largest component [$\mathrm{\AA}$]',
            None,
            0.5
        ),
        (
            'max_mid_diam',
            True,
            r'$d$ of largest component [$\mathrm{\AA}$]',
            None,
            0.5
        ),
        ('deltasa', False, r'$\Delta$ SAscore', (-10, 10), 0.5),
        ('deltasa', True, r'$\Delta$ SAscore', (-10, 10), 0.5),
    ]

    for pl in plots_to_do:
        col, stacked, xtitle, xlim, width = pl
        if stacked:
            fig, ax = pr.stacked_dist(
                data=target_data,
                col=col,
                xtitle=xtitle,
                xlim=xlim,
                width=width
            )

            fig.tight_layout()
            fig.savefig(
                f"stacked_{col}_{plot_suffix}.pdf",
                dpi=720,
                bbox_inches='tight'
            )
            fig, ax = pr.violinplot(
                data=target_data,
                col=col,
                xtitle=xtitle,
            )
            fig.tight_layout()
            fig.savefig(
                f"violin_{col}_{plot_suffix}.pdf",
                dpi=720,
                bbox_inches='tight'
            )
        else:
            fig, ax = pr.dist(
                X=target_data[col],
                xtitle=xtitle,
                xlim=xlim,
                width=width
            )
            fig.tight_layout()
            fig.savefig(
                f"dist_{col}_{plot_suffix}.pdf",
                dpi=720,
                bbox_inches='tight'
            )


if __name__ == "__main__":
    main()
