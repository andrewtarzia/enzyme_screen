#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to plot pI of all sequences in BRENDA flat files.

For my thesis.

Author: Andrew Tarzia

Date Created: 21 Dec 2018

"""
import numpy as np
import matplotlib.pyplot as plt
from ercollect.bulk_protein_analysis import read_seq_output


if __name__ == "__main__":
    import sys
    files = ['1__BRENDA_sequences_output.csv',
             '2__BRENDA_sequences_output.csv',
             '3__BRENDA_sequences_output.csv',
             '4__BRENDA_sequences_output.csv',
             '5__BRENDA_sequences_output.csv',
             '6__BRENDA_sequences_output.csv',
             ]
    pI_list = []
    for file in files:
        print(file)
        output = read_seq_output(file)
        for i, j in zip(list(output.pI), list(output.note)):
            if j == 'Swiss-Prot':
                pI_list.append(i)
        del output  # save some memory
        print('done - ', len(pI_list))
    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    # ax.hist(pI_list, facecolor='purple', alpha=0.6,
    #         histtype='stepfilled', density=True, edgecolor='k',
    #         bins=np.arange(0, 14 + 0.2, 0.5))
    width = 0.5
    X_bins = np.arange(0, 14 + 0.2, width)
    hist, bin_edges = np.histogram(a=pI_list,
                                   bins=X_bins,
                                   density=True)
    ax.bar(bin_edges[:-1],
           hist,
           align='edge',
           alpha=0.6, width=width,
           color='purple',
           edgecolor='k')
    ax.axvline(x=7, c='k', lw=2)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('pI', fontsize=16)
    ax.set_ylabel('proportion', fontsize=16)
    ax.set_xlim(0, 14)
    # legend
    # ax.legend(fontsize=16)
    fig.tight_layout()
    fig.savefig("dist_pI_full.pdf",
                dpi=720, bbox_inches='tight')
    sys.exit()
