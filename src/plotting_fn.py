#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting functions.

Author: Andrew Tarzia

Date Created: 15 Sep 2018

"""

import matplotlib.pyplot as plt


def EC_descriptions():
    """Dictionary of EC descriptions + colours.

    """
    top_tier = {
        '-': ('unknown', '#1469b5'),
        '1': ('oxidoreductases', '#FF7900'),
        '2': ('transferases', '#00B036'),
        '3': ('hydrolases', '#EB0000'),
        '4': ('lyases', '#A440BC'),
        '5': ('isomerases', '#945348'),
        '6': ('ligases', '#FA4BBE')
    }

    return top_tier


def hit_point_plot(hit_points, ET, center, radii, rotation):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # plot points
    # atom_positions = conf.GetPositions()
    # ax.scatter(atom_positions[:, 0], atom_positions[:, 1],
    #            atom_positions[:, 2],
    #            color='k', marker='o', s=100)
    ax.scatter(
        hit_points[:, 0], hit_points[:, 1], hit_points[:, 2],
        color='g', marker='x', edgecolor=None,
        s=50, alpha=0.5
    )

    # plot ellipsoid
    ET.plotEllipsoid(
        center, radii, rotation, ax=ax, plotAxes=False
    )
    ax.set_xlabel(r"$x$ [$\mathrm{\AA}$]", fontsize=16)
    ax.set_ylabel(r"$y$ [$\mathrm{\AA}$]", fontsize=16)
    ax.set_zlabel(r"$z$ [$\mathrm{\AA}$]", fontsize=16)
    # ax.set_xlim(-max(radii*2), max(radii*2))
    # ax.set_ylim(-max(radii*2), max(radii*2))
    # ax.set_zlim(-max(radii*2), max(radii*2))
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-10, 10)
    ax.set_aspect('equal', 'box')
    plt.axis('off')
    # ax.grid(False)
    # # ax.xaxis.pane.set_edgecolor('black')
    # # ax.yaxis.pane.set_edgecolor('black')
    # # ax.zaxis.pane.set_edgecolor('black')
    # ax.xaxis.set_major_locator(MultipleLocator(2))
    # ax.yaxis.set_major_locator(MultipleLocator(2))
    # ax.zaxis.set_major_locator(MultipleLocator(2))
    # ax.xaxis.pane.fill = False
    # ax.yaxis.pane.fill = False
    # ax.zaxis.pane.fill = False
    dist = [30, 30]
    angles = [-90, -180]
    for i, j in zip(dist, angles):
        ax.view_init(i, j)
        fig.tight_layout()
        fig.savefig(
            'temporary_'+str(i)+'_'+str(j)+'.pdf', dpi=720,
            bbox_inches='tight'
        )


def define_diff_categ_plot(ax, title, ytitle, xtitle, xlim, ylim):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_ylabel(ytitle, fontsize=16)
    # ax.legend(
    #     [y, n], ['aligned', 'not aligned'], loc=4, fancybox=True
    # )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticklabels(['diffuses', 'does not diffuse'])
    ax.set_xticks([0.25, 0.75])


def define_standard_plot(
    ax,
    ytitle,
    xtitle,
    title=None,
    xlim=None,
    ylim=None
):
    """
    Series of matplotlib pyplot settings to make all plots unitform.

    """

    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if title is not None:
        ax.set_title(title, fontsize=16)


def define_3d_plot(
    ax, title, xtitle, ytitle, ztitle, xlim, ylim, zlim
):
    """
    Series of matplotlib pyplot settings to make all plots unitform.
    """
    # Set number of ticks for x-axis
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.set_xlabel(xtitle, fontsize=16)
    ax.set_ylabel(ytitle, fontsize=16)
    ax.set_ylabel(ztitle, fontsize=16)
    # ax.legend(
    #     [y, n], ['aligned', 'not aligned'], loc=4, fancybox=True
    # )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_ylim(zlim)
