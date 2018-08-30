#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions to go from RDKIT molecule to an ellipsoid that encompasses it's
VDW cloud.

Author: Andrew Tarzia

Date Created: 27 Aug 2018

"""

import ellipsoid
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Geometry import rdGeometry


def get_molec_shape(mol, conf, confId, vdwScale=1.0,
                    boxMargin=2.0, spacing=0.2):
    """Get the shape of a conformer of a molecule as a grid representation.

    """
    box = Chem.ComputeConfBox(conf)
    sideLen = (box[1].x-box[0].x + 2*boxMargin,
               box[1].y-box[0].y + 2*boxMargin,
               box[1].z-box[0].z + 2*boxMargin)
    shape = rdGeometry.UniformGrid3D(2*sideLen[0],
                                     2*sideLen[1],
                                     2*sideLen[2],
                                     spacing=spacing)
    Chem.EncodeShape(mol, shape, confId=confId, ignoreHs=False,
                     vdwScale=vdwScale)
    return box, sideLen, shape


def coord_to_RDKITMOL(coords):
    """Function that converts 3D coordinates into an RDKIT Molecule

    """
    print("I have never done this - but my suggestion is that RDKIT should")
    print("have some function that does exactly this process for you.")
    print("ASE might do the job.")
    print("NOTE: I think you need to assign your coordinates to a CONF with ")
    print("CONFID = 0 within the MOLECULE object")
    return False


def calc_molecule_diameters(RDKIT_MOL, vdwScale=1.0,
                            boxMargin=4.0, spacing=1.0):
    """Calculate the diameters of an RDKIT Molecule.

    Keywords:
        RDKIT_MOL () - RDKIT Molecule you want to get the ellipsoid of
        vdwScale (float) - set the scale of the VDW radii (RDKIT VDW radii are
            over esimates - one search in their GIT repo issues section will
            show a thread about this) - defaults to 1.0
        boxMargin (float) - sets extent of the grid (check results for
            sensitivity to this)
        spacing (float) - sets the grid spacing for VDW cloud. (check results
            for sensitivity to this)
            Speed scales with this, memory should be ok.

    Returns:
        diameters (list) - three principal diameters of ellipsoid
            (min, mid, max)


    """
    box, sideLen, shape = get_molec_shape(RDKIT_MOL,
                                          vdwScale=vdwScale,
                                          boxMargin=boxMargin,
                                          spacing=spacing)

    # get ellipsoid fitting all points with value > 2
    # - i.e. within vdw shape
    hit_points = []
    for idx in range(shape.GetSize()):
        pt = shape.GetGridPointLoc(idx)
        value = shape.GetVal(idx)
        if value > 2:
            point = np.array([pt.x, pt.y, pt.z])
            hit_points.append(point)
    hit_points = np.asarray(hit_points)

    # find the ellipsoid that envelopes all hit points
    ET = ellipsoid.EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(hit_points, .01)

    diameters = sorted(np.asarray(radii)*2)

    # this can be modified to plot the ellipsoid and all hit points
    if False:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # plot points
        # atom_positions = conf.GetPositions()
        # ax.scatter(atom_positions[:,0], atom_positions[:,1],
        #            atom_positions[:,2],
        #            color='k', marker='o', s=100)
        ax.scatter(hit_points[:, 0], hit_points[:, 1], hit_points[:, 2],
                   color='g', marker='*', s=100)

        # plot ellipsoid
        ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False)

        ax.set_xlabel("$X$")
        ax.set_ylabel("$Y$")
        ax.set_zlabel("$Z$")
        ax.set_xlim(-max(radii*2), max(radii*2))
        ax.set_ylim(-max(radii*2), max(radii*2))
        ax.set_zlim(-max(radii*2), max(radii*2))
        # ax.set_aspect('equal','box')
        if input('save fig?') is 't':
            fig.tight_layout()
            fig.savefig(input("file name?"), dpi=720,
                        bbox_inches='tight')
        else:
            plt.show()

    return diameters


if __name__ == "__main__":
    # read in coordinates somehow
    coords = 0
    RDKIT_MOL = coord_to_RDKITMOL(coords)
    # get ellipsoid diameters
    min_diam, mid_diam, max_diam = calc_molecule_diameters(RDKIT_MOL,
                                                           vdwScale=1.0,
                                                           boxMargin=4.0,
                                                           spacing=1.0)
    print('ellipsoid diameters:')
    print('min:', min_diam, 'angstrom')
    print('mid:', mid_diam, 'angstrom')
    print('max:', max_diam, 'angstrom')
