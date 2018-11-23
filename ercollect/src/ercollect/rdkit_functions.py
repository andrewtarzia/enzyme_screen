#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Functions for RDKIT usage.

Author: Andrew Tarzia

Date Created: 16 Jul 2018

"""

import numpy as np
import os
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, Draw, PyMol
from rdkit.Chem.Descriptors3D import NPR1, NPR2, PMI1, PMI2, PMI3
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Geometry import rdGeometry
from rdkit import Geometry
import tempfile
from ercollect import ellipsoid
from matplotlib.ticker import MultipleLocator


def draw_svg_for_all_molecules(molecules, output_dir):
    """Draw an SVG image for each molecule in a dictionary of molecules.

    {name: SMILES}

    """
    for key, val in molecules.items():
        try:
            draw_smiles_to_svg(
                val, output_dir+key.replace(' ', '_')+'_2d.svg')
        except UnicodeEncodeError:
            pass


def calculate_all_MW(molecules):
    """Calculate the molecular weight of all molecules in DB dictionary.

    {name: SMILES}

    """
    for m, smile in molecules.items():
        # Read SMILES and add Hs
        mol = Chem.AddHs(Chem.MolFromSmiles(smile))
        MW = Descriptors.MolWt(mol)
        print(m, '---', smile, '---', 'MW =', MW, 'g/mol')


def draw_smiles_to_svg(smiles, filename):
    """Draw a single molecule to an SVG file with transparent BG.

    """
    mol = Chem.MolFromSmiles(smiles)
    # change BG to transperent
    # (https://sourceforge.net/p/rdkit/mailman/message/31637105/)
    o = DrawingOptions()
    o.bgColor = None
    Chem.Compute2DCoords(mol)
    Draw.MolToFile(mol, filename, fitImage=True, imageType='svg',
                   options=o)


def read_mol_txt_file(filename):
    """Function to read molecule SMILES and information from txt file.

    """
    data = pd.read_table(filename, delimiter=':')
    molecules = {}
    diameters = {}
    for i, row in data.iterrows():
        # try:
        #     name, smile, radius = line.rstrip().split(':')
        # except ValueError:
        #     print(line, 'had : in there twice, fix this naming or SMILE')
        #     print('skipped')
        name = row['molecule']
        smile = row['smile']
        diameter = row['diameter']
        molecules[name] = smile
        diameters[name] = diameter
    return data, molecules, diameters


def produce_quick_fig_mol(molecules, filename, labels=True, mpr=5, ims=200):
    """Produce a quick/dirty figure showing all the 2D coordinates of molecules
        in the data set.

    """
    DrawingOptions.bondLineWidth = 1.8
    DrawingOptions.atomLabelFontSize = 16
    mols = [Chem.MolFromSmiles(x) for x in molecules.values()]
    for m in mols: tmp = Chem.Compute2DCoords(m)
    # Draw.MolToFile(mols[0], output_dir+'mol1.png')
    if len(mols) > 20:
        im = 1
        M = []
        for i in mols:
            M.append(i)
            if len(M) == 20:
                if labels:
                    img = Draw.MolsToGridImage(M, molsPerRow=mpr,
                                               subImgSize=(ims, ims),
                                               legends=[x for x in molecules.keys()])
                else:
                    img = Draw.MolsToGridImage(M, molsPerRow=mpr,
                                               subImgSize=(ims, ims))
                out_name = filename.replace('.pdf', '_'+str(im)+'.pdf')
                print(out_name)
                img.save(out_name)
                im += 1
                M = []
        # final figure with remaining
        if len(M) > 0:
            if labels:
                img = Draw.MolsToGridImage(M, molsPerRow=mpr,
                                           subImgSize=(ims, ims),
                                           legends=[x for x in molecules.keys()])
            else:
                img = Draw.MolsToGridImage(M, molsPerRow=mpr,
                                           subImgSize=(ims, ims))
            out_name = filename.replace('.pdf', '_'+str(im)+'.pdf')
            print(out_name)
            img.save(out_name)

    else:
        out_name = filename
        M = mols
        if labels:
            img = Draw.MolsToGridImage(M, molsPerRow=mpr,
                                       subImgSize=(ims, ims),
                                       legends=[x for x in molecules.keys()])
        else:
            img = Draw.MolsToGridImage(M, molsPerRow=mpr,
                                       subImgSize=(ims, ims))
        img.save(out_name)


def get_inertial_prop(mol, cids):
    """Get inertial 3D descriptors for all conformers in mol.

    """
    # ratio 1 is I1/I3
    # ratio 2 is I2/I3
    sml_PMI, mid_PMI, lge_PMI = [], [], []
    ratio_1_, ratio_2_ = [], []
    for cid in cids:
        sml_PMI.append(PMI1(mol, confId=cid))
        mid_PMI.append(PMI2(mol, confId=cid))
        lge_PMI.append(PMI3(mol, confId=cid))
        ratio_1_.append(NPR1(mol, confId=cid))
        ratio_2_.append(NPR2(mol, confId=cid))

    return sml_PMI, mid_PMI, lge_PMI, ratio_1_, ratio_2_


def get_COMs(mol, cids):
    """Get COM of all conformers of mol.

    Code from:
    https://iwatobipen.wordpress.com/2016/08/16/scoring-3d-diversity-using-rdkit-rdkit/

    """
    coms = []
    numatoms = mol.GetNumAtoms()
    for confId in range(len(cids)):
        # print('conf:', confId)
        # print('number of atoms:', numatoms)
        conf = mol.GetConformer(confId)
        coords = np.array([list(conf.GetAtomPosition(atmidx)) for atmidx in range(numatoms)])
        # print('coords:')
        # print(coords)
        atoms = [atom for atom in mol.GetAtoms()]
        mass = Descriptors.MolWt(mol)
        # print('mass:', mass)
        centre_of_mass = np.array(np.sum(atoms[i].GetMass() * coords[i] for i in range(numatoms))) / mass
        # print(centre_of_mass)
        coms.append(centre_of_mass)

    return coms


def show_all_conformers(viewer, mol, cids):
    """Show all conformers in a pymol viewer.

    Code from: http://nbviewer.jupyter.org/gist/greglandrum/4316435

    """
    viewer.DeleteAll()
    for cid in cids:
        viewer.ShowMol(mol, confId=cid, name='Conf-%d' % cid, showOnly=False)


def show_axes(mol, confId, mol_coms, conf_axes):
    """Show the principal moments of inertia on a single conformation.

    """
    # show the principle axes on one conformer
    try:
        v = PyMol.MolViewer()
    except ConnectionRefusedError:
        print("run 'pymol -R' to visualise structures")
        print('visualisation will be skipped')
    v.server.do('reinitialize')
    v.server.do('delete all')

    # v.DeleteAll()
    v.ShowMol(mol, confId=confId, name='Conf-0', showOnly=False)
    # com
    v.server.sphere(list([float(i) for i in mol_coms[confId]]), 0.2, (2, 1, 0), 'COM')
    # principal axes
    v.server.sphere(list([float(i) for i in conf_axes[confId][0, :]]), 0.2, (2, 0, 0), 'AX1')
    v.server.sphere(list([float(i) for i in conf_axes[confId][1, :]]), 0.2, (1, 0, 1), 'AX2')
    v.server.sphere(list([float(i) for i in conf_axes[confId][2, :]]), 0.2, (0, 1, 0), 'AX3')
    v.GetPNG()


def show_shape(viewer, mol, cid, shape):
    """Show the encoded shape for a conformer (cid) of a molecule in viewer.

    """
    viewer.server.deleteAll()
    # write shape to file
    tmpFile = tempfile.mktemp('.grd')
    Geometry.WriteGridToFile(shape, tmpFile)
    viewer.ShowMol(mol, name='testMol', showOnly=True, confId=cid)
    viewer.server.loadSurface(tmpFile, 'testGrid', '', 2.5)


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


def def_point(x, y, z):
    """Define a 3D point in RDKIT

    """
    point = rdGeometry.Point3D()
    point.x = x
    point.y = y
    point.z = z

    return point


def define_vector(axis, max_side_len, vec_spacing=0.05):
    """Define vector to test shape values for.

    """
    vector_mag = np.arange(-max_side_len/2, max_side_len/2 + 0.1, vec_spacing)
    vectors = [i * axis for i in vector_mag]
    return vectors


def define_plane(AX1, AX2, pt, limit, vec_spacing):
    """Define a plane of interest based on two axes and a point in space.

    """
    v1 = np.asarray(define_vector(AX1, limit, vec_spacing=vec_spacing)) + np.asarray(pt)
    v2 = np.asarray(define_vector(AX2, limit, vec_spacing=vec_spacing)) + np.asarray(pt)
    return v1, v2


def get_boundary_idx(values):
    """Get index of shape grid boundaries along axis.

    """
    ind_1, ind_2 = None, None
    for i in np.arange(len(values)):
        # get the index of the first non zero value
        if values[i] > 2:
            ind_1 = i
            break
    for i in np.arange(len(values)-1, -1, -1):
        # get the index of the last non zero value
        if values[i] > 2:
            ind_2 = i
            break
    return ind_1, ind_2


def get_dist_and_values(vectors, com_pt, shape):
    """

    """
    values = []
    distances = []
    for p in vectors:
        pt = def_point(*p)
        distances.append(pt.Distance(com_pt))
        values.append(shape.GetValPoint(pt))

    return values, distances


def get_vdw_diameters(mol, cids, mol_coms, vdwScale=1.0, boxMargin=2.0,
                      spacing=0.2, vec_spacing=0.05, show=False, plot=False):
    """Get the extent of the VDW size of each conformer along its principle axes.


    """
    # over conformers
    conf_diameters = []
    conf_axes = []
    conf_moments = []
    for confId in cids:
        conf = mol.GetConformer(confId)
        print(confId)
        # try:
        #     Chem.CanonicalizeConformer(conf)
        # except RuntimeError:
        #     pass
        box, sideLen, shape = get_molec_shape(mol, conf, confId,
                                              vdwScale=vdwScale,
                                              boxMargin=boxMargin,
                                              spacing=spacing)
        if show is True:
            try:
                v = PyMol.MolViewer()
            except ConnectionRefusedError:
                print("run 'pymol -R' to visualise structures")
                print('visualisation will be skipped')
            show_shape(v, mol, confId, shape)
            v.server.do('set transparency=0.5')

        # get extent of shape along principle axes
        axes, moments = Chem.ComputePrincipalAxesAndMoments(conf,
                                                            ignoreHs=False)
        conf_axes.append(axes)
        conf_moments.append(moments)
        sml_PMI, mid_PMI, lge_PMI = moments
        COM = mol_coms[confId]
        com_pt = def_point(*COM)
        # define vector from COM along each principal axis
        max_side_len = max([i*2 for i in sideLen])
        diameters = []
        vectors = []
        vals = []
        for AX in [0, 1, 2]:
            axis = axes[AX, :]
            vector = define_vector(axis, max_side_len, vec_spacing)
            vectors.append(vector)
            values, distances = get_dist_and_values(vector, com_pt, shape)
            vals.append(values)
            ind_1, ind_2 = get_boundary_idx(values)
            print(ind_1, ind_2)
            if ind_1 is not None or ind_2 is not None:
                diameter = distances[ind_1] + distances[ind_2]
                diameters.append(diameter)
            else:
                try:
                    v = PyMol.MolViewer()
                except ConnectionRefusedError:
                    print("run 'pymol -R' to visualise structures")
                    print('visualisation will be skipped')
                show_shape(v, mol, confId, shape)
                print(max_side_len)
                print(com_pt)
                print(axis)
                v.server.do('set transparency=0.5')
                v.server.sphere(list([float(i) for i in com_pt]), 0.2, (2, 1, 0), 'COM')
                # principal axes
                v.server.sphere(list([float(i) for i in axis]), 0.2, (2, 0, 0), 'AX')
                import sys
                sys.exit()
        if len(diameters) == 3:
            conf_diameters.append(diameters)

        # only do plot for the first conformer
        if confId == 0 and plot is True:
            try:
                v = PyMol.MolViewer()
            except ConnectionRefusedError:
                pass
            show_shape(v, mol, confId, shape)
            v.server.do('set transparency=0.5')
            v.server.sphere(list([float(i) for i in mol_coms[confId]]), 0.2, (2, 1, 0), 'COM')
            # principal axes
            v.server.sphere(list([float(i) for i in conf_axes[confId][0, :]]), 0.2, (2, 0, 0), 'AX1')
            v.server.sphere(list([float(i) for i in conf_axes[confId][1, :]]), 0.2, (1, 0, 1), 'AX2')
            v.server.sphere(list([float(i) for i in conf_axes[confId][2, :]]), 0.2, (0, 1, 0), 'AX3')
            for ax in [0, 1, 2]:
                for i, vec in enumerate(vectors[ax]):
                    if vals[ax][i] <= 2:
                        C = (1, 1, 1)
                    else:
                        C = (1, 0, 0)
                    v.server.sphere(list([float(i) for i in vec]), 0.05, C, 'pt')

    return conf_diameters, conf_axes, conf_moments


def get_ellip_diameters(mol, cids, vdwScale=1.0, boxMargin=2.0,
                        spacing=0.2, show=False, plot=False):
    """Fit an ellipsoid to the points within the VDW cloud of a all conformers.

    Keywords:
        mol (RDKIT Molecule) - RDKIT molecule object to calculate ellipsoid for
        cids (list) - list of RDKIT conformer IDs of mol
        vdwScale (float) - Scaling factor for the radius of the atoms to
            determine the base radius used in the encoding
            - grid points inside this sphere carry the maximum occupancy
            default = 1.0 Angstrom
        boxMargin (float) - added margin to grid surrounding molecule
            default=4.0 Angstrom
        spacing (float) - grid spacing - default = 1.0 Angstrom
        show (bool) - show VDW cloud using PYMOL? - default = False
        plot (bool) - show ellipsoid around VDW? - default = False

    Returns:
        conf_diameters (list) - 3 principal diameters of ellipsoids for all
            conformers
        conf_axes (list) - 3 principal axes of mol for all conformers
        conf_axes (list) - 3 principal moments of mol for all conformers
    """
    # over conformers
    conf_diameters = []
    conf_axes = []
    conf_moments = []
    for confId in cids:
        conf = mol.GetConformer(confId)
        box, sideLen, shape = get_molec_shape(mol, conf, confId,
                                              vdwScale=vdwScale,
                                              boxMargin=boxMargin,
                                              spacing=spacing)
        if show is True:
            try:
                v = PyMol.MolViewer()
            except ConnectionRefusedError:
                print("run 'pymol -R' to visualise structures")
                print('visualisation will be skipped')
            show_shape(v, mol, confId, shape)
            v.server.do('set transparency=0.5')

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

        # get inertial properties of conformer
        axes, moments = Chem.ComputePrincipalAxesAndMoments(conf,
                                                            ignoreHs=False)
        conf_axes.append(axes)
        conf_moments.append(moments)
        # find the ellipsoid that envelopes all hit points
        ET = ellipsoid.EllipsoidTool()
        (center, radii, rotation) = ET.getMinVolEllipse(hit_points, .01)

        conf_diameters.append(sorted(np.asarray(radii)*2))

        # only do plot for the first conformer
        if confId == 0 and plot is True:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            # plot points
            # atom_positions = conf.GetPositions()
            # ax.scatter(atom_positions[:, 0], atom_positions[:, 1],
            #            atom_positions[:, 2],
            #            color='k', marker='o', s=100)
            ax.scatter(hit_points[:, 0], hit_points[:, 1], hit_points[:, 2],
                       color='g', marker='x', edgecolor=None,
                       s=50, alpha=0.5)

            # plot ellipsoid
            ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=False)
            ax.set_xlabel("$x$ [$\mathrm{\AA}$]", fontsize=16)
            ax.set_ylabel("$y$ [$\mathrm{\AA}$]", fontsize=16)
            ax.set_zlabel("$z$ [$\mathrm{\AA}$]", fontsize=16)
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
                fig.savefig('temporary_'+str(i)+'_'+str(j)+'.pdf', dpi=720,
                            bbox_inches='tight')

    return conf_diameters, conf_axes, conf_moments


def write_molecule_results(filename, cids, conf_diameters, ratio_1_, ratio_2_):
    """Write results for a molecule to file and PANDAS DataFrame.

    """
    results = pd.DataFrame(columns=['confId', 'diam1', 'diam2',
                                    'diam3', 'ratio_1', 'ratio_2'])
    results['confId'] = cids
    results['diam1'] = [sorted(i)[0] for i in conf_diameters]
    results['diam2'] = [sorted(i)[1] for i in conf_diameters]
    results['diam3'] = [sorted(i)[2] for i in conf_diameters]
    results['ratio_1'] = ratio_1_
    results['ratio_2'] = ratio_2_

    results.to_csv(filename, index=False)
    return results


def calc_molecule_diameter(name, smile, out_dir='./',
                           vdwScale=1.0, boxMargin=4.0, spacing=1.0,
                           MW_thresh=130, show_vdw=False, plot_ellip=False,
                           N_conformers=10):
    """Calculate the diameter of a single molecule.

    Keywords:
        name (str) - name of molecule
        smile (str) - Canononical SMILES of molecule
        out_dir (str) - directory to output molecule files
        vdwScale (float) - Scaling factor for the radius of the atoms to
            determine the base radius used in the encoding
            - grid points inside this sphere carry the maximum occupancy
            default = 1.0 Angstrom
        boxMargin (float) - added margin to grid surrounding molecule
            default=4.0 Angstrom
        spacing (float) - grid spacing - default = 1.0 Angstrom
        MW_thresh (float) - Molecular Weight maximum - default = 130 g/mol
        show_vdw (bool) - show VDW cloud using PYMOL? - default = False
        plot_ellip (bool) - show ellipsoid around VDW? - default = False
        N_conformers (int) - number of conformers to calculate diameter of
            default = 10

    Returns:
        res (DataFrame) - Dataframe of the properties of all conformers

    """
    # check if calculation has already been done
    # need to replace all ' ' and '/' in file names with something else
    out_file = out_dir+name.replace(' ', '_').replace('/', '__')+'_diam_result.csv'
    if os.path.isfile(out_file):
        print('calculation already done.')
        res = pd.read_csv(out_file)
        return res
    # Read SMILES and add Hs
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    MW = Descriptors.MolWt(mol)
    # use a MW threshold (default of 130 g/mol)
    # to avoid unneccesary calculations
    if MW > MW_thresh:
        print(name, 'molecule is too big - skipping....')
        return None
    # try based on RuntimeError from RDKit
    try:
        # 2D to 3D
        # with multiple conformers
        # use a set randomSeed so that running the code multiple times
        # gives the same series of conformers
        cids = Chem.EmbedMultipleConfs(mol=mol, numConfs=N_conformers,
                                       useExpTorsionAnglePrefs=True,
                                       useBasicKnowledge=True,
                                       randomSeed=1000)
        # quick UFF optimize
        for cid in cids:
            Chem.UFFOptimizeMolecule(mol, confId=cid)
    except RuntimeError:
        print('RDKit error. Skipping.')
        res = None
        return res
    _, _, _, ratio_1_, ratio_2_ = get_inertial_prop(mol, cids)
    conf_diameters, conf_axes, conf_moments = get_ellip_diameters(
                                                mol,
                                                cids,
                                                vdwScale=vdwScale,
                                                boxMargin=boxMargin,
                                                spacing=spacing,
                                                show=show_vdw,
                                                plot=plot_ellip)

    res = write_molecule_results(out_file, cids,
                                 conf_diameters, ratio_1_, ratio_2_)
    return res


def calc_molecule_diameters(molecules, out_dir='./',
                            vdwScale=1.0, boxMargin=4.0, spacing=1.0,
                            MW_thresh=130, show_vdw=False, plot_ellip=False,
                            N_conformers=10, rerun=True):
    """Calculate the diameters of a dictionary of molecules.

    Keywords:
        molecules (dict) - {molcule names (str): SMILES (str)}
        out_dir (str) - directory to output molecule files
        vdwScale (float) - Scaling factor for the radius of the atoms to
            determine the base radius used in the encoding
            - grid points inside this sphere carry the maximum occupancy
            default = 1.0 Angstrom
        boxMargin (float) - added margin to grid surrounding molecule
            default=4.0 Angstrom
        spacing (float) - grid spacing - default = 1.0 Angstrom
        MW_thresh (float) - Molecular Weight maximum - default = 130 g/mol
        show_vdw (bool) - show VDW cloud using PYMOL? - default = False
        plot_ellip (bool) - show ellipsoid around VDW? - default = False
        N_conformers (int) - number of conformers to calculate diameter of
            default = 10
        rerun (bool) - rerun previously done molecules? - default = True

    Returns:
        results are saved to files for analysis.

    """
    count = 0
    for name, smile in molecules.items():
        print('molecule:', name, ':', 'SMILES:', smile)
        out_file = out_dir+name.replace(' ', '_').replace('/', '__')+'_diam_result.csv'
        if rerun is False:
            if os.path.isfile(out_file) is True:
                continue
        res = calc_molecule_diameter(name, smile,
                                     out_dir=out_dir,
                                     vdwScale=vdwScale,
                                     boxMargin=boxMargin,
                                     spacing=spacing,
                                     MW_thresh=MW_thresh,
                                     show_vdw=show_vdw,
                                     plot_ellip=plot_ellip,
                                     N_conformers=N_conformers)
        count += 1
        print(count, 'out of', len(molecules), 'done')