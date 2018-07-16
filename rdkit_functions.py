#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for RDKIT usage.

Author: Andrew Tarzia

Date Created: 16 Jul 2018

License:


"""

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Descriptors3D import NPR1, NPR2, PMI1, PMI2, PMI3
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Geometry import rdGeometry
from rdkit.Chem import PyMol
from rdkit import Geometry
import tempfile


def draw_smiles_to_svg(smiles, filename):
    """Draw a single molecule to an SVG file.
    
    """
    mol = Chem.MolFromSmiles(smiles)
    Chem.Compute2DCoords(mol)
    Draw.MolToFile(mol, filename, fitImage=True, imageType='svg')


def read_mol_txt_file(filename):
    """Function to read molecule SMILES and information from txt file.
    
    """
    data = pd.read_table(filename, delimiter=':')
    molecules = {}
    diameters = {}
    for i, row in data.iterrows():
#         try:
#             name, smile, radius = line.rstrip().split(':')
#         except ValueError:
#             print(line, 'had : in there twice, fix this naming or SMILE')
#             print('skipped')
        name = row['molecule']
        smile = row['smile']
        diameter = row['diameter']
        molecules[name] = smile
        diameters[name] = diameter
    return data, molecules, diameters


def produce_quick_fig_mol(molecules, filename):
    """Produce a quick/dirty figure showing all the 2D coordinates of molecules in the data set.
    
    """
    DrawingOptions.bondLineWidth = 1.8
    DrawingOptions.atomLabelFontSize = 16
    mols = [Chem.MolFromSmiles(x) for x in molecules.values()]
    for m in mols: tmp = Chem.Compute2DCoords(m)
    # Draw.MolToFile(mols[0], output_dir+'mol1.png')
    img=Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200), legends=[x for x in molecules.keys()])
    img.save(filename)
    

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
    
    Code from: https://iwatobipen.wordpress.com/2016/08/16/scoring-3d-diversity-using-rdkit-rdkit/
    
    """
    coms = []
    numatoms = mol.GetNumAtoms()
    for confId in range(len(cids)):
        # print('conf:', confId)
        # print('number of atoms:', numatoms)
        conf = mol.GetConformer(confId)    
        coords =  np.array([list(conf.GetAtomPosition(atmidx)) for atmidx in range(numatoms)])
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
        viewer.ShowMol(mol, confId=cid, name='Conf-%d'%cid, showOnly=False)
        

def show_axes(mol, confId, mol_coms, conf_axes):
    """Show the principal moments of inertia on a single conformation.
    
    """
    # show the principle axes on one conformer
    try:
        v = PyMol.MolViewer()
    except ConnectionRefusedError:
        pass
    v.server.do('reinitialize')
    v.server.do('delete all')

    #v.DeleteAll()
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
    
    
def get_molec_shape(mol, conf, confId, vdwScale=1.0, boxMargin=2.0, spacing=0.2):
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
    Chem.EncodeShape(mol, shape, confId=confId, ignoreHs=False, vdwScale=vdwScale)
    return sideLen, shape


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


def get_boundary_idx(values):
    """Get index of shape grid boundaries along axis.
    
    """
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
        # try:
        #     Chem.CanonicalizeConformer(conf)
        # except RuntimeError:
        #     pass
        sideLen, shape = get_molec_shape(mol, conf, confId, vdwScale=vdwScale, 
                                         boxMargin=boxMargin,
                                         spacing=spacing)
        if show is True:
            try:
                v = PyMol.MolViewer()
            except ConnectionRefusedError:
                pass
            show_shape(v, mol, confId, shape)
            v.server.do('set transparency=0.5')

        # get extent of shape along principle axes
        axes, moments = Chem.ComputePrincipalAxesAndMoments(conf, ignoreHs=False)
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
            diameter = distances[ind_1] + distances[ind_2]
            diameters.append(diameter)
        conf_diameters.append(diameters)

        # # get ellipsoid
        # # equations from:
        # # Elliptic fit of objects in two and three dimensions
        # # by moment of inertia optimization 
        # A_radius = np.sqrt((5/2) * (mid_PMI + lge_PMI - sml_PMI))
        # B_radius = np.sqrt((5/2) * (sml_PMI + lge_PMI - mid_PMI))
        # C_radius = np.sqrt((5/2) * (sml_PMI + mid_PMI - lge_PMI))
        # print(sml_PMI, mid_PMI, lge_PMI)
        # print(2*A_radius, 2*B_radius, 2*C_radius)
        
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

