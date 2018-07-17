# achive of shadow calculation from 17-07-18

fig, ax = plt.subplots()
confId = 0
conf = mol.GetConformer(confId)
# Chem.CanonicalizeConformer(conf)  # align principal axes with X, Y, Z
atom_positions = conf.GetPositions()
ax.scatter(atom_positions[:, 0], atom_positions[:, 1], c='k', alpha=0.2)
box, sideLen, shape = rdkit_functions.get_molec_shape(mol, conf, confId, vdwScale=vdwScale, 
                                         boxMargin=boxMargin,
                                         spacing=spacing)
ax.scatter(box[1].x+boxMargin, box[1].y+boxMargin, c='b', marker='x')
ax.scatter(box[0].x-boxMargin, box[1].y+boxMargin, c='b', marker='x')
ax.scatter(box[0].x-boxMargin, box[0].y-boxMargin, c='b', marker='x')
ax.scatter(box[1].x+boxMargin, box[0].y-boxMargin, c='b', marker='x')
axes, moments = Chem.ComputePrincipalAxesAndMoments(conf, ignoreHs=False)
COM = mol_coms[confId]
ax.scatter(COM[0], COM[1], c='r')
max_side_len = max([i*2 for i in sideLen])

# define directions based on AX1, AX2, AX3
# use two AX to define plane, and 3rd AX to define direction of vectors
AX1 = np.array([1, 0, 0])
AX2 = np.array([0, 1, 0])
AX3 = np.array([0, 0, 1])
v1, v2 = rdkit_functions.define_plane(AX1, AX2, COM, max_side_len, vec_spacing=vec_spacing)
hit_points = []
for i in range(v1.shape[0]):
    for j in range(v2.shape[0]):
        # point in plane
        a = v1[i] + v2[j]
        ax.scatter(a[0], a[1], c='purple', alpha=0.2)
        # define a series of points that exist on a vector along 3rd AX
        # that goes through 'a'
        dist_a = np.arange(-max_side_len, max_side_len + 0.1, vec_spacing)
        vect_a = np.asarray([(i)*AX3 for i in dist_a]) + a
        
        # determine if the series of points hits VDW cloud
        for p in vect_a:
            pt = rdkit_functions.def_point(*p)
            if shape.GetValPoint(pt) > 2:
                hit_points.append(np.array([pt.x, pt.y, pt.z]))
                break
                
hit_points = np.asarray(hit_points)
ax.scatter(v1[:, 0], v1[:, 1], c='green')
ax.scatter(v2[:, 0], v2[:, 1], c='green')
ax.scatter(hit_points[:, 0], hit_points[:, 1], c='b')

ax.scatter(atom_positions[:, 0], atom_positions[:, 1], c='k', alpha=1)

define_parity_plot_variables(ax, title='',
                             xtitle='$X$',
                             ytitle='$Y$',
                             xlim=(-20, 20),
                             ylim=(-12, 12))