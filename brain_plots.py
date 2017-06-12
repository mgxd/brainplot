# To run, execute these steps:
# srun -p om_interactive -N1 -c2 --mem=8G --pty bash
# source activate mathiasg_vd_env
# module add openmind/xvfb-fix/0.1
# export QT_API=pyqt
# python plot_brain.py -i <input> -o <output> -c <path to conte atlas> -r <path to resting atlas>

# https://github.com/cgoldberg/xvfbwrapper

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """    
    
    import numpy as np
    import math

    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def useZstat(zstat, file_path_name_save, file_path_conte, file_path_name_resting_atlas):
    """Plot and save the image.
    
    Arguments
    ---------
    zstat : string
        Full file path and name to nii to plot.
    
    file_path_name_save : string
        Full file path and name to png output.  Output dir will be created if it doesn't exist.
    
    file_path_conte : string
        Full file path to Conte atlas
        
    file_path_name_resting_atlas : string
    
    Returns
    -------
    None.  Normal error message:  
    pixdim[1,2,3] should be non-zero; setting 0 dims to 1
    plot_brain.py: Fatal IO error: client killed
    
    Example
    -------
    python plot_brain.py -i /groupAnalysis/l2/zstat1_threshold.nii.gz -o /plots/l2test.png -c /git/bp2/32k_ConteAtlas_v2 -r rfMRI_REST1_LR_Atlas.dtseries.nii
    
    MIT OM Specific Tip
    -------------------
    Call this function from a shell script to run headerless BUT requires:
    source activate mathiasg_vd_env
    export QT_API=pyqt
    module add openmind/xvfb-fix/0.1

    #file_path_name=$1
    #file_path_name_save=$2
    #file_path_conte=$3
    #file_path_name_resting_atlas=$4
    python plot_brain.py \
    -i $1 \
    -o $2 \
    -c $3 \
    -r $4
    
    """

    import matplotlib.pyplot as plt
    import os
    from glob import glob
    import numpy as np
    import nibabel as nb
    import nibabel.gifti as gifti

    # Crucial:  xvfb must be imported and started before importing mayavi
    from xvfbwrapper import Xvfb
    vdisplay = Xvfb()
    vdisplay.start()

    # Crashes on this line if run with plain python (not xvfb-run ... python) and if xvfbwrapper is after it.
    from mayavi import mlab
    from tvtk.api import tvtk
    import math

    mlab.options.offscreen = True #offscreen window for rendering

    img = nb.load(file_path_name_resting_atlas)
    #img = nb.load('/Users/MathiasMacbook/Desktop/rfMRI_REST1_LR_Atlas.dtseries.nii')
    mim = img.header.matrix.mims[1]
    #for idx, bm in enumerate(mim.brainModels):
    #    print((idx, bm.indexOffset, bm.brainStructure))
    bm1 = mim.brainModels[0]
    lidx = bm1.vertexIndices.indices
    bm2 = mim.brainModels[1]
    ridx = bm1.surfaceNumberOfVertices + bm2.vertexIndices.indices
    bidx = np.concatenate((lidx, ridx))

    axis = [0, 0, 1]
    theta = np.pi

    inflated = True
    split_brain = True

    surf = gifti.read(file_path_conte + '/Conte69.L.midthickness.32k_fs_LR.surf.gii') 
    verts_L_data = surf.darrays[0].data
    faces_L_data = surf.darrays[1].data

    surf = gifti.read(file_path_conte + '/Conte69.R.midthickness.32k_fs_LR.surf.gii') 
    verts_R_data = surf.darrays[0].data
    faces_R_data = surf.darrays[1].data

    if inflated:
        surf = gifti.read(file_path_conte + '/Conte69.L.inflated.32k_fs_LR.surf.gii')
        verts_L_display = surf.darrays[0].data
        faces_L_display = surf.darrays[1].data
        surf = gifti.read(file_path_conte + '/Conte69.R.inflated.32k_fs_LR.surf.gii')
        verts_R_display = surf.darrays[0].data
        faces_R_display = surf.darrays[1].data
    else:
        verts_L_display = verts_L_data.copy()
        verts_R_display = verts_R_data.copy()
        faces_L_display = faces_L_data.copy()
        faces_R_display = faces_R_data.copy()

    verts_L_display[:, 0] -= max(verts_L_display[:, 0])
    verts_R_display[:, 0] -= min(verts_R_display[:, 0])
    verts_L_display[:, 1] -= (max(verts_L_display[:, 1]) + 1)
    verts_R_display[:, 1] -= (max(verts_R_display[:, 1]) + 1)

    faces = np.vstack((faces_L_display, verts_L_display.shape[0] + faces_R_display))

    if split_brain:
        verts2 = rotation_matrix(axis, theta).dot(verts_R_display.T).T
    else:
        verts_L_display[:, 1] -= np.mean(verts_L_display[:, 1])
        verts_R_display[:, 1] -= np.mean(verts_R_display[:, 1])
        verts2 = verts_R_display

    verts_rot = np.vstack((verts_L_display, verts2))
    verts = np.vstack((verts_L_data, verts_R_data))
    #print verts.shape
    #print faces.shape

    if not os.path.exists(os.path.split(file_path_name_save)[0]):
        os.makedirs(os.path.split(file_path_name_save)[0])    

    print('use zstat')
    img = nb.load(zstat)
    print('loaded img')
    
    threshold = 2.3 # 1000, lower limit
    display_threshold = 6 #8000, upper limit

    data = img.get_data()
    aff = img.affine
    indices = np.round((np.linalg.pinv(aff).dot(np.hstack((verts, 
                                              np.ones((verts.shape[0], 1)))).T))[:3, :].T).astype(int)
    scalars2 = data[indices[:, 0], indices[:, 1], indices[:, 2]]
    scalars2[np.abs(scalars2) < threshold] = 0.
    scalars = np.zeros(verts.shape[0])
    scalars[bidx] = scalars2[bidx]

    negative = positive = False
    if np.any(scalars < 0):
        negative = True
    if np.any(scalars > 0):
        positive = True

    nlabels = 2
    vmin = 0
    vmax = 0
    if negative and positive:
        maxval = max(-scalars.min(), scalars.max())
        if maxval > display_threshold:
            maxval = display_threshold
        vmin = -maxval
        vmax = maxval
        nlabels = 3
        vmin = -display_threshold ######
        vmax = display_threshold ######
    elif negative:
        vmin = scalars.min()
        if vmin < -display_threshold:
            vmin = -display_threshold
        vmax = 0
        vmin = -display_threshold ######
    elif positive:
        vmax = scalars.max()
        if vmax > display_threshold:
            vmax = display_threshold
        vmin = 0
        vmax = display_threshold ######
    #print zstat
    
    dual_split = True

    fig1 = mlab.figure(1, bgcolor=(0, 0, 0))
    mlab.clf()
    mesh = tvtk.PolyData(points=verts_rot, polys=faces)
    mesh.point_data.scalars = scalars
    mesh.point_data.scalars.name = 'scalars'
    surf = mlab.pipeline.surface(mesh, colormap='autumn', vmin=vmin, vmax=vmax)
    if dual_split:
        verts_rot_shifted = verts_rot.copy()
        verts_rot_shifted = rotation_matrix(axis, theta).dot(verts_rot_shifted.T).T
        verts_rot_shifted[:, 2] -= (np.max(verts_rot_shifted[:, 2]) - np.min(verts_rot_shifted[:, 2]))
        verts_rot_shifted[:, 0] -= np.max(verts_rot_shifted[:, 0])
        mesh2 = tvtk.PolyData(points=verts_rot_shifted, polys=faces)
        mesh2.point_data.scalars = scalars
        mesh2.point_data.scalars.name = 'scalars'
        surf2 = mlab.pipeline.surface(mesh2, colormap='autumn', vmin=vmin, vmax=vmax)
    colorbar = mlab.colorbar(surf, nb_labels=nlabels) #, orientation='vertical')
    lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()

    if negative and positive:
        half_index = lut.shape[0] / 2
        index =  int(half_index * threshold / vmax)
        lut[(half_index - index + 1):(half_index + index), :] = 192
        lut[(half_index + index):, :] = 255 * plt.cm.autumn(np.linspace(0, 255, half_index - index).astype(int))
        lut[:(half_index - index), :] = 255 * plt.cm.cool(np.linspace(0, 255, half_index - index).astype(int))
    elif negative:
        index =  int(lut.shape[0] * threshold / abs(vmin))
        lut[(lut.shape[0] - index):, :] = 192
        lut[:(lut.shape[0] - index), :] = 255 * plt.cm.cool(np.linspace(0, 255, lut.shape[0] - index).astype(int))
    elif positive:
        index =  int(lut.shape[0] * threshold / vmax)
        lut[:index, :] = 192
        lut[index:, :] = 255 * plt.cm.autumn(np.linspace(0, 255, lut.shape[0] - index).astype(int))
    lut[:, -1] = 255

    surf.module_manager.scalar_lut_manager.lut.table = lut
    if dual_split:
        surf2.module_manager.scalar_lut_manager.lut.table = lut
    surf.module_manager.scalar_lut_manager.show_scalar_bar = False
    surf.module_manager.scalar_lut_manager.show_legend = False
    surf.module_manager.scalar_lut_manager.label_text_property.font_size = 10
    surf.module_manager.scalar_lut_manager.show_scalar_bar = True
    surf.module_manager.scalar_lut_manager.show_legend = True
    mlab.draw()

    translate = [0, 0, 0]
    if inflated:
        zoom = -700
    else:
        zoom = -600
    if dual_split:
        if inflated:
            translate = [0,   0, -104.01467148]
        else:
            translate = [0,  0, -54.76305802]        
        if inflated:
            zoom = -750
        else:
            zoom = -570
    
    #mlab.view(0, 90.0, zoom, translate)
    mlab.view(9, 90.0)

    print(file_path_name_save)
   
    mlab.savefig(file_path_name_save, figure=fig1, magnification=5)

    vdisplay.stop()    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='brain_plots.py',
                                     description=__doc__)
    parser.add_argument('-i', '--input', required=True, help='input full file path and name .nii.gz')
    parser.add_argument('-o', '--output', required=True, help='output full file pand and name .png')
    parser.add_argument('-c', '--file_path_conte', help='file path to conte atlas folder')
    parser.add_argument('-r', '--file_path_name_resting_atlas', help='resting atlas nii')
    args = parser.parse_args()
  
    useZstat(args.input, args.output, args.file_path_conte, args.file_path_name_resting_atlas)
