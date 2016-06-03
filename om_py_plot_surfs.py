# xvfb-run --server-args="-screen 0 1024x768x24" python om_py_plot_surfs.py
import sip
sip.setapi('QDate', 2)
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QTime', 2)
sip.setapi('QUrl', 2)
sip.setapi('QVariant', 2)
sip.setapi('QDateTime', 2)
import matplotlib.pyplot as plt
import os
import numpy as np
import nibabel as nb
import nibabel.gifti as gifti
from mayavi import mlab
from tvtk.api import tvtk
import math
import argparse

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
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


def make_plot(stat, task, contrast, num, outdir, inflated,
              split_brain, dual_split, threshold,
              display_threshold, atlas_dir):
    # load from here until can include
    try:
        img = nb.load('/om/user/mathiasg/rfMRI_REST1_LR_Atlas.dtseries.nii')
    except:
        print('File missing - message mathiasg@mit.edu')
        raise FileNotFoundError
    try:
        bm1 = mim.brain_models[0]
        lidx = bm1.vertex_indices.indices
        bm2 = mim.brain_models[1]
        ridx = bm1.surface_number_of_vertices + bm2.vertex_indices.indices
    except AttributeError: #older citfi version
        bm1 = mim.brainModels[0]
        lidx = bm1.vertexIndices.indices
        bm2 = mim.brainModels[1]
        ridx = bm1.surfaceNumberOfVertices + bm2.vertexIndices.indices
    bidx = np.concatenate((lidx, ridx))
    axis = [0, 0, 1]
    theta = np.pi
    try:
        surf = gifti.read(os.path.join(atlas_dir,'Conte69.L.midthickness.32k_fs_LR.surf.gii'))
    except:
        print('Atlas not found - pass in path with flag -a')
        raise FileNotFoundError
    verts_L_data = surf.darrays[0].data
    faces_L_data = surf.darrays[1].data
    surf = gifti.read(os.path.join(atlas_dir,'Conte69.R.midthickness.32k_fs_LR.surf.gii'))
    verts_R_data = surf.darrays[0].data
    faces_R_data = surf.darrays[1].data
    if inflated:
        surf = gifti.read(os.path.join(atlas_dir,'Conte69.L.inflated.32k_fs_LR.surf.gii'))
        verts_L_display = surf.darrays[0].data
        faces_L_display = surf.darrays[1].data
        surf = gifti.read(os.path.join(atlas_dir,'Conte69.R.inflated.32k_fs_LR.surf.gii'))
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
    #load stat
    img = nb.load(stat)
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
        vmin = -display_threshold
        vmax = display_threshold
    elif negative:
        vmin = scalars.min()
        if vmin < -display_threshold:
            vmin = -display_threshold
        vmax = 0
        vmin = -display_threshold
    elif positive:
        vmax = scalars.max()
        if vmax > display_threshold:
            vmax = display_threshold
        vmin = 0
        vmax = display_threshold
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
    mlab.view(0, 90.0, zoom, translate)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #os.chdir(outdir)
    outname = '%s-%s-%s.png' % (task,contrast,num)
    mlab.savefig(os.path.join(outdir,outname), figure=fig1, magnification=5)

def plot_stats(base, tasks, outdir, atlas_dir, inflated=True,
               split_brain=True, dual_split=True, threshold=2.3,
               display_threshold=6):
    for task in tasks:
        print('-----%s-----'%task)
        taskdir = os.path.join(base,task)
        for contrast in os.listdir(taskdir):
            cons = os.path.join(taskdir,contrast,'stats')
            for x in os.listdir(cons):
                subpath = os.path.join(cons,x,'zstat1.nii.gz')
                print("Converting:\n" + subpath)
                make_plot(subpath, task, contrast, x[-1], outdir,
                          inflated, split_brain, dual_split, 
                          threshold, display_threshold, atlas_dir)
                print("Finished!\n")

if __name__ == '__main__':
    docstr = '\n'.join((__doc__,
"""
           Example:
           python om_py_plot_surfs.py -d mydata/zstats -t task001 task002 task003 -th 2.5
s3
"""))
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data_dir',
                        dest='data_dir',
                        required=True,
                        help='''location of the data to plot''')
    parser.add_argument('-t', '--tasks', dest='tasks', required=True,
                        type=str, nargs='+', help='''list of tasks to get
                        contrasts''')
    parser.add_argument('-o', '--outdir', dest='outputdir',
                        default=os.getcwd(),
                        help='''output directory for resulting images''')
    parser.add_argument('-a', '--atlas', dest='atlas_dir',
                        default=os.path.abspath('32k_ConteAtlas_v2'),
                        help='''brain atlas directory, default cwd''')
    parser.add_argument('-i', '--inflated', dest='inflated',
                        default=True, action='store_false',
                        help='''disable inflated brain image''')
    parser.add_argument('-s', '--split', dest='split_brain',
                        default=True, action='store_false',
                        help='''disable split brain image''')
    parser.add_argument('-ss', '--duosplit', dest='dual_split',
                        default=True, action='store_false',
                        help='''disable dualsplit brain image''')
    parser.add_argument('-th', '--threshold', dest='threshold',
                        type=float, default=2.3,
                        help='''set threshold value (default=2.3) - must
                        be a float''')
    parser.add_argument('-dt', '--displaythresh', dest='display_threshold',
                        type=int, default=6,
                        help='''set min/max for thresholded values - must
                        be an int''')
    args = parser.parse_args()
    plot_stats(args.data_dir, args.tasks,
               os.path.abspath(args.outputdir),
               atlas_dir=os.path.abspath('32k_ConteAtlas_v2')
               inflated=args.inflated,
               split_brain=args.split_brain,
               dual_split=args.dual_split,
               threshold=args.threshold,
               display_threshold=args.display_threshold)
