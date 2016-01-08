#om_py_plots_surfs.py

#generates a variety of types of brains on openmind

#1-8-16: Use the following command to run
####  export QT_API=pyqt

# To run windowless:

####xvfb-run --server-args="-screen 0 1024x768x24" python om_py_plot_surfs.py

import matplotlib.pyplot as plt
import os
from glob import glob
import numpy as np
import nibabel as nb
import nibabel.gifti as gifti
from mayavi import mlab
from tvtk.api import tvtk
import math

# Change these to fit your data

base = '/om/project/voice/processedData/l1analysis/l1output_20151202' ##

model = 'model200' ##

tasks = ["task001","task002","task003","task005","task006","task007"]

#where the images will save
output_dir = '/om/project/voice/processedData/plots/speech_baseline/cool' ##



img = nb.load('/om/user/mathiasg/rfMRI_REST1_LR_Atlas.dtseries.nii')
mim = img.header.matrix.mims[1]
bm1 = mim.brainModels[0]
lidx = bm1.vertexIndices.indices
bm2 = mim.brainModels[1]
ridx = bm1.surfaceNumberOfVertices + bm2.vertexIndices.indices
bidx = np.concatenate((lidx, ridx))


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

axis = [0, 0, 1]
theta = np.pi

inflated = True ## Semi-inflated and Inflated
split_brain = True ## Split or whole brain

surf = gifti.read('32k_ConteAtlas_v2/Conte69.L.midthickness.32k_fs_LR.surf.gii') #inflated.32k_fs_LR.surf.gii')

verts_L_data = surf.darrays[0].data
faces_L_data = surf.darrays[1].data

surf = gifti.read('32k_ConteAtlas_v2/Conte69.R.midthickness.32k_fs_LR.surf.gii') #inflated.32k_fs_LR.surf.gii')
verts_R_data = surf.darrays[0].data
faces_R_data = surf.darrays[1].data

if inflated:
    surf = gifti.read('32k_ConteAtlas_v2/Conte69.L.inflated.32k_fs_LR.surf.gii')
    verts_L_display = surf.darrays[0].data
    faces_L_display = surf.darrays[1].data
    surf = gifti.read('32k_ConteAtlas_v2/Conte69.R.inflated.32k_fs_LR.surf.gii')
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

# Loads zstat for subject
def useZstat(zstat,task,subj):
    img = nb.load(zstat)
    threshold = 2.3 # activation threshold
    display_threshold = 6 # max/min for uniformity

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
    elif negative:
        vmin = scalars.min()
        if vmin < -display_threshold:
            vmin = -display_threshold
        vmax = 0
    elif positive:
        vmax = scalars.max()
        if vmax > display_threshold:
            vmax = display_threshold
        vmin = 0
    print zstat
    
    
    dual_split = True ## further splitting

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
    surf.module_manager.scalar_lut_manager.show_scalar_bar = True #displays bar
    surf.module_manager.scalar_lut_manager.show_legend = True #displays legend
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
    
    os.chdir(output_dir)

    figname = '%s_%s.png' % (subj,task) #savename

    #figname = '%s_%s_%s.png' % (subj,task,num)
    
    mlab.savefig(figname, figure=fig1, magnification=5)

# runs through all subjects
# export QT_API=pyqt

#for task in tasks:
#    print '-----%s-----' % (task)
#    taskdir = os.path.join(base,task)
#    contrasts = next(os.walk(taskdir))[1]
#    
#    for contrast in contrasts:
#        if 'raw' in contrast:
#            continue
#        for num in ['1','2','3','4']:
#            contrastnum = 'contrast_%s' % (num)        
#            subpath = (os.path.join(taskdir,contrast,'stats',contrastnum,'zstat1_threshold.nii.gz'))
#            print subpath
#            useZstat(subpath)
#            visualize(task,contrast,num)

for task in tasks:
    print '-------%s-------' % (task.upper())
    taskdir = os.path.join(base,model,task,model,task)
    #subjs = next(os.walk(taskdir))[1] #takes all dirs in next directory
    subjs = ['voice999'] #specific subjs

    for subj in subjs:
        if subj[:5] != 'voice': #eliminate non subject dir 
            continue
        #for num in ['01','02']: #contrast
        #zNum = 'zstat_%s' % (num)        
        subpath = os.path.join(taskdir,subj,'zstats','mni','zstat01.nii.gz') ## depends on contrast
        useZstat(subpath,task,subj)
