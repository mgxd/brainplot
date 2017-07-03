# vi: set ft=python sts=4 ts=4 sw=4 et:

# To run, execute these steps:
# srun -p om_interactive -N1 -c2 --mem=8G --pty bash
# module add openmind/xvfb-fix/0.1
# python plot_brain.py <input_statistic>

# https://github.com/cgoldberg/xvfbwrapper

import os
from glob import glob
import math
import atexit
os.environ['QT_API'] = 'pyqt'

import matplotlib.pyplot as plt
import nibabel as nb
import nibabel.gifti as gifti
import numpy as np

from xvfbwrapper import Xvfb
vdisplay = Xvfb()
vdisplay.start()
from mayavi import mlab
from tvtk.api import tvtk

@atexit.register
def close_xvfb():
    """ Closes virtual display when exiting; may raise Fatal IO error """
    vdisplay.stop()


def rotation_matrix(axis=[0,0,1], theta=np.pi):
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

def gen_brain_indices(rest_atlas):
    """ Combine left and right indices from resting atlas """
    img = nb.load(rest_atlas)
    try:
        models = img.header.matrix[1].brain_models
        bm1 = next(models)
        lidx = np.array(bm1.vertex_indices)
        bm2 = next(models)
        ridx = bm1.surface_number_of_vertices + np.array(bm2.vertex_indices)
    except TypeError:
        Warning('Using deprecated CIFTI support.')
        mim = img.header.matrix.mims[1]
        bm1 = mim.brainModels[0]
        lidx = bm1.vertexIndices.indices
        bm2 = mim.brainModels[1]
        ridx = bm1.surfaceNumberOfVertices + bm2.vertexIndices.indices
    return np.concatenate((lidx, ridx))

def load_conte_mesh(conte_atlas, inflation):
    """ Load conte69 mesh: load midthickness nomatter what and then other """
    l_atlas = os.path.join(conte_atlas, 'Conte69.L.{}.32k_fs_LR.surf.gii')
    r_atlas = l_atlas.replace('Conte69.L.', 'Conte69.R.')

    l_surf = gifti.read(l_atlas.format('midthickness'))
    verts_L_data, faces_L_data = l_surf.darrays[0].data, l_surf.darrays[1].data

    r_surf = gifti.read(r_atlas.format('midthickness'))
    verts_R_data, faces_R_data = r_surf.darrays[0].data, r_surf.darrays[1].data

    if inflation != 'low':
        if inflation == 'medium':
            inf = 'inflated'
        elif inflation == 'high':
            inf = 'very_inflated'

        l_surf2 = gifti.read(l_atlas.format(inf))
        verts_L_display = l_surf2.darrays[0].data
        faces_L_display = l_surf2.darrays[1].data

        r_surf2 = gifti.read(r_atlas.format(inf))
        verts_R_display = r_surf2.darrays[0].data
        faces_R_display = r_surf2.darrays[1].data

    else:
        verts_L_display = verts_L_data.copy()
        verts_R_display = verts_R_data.copy()
        faces_L_display = faces_L_data.copy()
        faces_R_display = faces_R_data.copy()

    return (verts_L_data, faces_L_data, verts_R_data, faces_R_data,
           verts_L_display, verts_R_display, faces_L_display, faces_R_display)

def plot_stat(args, conte_atlas, rest_atlas):
    """Plot and save the image.

    Arguments
    ---------
    args : Namespace Object
        Contains the following argparse attributes:
            - conte_atlas
            - imagesize
            - in_stat
            - outfile
            - resting_atlas
            - threshold
            - view
            - windowed
            - inflation

    conte_atlas : string
        Path to Conte69 atlas

    rest_atlas : string
        Path to Connectome Resting atlas

    Returns
    -------
    None

    """

    IMAGETYPES = ['nii.gz', 'nii', 'gii']

    mlab.options.offscreen = args.windowed

    bidx = gen_brain_indices(rest_atlas)

    (verts_L_data, faces_L_data, verts_R_data,
    faces_R_data, verts_L_display, verts_R_display,
    faces_L_display, faces_R_display) = load_conte_mesh(conte_atlas,
                                                        args.inflation)

    split_brain = True
    dual_split = False
    if args.view != 'lat':
        split_brain, dual_split = False, False

    verts_L_display[:, 0] -= max(verts_L_display[:, 0])
    verts_R_display[:, 0] -= min(verts_R_display[:, 0])
    verts_L_display[:, 1] -= (max(verts_L_display[:, 1]) + 1)
    verts_R_display[:, 1] -= (max(verts_R_display[:, 1]) + 1)

    faces = np.vstack((faces_R_display, verts_R_display.shape[0] +
                       faces_L_display))

    if split_brain:
        verts2 = rotation_matrix().dot(verts_R_display.T).T
    else:
        verts_L_display[:, 1] -= np.mean(verts_L_display[:, 1])
        verts_R_display[:, 1] -= np.mean(verts_R_display[:, 1])
        verts2 = verts_R_display

    verts_rot = np.vstack((verts_L_display, verts2))
    verts = np.vstack((verts_L_data, verts_R_data))

    if not args.outfile:
        try:
            outfile = os.path.abspath(
               [os.path.basename(args.in_stat).replace(img, 'png')
                for img in IMAGETYPES if args.in_stat.endswith(img)][0])
        except IndexError:
            raise AttributeError('Stat file {} not supported. Supported '
                'extensions: {}'.format(in_stat, ', '.join(IMAGETYPES)))
    else:
        outfile = os.path.abspath(args.outfile)

    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    img = nb.load(args.in_stat) # load statistic

    threshold = args.threshold
    display_threshold = 6

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
        nlabels = 3
        vmin = -display_threshold
        vmax = display_threshold
    elif negative:
        vmin = scalars.min()
        vmin = -display_threshold
    elif positive:
        vmax = scalars.max()
        vmax = display_threshold

    fig1 = mlab.figure(1, bgcolor=(0, 0, 0))
    mlab.clf()
    mesh = tvtk.PolyData(points=verts_rot, polys=faces)
    mesh.point_data.scalars = scalars
    mesh.point_data.scalars.name = 'scalars'
    surf = mlab.pipeline.surface(mesh, colormap='autumn', vmin=vmin, vmax=vmax)
    if dual_split:
        verts_rot_shifted = rotation_matrix().dot(verts_rot.copy().T).T
        verts_rot_shifted[:, 2] -= (np.max(verts_rot_shifted[:, 2]) -
                                    np.min(verts_rot_shifted[:, 2]))
        verts_rot_shifted[:, 0] -= np.max(verts_rot_shifted[:, 0])
        mesh2 = tvtk.PolyData(points=verts_rot_shifted, polys=faces)
        mesh2.point_data.scalars = scalars
        mesh2.point_data.scalars.name = 'scalars'
        surf2 = mlab.pipeline.surface(mesh2, colormap='autumn', vmin=vmin,
                                      vmax=vmax)
    colorbar = mlab.colorbar(surf, nb_labels=nlabels)
    lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()

    if negative and positive:
        half_index = lut.shape[0] // 2
        index =  int(half_index * threshold / vmax)
        lut[(half_index - index):(half_index + index), :] = 192
        lut[(half_index + index):, :] = 255 * plt.cm.autumn(
                          np.linspace(0, 255, half_index - index).astype(int))
        lut[:(half_index - index), :] = 255 * plt.cm.winter(
                          np.linspace(0, 255, half_index - index).astype(int))
    elif negative:
        index =  int(lut.shape[0] * threshold / abs(vmin))
        lut[(lut.shape[0] - index):, :] = 192
        lut[:(lut.shape[0] - index), :] = 255 * plt.cm.winter(
                          np.linspace(0, 255, lut.shape[0] - index).astype(int))
    elif positive:
        index =  int(lut.shape[0] * threshold / vmax)
        lut[:index, :] = 192
        lut[index:, :] = 255 * plt.cm.autumn(
                          np.linspace(0, 255, lut.shape[0] - index).astype(int))
    lut[:, -1] = 255

    surf.module_manager.scalar_lut_manager.lut.table = lut
    if dual_split:
        surf2.module_manager.scalar_lut_manager.lut.table = lut
    surf.module_manager.scalar_lut_manager.show_scalar_bar = False
    surf.module_manager.scalar_lut_manager.show_legend = False
    surf.module_manager.scalar_lut_manager.label_text_property.font_size = 7
    surf.module_manager.scalar_lut_manager.show_scalar_bar = True
    surf.module_manager.scalar_lut_manager.show_legend = True
    mlab.draw()

    # viewing
    if args.inflation == 'low':
        zoom = -600
    else:
        zoom = -700

    if split_brain:
        if args.inflation == 'low':
            zoom = -750
        else:
            zoom = -800
            if dual_split:
                zoom = -950

    if args.view == 'lat':
        x, y = 0, 90
    elif args.view == 'sup':
        x, y = 0, 180
    elif args.view == 'inf':
        x, y = 0, 0

    mlab.view(x, y, zoom)
    mlab.savefig(outfile, figure=fig1, magnification=args.imagesize)

def main():
    conte_atlas = os.path.abspath('Conte69_Atlas')
    #rest_atlas = os.path.abspath('/om/user/mathiasg/scripts/templates/'
    rest_atlas = os.path.abspath('/home/mathiasg/code/datasets/'
                                 'rfMRI_REST1_LR_Atlas.dtseries.nii')

    import argparse
    def existing_file(filename):
        filename = os.path.abspath(filename)
        if not os.path.exists(filename):
            raise argparse.ArgumentTypeError('{} is not an existing '
                                             'file!'.format(filename))
        return filename

    parser = argparse.ArgumentParser(prog='brain_plots.py',
                                     description=__doc__)
    parser.add_argument('in_stat', type=existing_file,
                        help='input full file path')
    parser.add_argument('-o', '--outfile',
                        help='output full file path')
    parser.add_argument('-c', '--conte_atlas',
                        help='Conte69 32k mesh atlas folder')
    parser.add_argument('-r', '--resting_atlas',
                        help='Connectome resting atlas')
    parser.add_argument('-t', '--threshold', type=float, default=2.3,
                        help='set threshold (default 2.3)')
    parser.add_argument('-s', '--imagesize', type=int, default=2,
                        choices=range(1,6),
                        help='set image size; 1-smallest, 5-largest')
    parser.add_argument('-v', '--view', default='lat',
                        choices=['lat', 'sup', 'inf'], # porque no los tres?
                        help='view of brain: lateral, superior, inferior')
    parser.add_argument('-w', '--windowed', action='store_false',
                        help='run windowed')
    parser.add_argument('-i', '--inflation', default='medium',
                        choices=['low', 'medium', 'high'],
                        help='inflation level of generated brain mesh: low, '
                        'medium, or high')
    args = parser.parse_args()

    if args.conte_atlas:
        conte_atlas = os.path.abspath(args.conte_atlas)

    if args.resting_atlas:
        rest_atlas = os.path.abspath(args.resting_atlas)

    plot_stat(args, conte_atlas, rest_atlas)

if __name__ == '__main__':
    main()
