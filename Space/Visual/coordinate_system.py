from __future__ import division

import numpy as np
from Space.Coordinates import transforms as gt
from mayavi import mlab
from tvtk.api import tvtk

from Space.Coordinates.coordinates import Cartesian
from Space.Figures import figures


def CS_arrows(CS, offset=0, scale=1.0):
    points = []
    lengths = []
    for i in range(3):
        points.append(CS.origin + scale * CS.basis[:, i] * offset)
        lengths.append(CS.basis[i, :] * scale)
    points = np.array(points)
    lengths = np.array(lengths)
    return points, lengths


def draw_CS_axes(fig, CS, offset=0, scale=1.0, draw_labels=True):
    points, lengths = CS_arrows(CS, offset=offset, scale=scale)
    mlab.figure(fig, bgcolor=fig.scene.background)
    arrows = mlab.quiver3d(points[:, 0], points[:, 1], points[:, 2],
                           lengths[0, :], lengths[1, :], lengths[2, :],
                           scalars=np.array([3, 2, 1]), mode='arrow')
    arrows.glyph.color_mode = 'color_by_scalar'
    arrows.glyph.glyph.scale_factor = scale
    data = arrows.parent.parent
    data.name = CS.name
    glyph_scale = arrows.glyph.glyph.scale_factor * 1.1
    label_col = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    labels = []
    if draw_labels:
        for i in range(3):
            labels.append(mlab.text3d(points[i, 0] + glyph_scale * CS.basis[0, i],
                                      points[i, 1] + glyph_scale * CS.basis[1, i],
                                      points[i, 2] + glyph_scale * CS.basis[2, i],
                                      CS.labels[i], color=label_col[i], scale=0.1*scale))
    return arrows, labels


def draw_CS_box(fig, CS, offset=0.5, scale=1.0, draw_axes=True, draw_labels=True):
    mlab.figure(fig, bgcolor=fig.scene.background)
    cube_points = figures.generate_cube_points(scale, scale, scale,
                                               origin=np.array([scale/2, scale/2, scale/2]))
    cube = tvtk.StructuredGrid(dimensions=(2, 2, 2))
    cube.points = CS.to_global(cube_points)
    euler_color = gt.euler_color(CS.euler_angles)
    cube_surface = mlab.pipeline.surface(cube, color=euler_color)
    cube_surface.actor.property.edge_visibility = 1
    cube_surface.actor.property.edge_color = euler_color
    #for prop in cube_surface.parent.parent.all_trait_names():
    #    print prop
    #cube_surface.parent.parent.data = cube
    #print cube
    #print cube_surface.mlab_source
    arrows, labels = None, None
    if draw_axes:
        arrows, labels = draw_CS_axes(fig, CS, offset=offset, scale=scale, draw_labels=draw_labels)
    return cube_surface, arrows, labels


def update_CS_axes(CS, arrows, labels, offset=0, scale=1.0):
    points, lengths = CS_arrows(CS, offset=offset, scale=scale)
    data = arrows.parent.parent
    data.mlab_source.points = points
    data.mlab_source.u = lengths[0, :]
    data.mlab_source.v = lengths[1, :]
    data.mlab_source.w = lengths[2, :]
    glyph_scale = arrows.glyph.glyph.scale_factor * 1.1
    for i in range(len(labels)):
        labels[i].position = np.array([points[i, 0] + glyph_scale * CS.basis[0, i],
                                       points[i, 1] + glyph_scale * CS.basis[1, i],
                                       points[i, 2] + glyph_scale * CS.basis[2, i]])
        labels[i].scale = np.ones(3) * 0.1 * scale
    return arrows, labels


def update_CS_box(CS, cube_surface, arrows, labels, offset=0.5, scale=1.0):
    cube_points = figures.generate_cube_points(scale, scale, scale,
                                               origin=np.array([scale/2, scale/2, scale/2]))
    cube = tvtk.StructuredGrid(dimensions=(2, 2, 2))
    cube.points = CS.to_global(cube_points)
    euler_color = gt.euler_color(CS.euler_angles)
    cube_surface.parent.parent.data = cube
    cube_surface.actor.property.edge_visibility = 1
    cube_surface.actor.property.edge_color = euler_color
    cube_surface.actor.property.color = euler_color
    if arrows is None:
        return cube_surface, arrows, labels
    else:
        arrows, labels = update_CS_axes(CS, arrows, labels, offset=offset, scale=scale)
    return cube_surface, arrows, labels

