'''
Created on Mar 24, 2014

@author: Anton Bondarenko
'''

from __future__ import division

import numpy as np
from Space.Coordinates import transforms as gt
from mayavi import mlab
from tvtk.api import tvtk

from Space.Coordinates.coordinates import Cartesian
from Space.Figures import figures

Scale = 100.0

def draw_CS_axes(fig, CS, scale=1.0):
    points = np.array([CS.origin, CS.origin, CS.origin])
    mlab.figure(fig, bgcolor=fig.scene.background)
    points_v = mlab.quiver3d(points[:,0], points[:,1], points[:,2],
                          CS.basis[0,:]*scale, CS.basis[1,:]*scale, CS.basis[2,:]*scale,
                          scalars=np.array([3,2,1]), mode='arrow')
    points_v.glyph.color_mode = 'color_by_scalar'
    points_v.glyph.glyph.scale_factor = scale
    data = points_v.parent.parent
    data.name = CS.name
    glyph_scale = points_v.glyph.glyph.scale_factor * 1.1
    label_col = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for i in range(3):
        mlab.text3d(points[i, 0] + glyph_scale * CS.basis[0, i],
                    points[i, 1] + glyph_scale * CS.basis[1, i],
                    points[i, 2] + glyph_scale * CS.basis[2, i],
                    CS.labels[i], color=label_col[i], scale=0.1*scale)

def draw_CS_box(fig, CS, scale=1.0, draw_axes=True):
    mlab.figure(fig, bgcolor=fig.scene.background)
    cube_points = figures.generate_cube_points(scale, scale, scale,
                                               origin=np.array([scale/2, scale/2, scale/2]))
    cube = tvtk.StructuredGrid(dimensions=(2, 2, 2))
    cube.points = CS.to_global(cube_points)
    euler_color = gt.euler_color(CS.euler_angles)        
    cube_s = mlab.pipeline.surface(cube, color=euler_color)
    cube_s.actor.property.edge_visibility = 1
    cube_s.actor.property.edge_color = euler_color
    if not draw_axes:
        return
    points = [CS.origin]
    for i in range(3):
        points.append(scale/2 * CS.basis[:, i] + CS.origin)
    points = np.array(points)
    points_v = mlab.quiver3d(points[1:,0], points[1:,1], points[1:,2],
                    CS.basis[0,:]*scale, CS.basis[1,:]*scale, CS.basis[2,:]*scale,
                    scalars=np.array([3,2,1]), mode='arrow')
    points_v.glyph.color_mode = 'color_by_scalar'
    points_v.glyph.glyph.scale_factor = scale/2
    glyph_scale = points_v.glyph.glyph.scale_factor * 1.1
    label_col = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for i in range(3):
        mlab.text3d(points[i+1, 0] + glyph_scale * CS.basis[0, i],
                    points[i+1, 1] + glyph_scale * CS.basis[1, i],
                    points[i+1, 2] + glyph_scale * CS.basis[2, i],
                    CS.labels[i], color=label_col[i], scale=0.1*scale)

def draw_polepiece(fig, CS, angle=55, hole=0.1, R_min=0.2, height=1, draw_axes=False):
    CS_scaled = Cartesian(origin=np.copy(CS.origin)*Scale, basis=np.copy(CS.basis))
    #CS_scaled.origin *= Scale
    #hole *= Scale
    #R_min *= Scale
    #height *=Scale
    dims = (2, 64, 2)
    phi = np.linspace(0, 2 * np.pi, dims[1], endpoint=True)
    z = np.linspace(0, height, dims[2])
    points = figures.generate_polepiece(phi, z, gt.deg_to_rad(angle), hole, R_min) * Scale
    sgrid = tvtk.StructuredGrid(dimensions=(dims[1], dims[0], dims[2]))
    sgrid.points = CS_scaled.to_global(points)
    mlab.figure(fig, bgcolor=fig.scene.background)
    polepiece = mlab.pipeline.add_dataset(sgrid, 'Polepiece')
    mlab.pipeline.surface(polepiece)
    if draw_axes:
        scale = hole
        draw_CS_axes(fig, CS_scaled, scale * Scale)

def draw_sample(fig, CS, size, draw_axes=False, draw_holder=False):
    CS_scaled = Cartesian(origin=np.copy(CS.origin)*Scale, basis=np.copy(CS.basis))
    #CS_scaled = CS
    #CS_scaled.origin *= Scale
    #size *= Scale
    sample_points = figures.generate_cube_points(size[0], size[1], size[2],
                                                 origin=np.array([size[0]/2, size[1]/2, size[2]]))
    sgrid = tvtk.StructuredGrid(dimensions=(2, 2, 2))
    sgrid.points = CS_scaled.to_global(sample_points * Scale)
    mlab.figure(fig, bgcolor=fig.scene.background)
    sample = mlab.pipeline.add_dataset(sgrid, 'Sample')
    mlab.pipeline.surface(sample)
    if draw_axes:
        scale = min(size[0], size[1]) / 2
        draw_CS_axes(fig, CS_scaled, scale * Scale)

def draw_detector(fig, CS, screen_width=1, vh_ratio=0.75,
                  pattern_center=[0.5, 0.5], draw_axes=False):
    CS_scaled = Cartesian(origin=np.copy(CS.origin)*Scale, basis=np.copy(CS.basis))
    #CS_scaled = CS
    #CS_scaled.origin *= Scale
    #screen_width *= Scale
    d_h = 0.02 * screen_width
    pc_x = screen_width * pattern_center[0]
    pc_y = screen_width * pattern_center[1] * vh_ratio
    detector_points = figures.generate_cube_points(screen_width, screen_width * vh_ratio, d_h,
                                                   origin=[pc_x, pc_y, d_h])
    sgrid = tvtk.StructuredGrid(dimensions=(2, 2, 2))
    sgrid.points = CS_scaled.to_global(detector_points * Scale)
    mlab.figure(fig, bgcolor=fig.scene.background)
    detector = mlab.pipeline.add_dataset(sgrid, 'Detector Screen')
    detector_surface = mlab.pipeline.surface(detector)
    detector_surface.actor.property.color=(1, 1, 0)
    if draw_axes:
        scale = min(screen_width, screen_width * vh_ratio) / 4
        draw_CS_axes(fig, CS_scaled, scale * Scale)

def draw_beamline(fig, CS_pp, CS_sample, CS_det, pp_pos=np.array([0, 0, 0]), sample_pos=np.array([0, 0, 0])):
    '''
    '''
    CS_pp_scaled = Cartesian(origin=np.copy(CS_pp.origin)*Scale, basis=np.copy(CS_pp.basis))
    CS_sample_scaled = Cartesian(origin=np.copy(CS_sample.origin)*Scale, basis=np.copy(CS_sample.basis))
    CS_det_scaled = Cartesian(origin=np.copy(CS_det.origin)*Scale, basis=np.copy(CS_det.basis))
    
    points = np.vstack((CS_pp_scaled.to_global(pp_pos), CS_sample_scaled.to_global(sample_pos)))
    mlab.figure(fig, bgcolor=fig.scene.background)
    mlab.plot3d(points[:, 0], points[:, 1], points[:, 2], color=(0, 1, 0), line_width=5, name='E-Beam')
    incident_beam = CS_sample_scaled.to_local(points[0] - points[1])
    
    reflected_beam_polar = CS_sample_scaled.to_polar(incident_beam)
    reflected_beam_polar[2] += np.pi
    reflected_beam_polar[2] -= 2*np.pi*(reflected_beam_polar[2] // (2*np.pi))
    reflected_beam = CS_sample_scaled.to_cartesian(reflected_beam_polar) + sample_pos
    reflected_beam = np.vstack((CS_sample_scaled.to_global(sample_pos), CS_sample_scaled.to_global(reflected_beam)))
    mlab.plot3d(reflected_beam[:, 0], reflected_beam[:, 1], reflected_beam[:, 2], color=(0, 1, 0), line_width=5, name='Reflected beam')
    
    dd_0_points = np.vstack((CS_det_scaled.origin, CS_sample_scaled.origin))
    mlab.plot3d(dd_0_points[:, 0], dd_0_points[:, 1], dd_0_points[:, 2], color=(1, 1, 1), line_width=5, name='DD scan center')
    
    dd_vec = CS_det_scaled.to_local(sample_pos, CS_sample_scaled)
    dd_len = dd_vec[2]
    dd_vec = -dd_len * CS_det_scaled.basis[:, 2] + CS_sample_scaled.to_global(sample_pos)
    dd_points_current = np.vstack((dd_vec, CS_sample_scaled.to_global(sample_pos)))
    mlab.plot3d(dd_points_current[:, 0], dd_points_current[:, 1], dd_points_current[:, 2], color=(0, 0, 1), line_width=5, name='DD current')
    
    dd_polar = CS_sample_scaled.to_local(dd_vec)
    dd_polar = CS_sample_scaled.to_polar(dd_polar)
    norm_line = sample_pos + np.array([0, 0, dd_len/np.cos(dd_polar[1])])
    norm_line = np.vstack((CS_sample_scaled.to_global(sample_pos), CS_sample_scaled.to_global(norm_line)))
    mlab.plot3d(norm_line[:, 0], norm_line[:, 1], norm_line[:, 2], color=(1, 1, 1), line_width=5, name='Perpendicular')
