from __future__ import division, print_function
import numpy as np

from BDSpace.Coordinates import Cartesian
from BDSpace import Space
from BDSpace.Figure import Figure
from BDSpace.Figure.Sphere import Sphere
from BDSpace.Curve import Curve
from BDSpace.Field import HyperbolicPotentialSphericalConservativeField, SuperposedField

import BDSpaceVis as Visual

from mayavi import mlab

space = Space('Two charged balls')
pos_ball_position = np.array([5.0, 5.0, 0.0], dtype=np.double)
neg_ball_position = np.array([-5.0, -5.0, 0.0], dtype=np.double)

pos_ball_charge = 5.0
pos_ball_radius = 1.0
pos_charged_ball = Sphere(name='Pos Charged Ball', r_outer=pos_ball_radius)
pos_charged_ball.coordinate_system.origin = pos_ball_position

neg_ball_charge = -5.0
neg_ball_radius = 1.0
neg_charged_ball = Sphere(name='Neg Charged Ball', r_outer=neg_ball_radius)
neg_charged_ball.coordinate_system.origin = neg_ball_position

space.add_element(pos_charged_ball)
space.add_element(neg_charged_ball)

pos_electrostatic_field = HyperbolicPotentialSphericalConservativeField(name='Pos Charged ball field',
                                                                        field_type='electrostatic',
                                                                        a=pos_ball_charge, r=pos_ball_radius)
pos_charged_ball.add_element(pos_electrostatic_field)

neg_electrostatic_field = HyperbolicPotentialSphericalConservativeField(name='Neg Charged ball field',
                                                                        field_type='electrostatic',
                                                                        a=neg_ball_charge, r=neg_ball_radius)
neg_charged_ball.add_element(neg_electrostatic_field)

fields_superposition = SuperposedField('Field of two charged balls', [pos_electrostatic_field,
                                                                      neg_electrostatic_field])

space.add_element(fields_superposition)

fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure

pos_charged_ball_vis = Visual.FigureView(fig, pos_charged_ball)
pos_charged_ball_vis.set_color((1, 0, 0))
pos_charged_ball_vis.draw()

neg_charged_ball_vis = Visual.FigureView(fig, neg_charged_ball)
neg_charged_ball_vis.set_color((0, 0, 1))
neg_charged_ball_vis.draw()

fields_superposition_vis = Visual.FieldView(fig, fields_superposition, scalar_field_visible=False)

# pos_ball_field_vis = Visual.FieldView(fig, pos_electrostatic_field)
# neg_ball_field_vis = Visual.FieldView(fig, neg_electrostatic_field)

grid = np.mgrid[-10:10:40j, -10:10:40j, -5:5:10j]

fields_superposition_vis.set_grid(grid)
fields_superposition_vis.set_cs_visible(False)
fields_superposition_vis.draw()

# pos_ball_field_vis.set_grid(grid)
# pos_ball_field_vis.set_cs_visible(False)
# pos_ball_field_vis.draw()

# neg_ball_field_vis.set_grid(grid)
# neg_ball_field_vis.set_cs_visible(False)
# neg_ball_field_vis.draw()

mlab.show()