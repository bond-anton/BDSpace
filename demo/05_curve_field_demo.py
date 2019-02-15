import numpy as np

from BDSpace.Coordinates import Cartesian
from BDSpace import Space
from BDSpace.Curve.Parametric import Helix
from BDSpace.Field import HyperbolicPotentialCurveConservativeField

import BDSpaceVis as Visual
from mayavi import mlab

space = Space('Charged helix')
coordinate_system = Cartesian()
# coordinate_system.rotate_axis_angle(np.ones(3, dtype=np.double), np.deg2rad(45))

left_helix = Helix(name='Left Helix', coordinate_system=coordinate_system,
                   radius=2, pitch=5, start=0, stop=np.pi * 2, right=False)
helix_r = 0.5
print('Helix length:', left_helix.length())
pos_electrostatic_field = HyperbolicPotentialCurveConservativeField(name='Pos Charged Helix field',
                                                                    field_type='electrostatic',
                                                                    curve=left_helix, r=helix_r)
pos_electrostatic_field.a = 1.0
pos_electrostatic_field.curve.precision = 1.0e-2

space.add_element(left_helix)
space.add_element(pos_electrostatic_field)

fig = mlab.figure('CS demo', bgcolor=(0.0, 0.0, 0.0))  # Create the mayavi figure
Visual.draw_coordinate_system_axes(fig, coordinate_system)

left_helix_view = Visual.CurveView(fig=fig, curve=left_helix)
left_helix_view.set_color((0.9, 0.9, 0.9), 0.7)
left_helix_view.set_cs_visible(True)
left_helix_view.draw()

pos_field_vis = Visual.FieldView(fig, pos_electrostatic_field)

grid = np.mgrid[-10:10:40j, -10:10:40j, -5:5:20j]
pos_field_vis.set_grid(grid)
pos_field_vis.set_cs_visible(True)
pos_field_vis.draw()

mlab.show()
