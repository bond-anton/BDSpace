#!/bin/env python

import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space import Visual

# Create cartesian coordinate system

CS = Cartesian()  # if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system

# to visualise the coordinate system basis the module Visual is used

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure
#Visual.draw_CS_axes(fig, CS)
cube_surface, arrows, labels = Visual.draw_CS_box(fig, CS)
CS.rotate_axis_angle(np.ones(3), np.deg2rad(45))
cube_surface, arrows, labels = Visual.update_CS_box(CS, cube_surface, arrows, labels)
mlab.show()  # start mayavi

# Now we will do simple rotation of the coordinate system around its [111] axis
CS.rotate_axis_angle(np.ones(3), np.deg2rad(45))  # this is inplace transform
fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure
arrows, labels = Visual.draw_CS_axes(fig, CS, draw_labels=True)
CS.rotate_axis_angle(np.ones(3), np.deg2rad(45))
CS.rotate_axis_angle(np.ones(3), np.deg2rad(90))
Visual.update_CS_axes(CS, arrows, labels, offset=0)

#mlab.draw(fig)
mlab.show()  # see the result
