#!/bin/env python

import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space import Visual

# Create cartesian coordinate system

CS = Cartesian()  # if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system

# to visualise the coordinate system basis the module Visual is used

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure
Visual.draw_CS_axes(fig, CS)
mlab.draw(fig)  # start mayavi
mlab.show()

# Now we will do simple rotation of the coordinate system around its [111] axis
CS.rotate_axis_angle(np.ones(3), np.deg2rad(45))  # this is inplace transform
fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure
Visual.draw_CS_axes(fig, CS)
mlab.draw(fig)
mlab.show()  # see the result
