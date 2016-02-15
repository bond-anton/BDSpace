#!/bin/env python

import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space.Coordinates import transforms as gt
from Space import Visual

# Create cartesian coordinate system

CS = Cartesian()  # if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system

# to visualise the coordinate system basis the module Visual is used

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure

@mlab.animate(delay=100)
def anim():
    #f = mlab.gcf()
    cs_box, arrows, labels = Visual.draw_CS_box(fig, CS)
    while True:
        CS.rotate_axis_angle(np.ones(3), np.deg2rad(5))  # this is inplace transform
        cs_box, arrows, labels = Visual.update_CS_box(CS, cs_box, arrows, labels)
        yield

anim()
mlab.show()
