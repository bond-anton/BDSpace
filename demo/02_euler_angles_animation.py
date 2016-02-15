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
cs_box = Visual.draw_CS_box(fig, CS)
CS.rotate_axis_angle(np.ones(3), np.deg2rad(45))  # this is inplace transform
euler_color = gt.euler_color(CS.euler_angles)
cs_box.actor.property.edge_visibility = 1
cs_box.actor.property.color = euler_color
cs_box.actor.property.edge_color = euler_color

@mlab.animate(delay=100)
def anim():
    #f = mlab.gcf()
    while True:
        CS.rotate_axis_angle(np.ones(3), np.deg2rad(5))  # this is inplace transform
        euler_color = gt.euler_color(CS.euler_angles)
        cs_box.actor.property.edge_visibility = 1
        cs_box.actor.property.color = euler_color
        cs_box.actor.property.edge_color = euler_color
        yield

anim()
mlab.show()
