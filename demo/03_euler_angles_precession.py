#!/bin/env python

import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space.Coordinates import transforms as gt
from Space import Visual

# Create cartesian coordinate system

# if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system
CS_1 = Cartesian(origin=np.array([0, 0, 0]))
CS_2 = Cartesian(origin=np.array([0, 0, 0]))
step_prec = 1.0  # in degrees
step_rot = 10.0  # in degrees
direction_prec = 1
direction_rot = -1
Phi = 30 # tilt in degrees
CS_1.set_eulers(CS_1.euler_angles + np.array([0, np.deg2rad(Phi), 0]))  # tilt the CS
# to visualise the coordinate system basis the module Visual is used

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure

@mlab.animate(delay=10)
def anim():
    cs_box_1, arrows_1, labels_1 = Visual.draw_CS_box(fig, CS_1)
    arrows_2, labels_2 = Visual.draw_CS_axes(fig, CS_2, scale=2)
    while True:
        delta_eulers = np.array([direction_prec * np.deg2rad(step_prec), 0, direction_rot * np.deg2rad(step_rot)])
        CS_1.set_eulers(CS_1.euler_angles + delta_eulers)
        cs_box_1, arrows_1, labels_1 = Visual.update_CS_box(CS_1, cs_box_1, arrows_1, labels_1)
        yield

anim()
mlab.show()
