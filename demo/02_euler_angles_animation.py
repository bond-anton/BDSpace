#!/bin/env python

import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space.Coordinates import transforms as gt
from Space import Visual

# Create cartesian coordinate system

# if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system
CS_1 = Cartesian(origin=np.array([0, 0, 0]))
CS_2 = Cartesian(origin=np.array([3, 0, 0]))
CS_3 = Cartesian(origin=np.array([6, 0, 0]))
CS_4 = Cartesian(origin=np.array([0, 3, 0]))
CS_5 = Cartesian(origin=np.array([3, 3, 0]))
CS_6 = Cartesian(origin=np.array([6, 3, 0]))
step = 1.0  # in degrees
# to visualise the coordinate system basis the module Visual is used

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure

@mlab.animate(delay=1000)
def anim():
    cs_box_1, arrows_1, labels_1 = Visual.draw_CS_box(fig, CS_1)
    cs_box_2, arrows_2, labels_2 = Visual.draw_CS_box(fig, CS_2)
    cs_box_3, arrows_3, labels_3 = Visual.draw_CS_box(fig, CS_3)
    cs_box_4, arrows_4, labels_4 = Visual.draw_CS_box(fig, CS_4)
    cs_box_5, arrows_5, labels_5 = Visual.draw_CS_box(fig, CS_5)
    cs_box_6, arrows_6, labels_6 = Visual.draw_CS_box(fig, CS_6)
    direction = 1
    while True:
        CS_1.rotate_axis_angle(np.array([0, 1, 0]), np.deg2rad(step))  # this is inplace transform
        CS_2.rotate_axis_angle(np.array([1, 0, 0]), np.deg2rad(step))  # this is inplace transform
        CS_3.rotate_axis_angle(np.array([0, 0, 1]), np.deg2rad(step))  # this is inplace transform
        CS_4.set_eulers(CS_4.euler_angles + np.array([0, 0, np.deg2rad(step)]))
        CS_5.set_eulers(CS_5.euler_angles + direction * np.array([0, np.deg2rad(step), 0]))
        CS_6.set_eulers(CS_6.euler_angles - np.array([np.deg2rad(step), 0, 0]))

        if direction == 1 and abs(np.pi - CS_5.euler_angles[1]) < np.deg2rad(step):
            direction *= -1
        elif direction == -1 and abs(CS_5.euler_angles[1]) < np.deg2rad(step):
            direction *= -1

        cs_box_1, arrows_1, labels_1 = Visual.update_CS_box(CS_1, cs_box_1, arrows_1, labels_1)
        cs_box_2, arrows_2, labels_2 = Visual.update_CS_box(CS_2, cs_box_2, arrows_2, labels_2)
        cs_box_3, arrows_3, labels_3 = Visual.update_CS_box(CS_3, cs_box_3, arrows_3, labels_3)
        cs_box_4, arrows_4, labels_4 = Visual.update_CS_box(CS_4, cs_box_4, arrows_4, labels_4)
        cs_box_5, arrows_5, labels_5 = Visual.update_CS_box(CS_5, cs_box_5, arrows_5, labels_5)
        cs_box_6, arrows_6, labels_6 = Visual.update_CS_box(CS_6, cs_box_6, arrows_6, labels_6)
        yield

anim()
mlab.show()
