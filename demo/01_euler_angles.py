#!/bin/env python

import numpy as np
from mayavi import mlab

from Space.Coordinates import Cartesian
from Space import Visual

# Euler's angles are used in the proper notation Z (phi1) - X' (Phi) - Z" (phi2)
# phi1 and phi2 are defined to have modulo 2*pi radians [-pi; pi] or [0; 2*pi]
# Phi is defined to have modulo pi radians [-pi/2; pi/2] or [0; pi]

M = 10
N = 5

# Here we use degrees and numpy deg2rad for conversion

scale = min(360 / (2 * M), 180 / (2 * N)) # scale factor for mayavi scene

fig = mlab.figure('CS demo', bgcolor=(0, 0, 0))  # Create the mayavi figure

for phi1 in np.linspace(0, 360, M, endpoint=True):
    for Phi in np.linspace(0, 180, N, endpoint=True):
        for phi2 in np.linspace(0, 360, M, endpoint=True):
            euler_angles = np.deg2rad(np.array([phi1, Phi, phi2]))
            # Create cartesian coordinate system
            CS = Cartesian(origin=np.array([phi1, Phi, phi2]), labels=['i1', 'i2', 'i3'])
            # Set CS orientation using Euler's angles
            CS.set_eulers(euler_angles)
            # CS_box visualize CS as a cube colored according to Euler's angles
            Visual.draw_CS_box(fig, CS, scale=scale, draw_axes=True)
# mlab.outline(extent=[0, 360, 0, 180, 0, 360])  # uncomment to draw white outline
mlab.show()
