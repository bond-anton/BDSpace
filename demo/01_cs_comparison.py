from __future__ import division, print_function
import numpy as np
try:
    from matplotlib import pyplot as plt
    use_mpl = True
except ImportError:
    use_mpl = True
from BDSpace.Coordinates import Cartesian


def error(cs_1, cs_2):
    delta = cs_1.basis - cs_2.basis
    return np.sqrt(np.sum(delta ** 2))

# Create cartesian coordinate system

# if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system
cs_1 = Cartesian(origin=np.array([0, 0, 0]), euler_angles_convention='Bunge')
cs_2 = Cartesian(origin=np.array([0, 0, 0]), euler_angles_convention='Bunge')

print(cs_1)

turns = 10
axis = np.array([1, 1, 1])
steps = 1  # per turn
max_steps_order = 4
errors = []

for order in range(max_steps_order + 1):
    steps = 10**order # per turn
    print("Processing %g steps per turn" % steps)
    step = 2 * np.pi / steps
    print("  Angle increment is %g rad (%g deg)" % (step, np.rad2deg(step)))
    laps_errors = []
    print("    Lap:",)
    for i in range(turns):
        print(i+1,)
        for k in range(steps):
            cs_2.rotate_axis_angle(axis, step)
        laps_errors.append(error(cs_1, cs_2))
    print("done.\n")
    errors.append(np.mean(laps_errors))

print(errors)

if use_mpl:
    plt.semilogy(np.arange(max_steps_order + 1), errors, 'r-o')
    plt.show()
