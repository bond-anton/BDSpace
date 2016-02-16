import numpy as np
from matplotlib import pyplot as plt
from Space.Coordinates import Cartesian


def error(CS_1, CS_2):
    delta = CS_1.basis - CS_2.basis
    return np.sqrt(np.sum(delta**2))

# Create cartesian coordinate system

# if you don't pass arguments the basis coincide with 'Absolute' (mayavi) coordinate system
CS_1 = Cartesian(origin=np.array([0, 0, 0]))
CS_2 = Cartesian(origin=np.array([0, 0, 0]))

turns = 10
axis = np.array([1, 1, 1])
steps = 1  # per turn
max_steps_order = 4
errors = []

for order in range(max_steps_order + 1):
    steps = 10**order # per turn
    print "Processing %g steps per turn" % steps
    step = 2 * np.pi / steps
    print "  Angle increment is %g rad (%g deg)" % (step, np.rad2deg(step))
    laps_errors = []
    for i in range(turns):
        for k in range(steps):
            CS_2.rotate_axis_angle(axis, step)
        laps_errors.append(error(CS_1, CS_2))
    errors.append(np.mean(laps_errors))

print errors

plt.semilogy(np.arange(max_steps_order + 1), errors, 'r-o')
plt.show()
