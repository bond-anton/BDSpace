from __future__ import division, print_function
import numpy as np
from Space.Coordinates import Cartesian

coordinate_system = Cartesian(origin=np.array([0, 0, 0]), euler_angles_convention='canova')
print(coordinate_system)
