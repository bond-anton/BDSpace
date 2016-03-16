from __future__ import division
import numpy as np

from Space.Coordinates.transforms import reduce_angle
from Space.Figure import Figure


class ToricWedge(Figure):

    def __init__(self, name='Toric wedge', coordinate_system=None,
                 phi=np.pi/2, theta=np.array([0.0, np.pi/2]),
                 r_torus=1.0, r_tube=np.array([0, 0.25])):
        self.r_torus = float(r_torus)
        self.r_tube = np.array(r_tube, dtype=np.float)
        reduced_theta = reduce_angle(theta, keep_sign=True)
        theta_min = min(reduced_theta)
        theta_max = max(reduced_theta)
        if theta_max - theta_min >= 2*np.pi:
            theta_min = 0.0
            theta_max = 2 * np.pi
        self.theta = np.array([theta_min, theta_max])
        self.phi = reduce_angle(phi)

        super(ToricWedge, self).__init__(name, coordinate_system=coordinate_system)
